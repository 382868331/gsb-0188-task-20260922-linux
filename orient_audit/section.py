"""平面剖切：有理平面与四面体网格相交的精确截面。

在 :mod:`orient_audit` 的定向审计能力之上新增 :func:`section`：

- 平面 ``n·x = d`` 的法向与偏置全部为有理数，交点坐标用
  :class:`fractions.Fraction` 精确给出（拒绝 float）；
- **顶点统一分类**：每个原顶点相对平面的符号（正/负/在平面上）全局只
  计算一次；
- **共享边交点以原顶点对为身份**：交点 id 为 ``("edge", (i, j))``，
  相邻单元在共享棱上得到的是同一个点对象，不靠浮点坐标焊接；
  在平面上的原顶点 id 为 ``("vertex", i)``；
- 每个非退化单元给出三角形或四边形片；四边形按全局交点 ID 选择一条
  确定的对角线拆成两个三角形；
- 完全共面的原面全局只保留一份，并记录所有来源单元；
- 仅接触平面的点/棱作为零面积结果单独列出，不混入面片；
- 所有面片朝向由给定法向统一决定（叉积·法向 > 0）；负向输入单元先按
  原修正规则（只交换最后两个顶点）处理，零体积单元拒绝剖切。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .audit import repair_orientation
from .errors import MeshInputError, ZeroVolumeError
from .exact import parse_coord
from .mesh import FACE_LOCAL, TetMesh

# id 种类标签（元组首元素），同时作为可读身份的一部分。
_VERTEX = "vertex"
_EDGE = "edge"

PointId = tuple  # ("vertex", i) 或 ("edge", (i, j))
Coord = tuple[Fraction, Fraction, Fraction]


# ---------------------------------------------------------------- 平面输入


@dataclass(frozen=True)
class Plane:
    """有理平面 ``n·x = offset``（不可变）。

    属性:
        normal: 三元组 :class:`Fraction`，构造时保证非零向量。
        offset: :class:`Fraction`。
    """

    normal: tuple[Fraction, Fraction, Fraction]
    offset: Fraction

    @classmethod
    def build(cls, normal, offset=0) -> "Plane":
        """解析并校验平面。

        - ``normal`` 必须是 3 个有理分量（int/Fraction/有理字符串），
          拒绝 float、bool，错误可定位（``normal[k]``）；
        - 零法向向量以 :class:`MeshInputError` 拒绝（位置 ``normal``）；
        - ``offset`` 为有理标量（位置 ``offset``）。
        """

        if isinstance(normal, str) or not isinstance(normal, (list, tuple)):
            raise MeshInputError("normal 必须是长度为 3 的有理序列", "normal")
        if len(normal) != 3:
            raise MeshInputError(
                f"normal 必须恰好有 3 个分量，得到 {len(normal)} 个", "normal"
            )
        n = tuple(parse_coord(c, f"normal[{k}]") for k, c in enumerate(normal))
        d = parse_coord(offset, "offset")
        if n == (Fraction(0), Fraction(0), Fraction(0)):
            raise MeshInputError("平面法向不能是零向量", "normal")
        return cls(n, d)  # type: ignore[arg-type]


# ---------------------------------------------------------------- 结果类型


@dataclass(frozen=True)
class SectionPoint:
    """截面上的一个点：全局身份 + 精确有理坐标。

    属性:
        id: ``("vertex", i)`` 表示平面上的原顶点 i；
            ``("edge", (i, j))``（i<j）表示原棱 (i,j) 与平面的交点。
        coord: 精确有理坐标。
    """

    id: PointId
    coord: Coord

    @property
    def kind(self) -> str:
        return self.id[0]


@dataclass(frozen=True)
class SectionPolygon:
    """单个非退化单元与平面相交的面积片（三角形或四边形）。

    属性:
        tet_index: 来源单元下标（面片与单元一一对应；共面原面除外）。
        vertices: 按绕序排列的全局点 id（3 或 4 个），朝向已按给定法向
                  统一为叉积·法向 > 0。
        triangles: 三角形剖分（三角形片 1 个，四边形片 2 个）；
                   四边形对角线由全局点 ID 确定性选择。
        area_coefficient: 面积系数 α（Fraction）：面积 = α·‖n‖，
                          α = n·(q−p)×(r−p)/(2 n·n)，朝向修正后恒正。
    """

    tet_index: int
    vertices: tuple[PointId, ...]
    triangles: tuple[tuple[PointId, PointId, PointId], ...]
    area_coefficient: Fraction

    @property
    def is_quad(self) -> bool:
        return len(self.vertices) == 4


@dataclass(frozen=True)
class CoplanarFacePatch:
    """与平面完全共面的原三角面（全局唯一一份）。

    属性:
        key: 排序后的三个原顶点下标（无序三点键）。
        vertices: 按给定法向定向的三个 ``("vertex", ...)`` 点 id。
        triangles: 恒为单个三角形。
        tet_indices: 引用该面的全部单元下标（升序）。
        area_coefficient: 面积系数 α（恒正）。
    """

    key: tuple[int, int, int]
    vertices: tuple[PointId, PointId, PointId]
    triangles: tuple[tuple[PointId, PointId, PointId], ...]
    tet_indices: tuple[int, ...]
    area_coefficient: Fraction


@dataclass(frozen=True)
class ContactPoint:
    """仅以单个顶点接触平面的零面积结果。

    属性:
        id: 接触点全局 id（必为平面上的原顶点）。
        tet_indices: 只在该点接触（其余顶点严格同侧）的全部单元。
    """

    id: PointId
    tet_indices: tuple[int, ...]


@dataclass(frozen=True)
class ContactEdge:
    """仅以一条棱接触平面的零面积结果。

    属性:
        endpoints: 接触棱两端点 id（按全局 ID 排序，均为平面上的原顶点）。
        tet_indices: 只在该棱接触（另外两顶点严格同侧）的全部单元。
    """

    endpoints: tuple[PointId, PointId]
    tet_indices: tuple[int, ...]


@dataclass(frozen=True)
class SectionReport:
    """一次平面剖切的完整结果（不可变）。

    属性:
        plane: 输入平面。
        points: 本截面用到的全部全局点（含面片顶点与接触点），按全局 ID
                排序；相同身份只出现一次。
        polygons: 非退化单元的三角/四边形面片。
        coplanar_faces: 完全共面原面（每个三点键一份，记录全部来源单元）。
        contact_points: 零面积点接触。
        contact_edges: 零面积棱接触。
        total_area_coefficient: 全部面片（含共面原面）的面积系数和。
    """

    plane: Plane
    points: tuple[SectionPoint, ...]
    polygons: tuple[SectionPolygon, ...]
    coplanar_faces: tuple[CoplanarFacePatch, ...]
    contact_points: tuple[ContactPoint, ...]
    contact_edges: tuple[ContactEdge, ...]
    total_area_coefficient: Fraction

    def point_coord(self, pid: PointId) -> Coord:
        """按全局身份查精确坐标（身份不存在则 KeyError）。"""

        return self._coord_map[pid]

    @property
    def _coord_map(self) -> dict[PointId, Coord]:
        return {p.id: p.coord for p in self.points}

    def total_area_float(self) -> float:
        """总截面积的一次 double 估算：``float(α)·sqrt(float(n·n))``。

        唯一的浮点出口；精确判定从不使用它。误差为若干 ULP 量级
        （一次 Fraction→double 舍入 + 一次 sqrt 舍入），见 README。
        """

        import math

        n = self.plane.normal
        nn = float(n[0] * n[0] + n[1] * n[1] + n[2] * n[2])
        return float(self.total_area_coefficient) * math.sqrt(nn)


# ---------------------------------------------------------------- 主入口


def section(mesh: TetMesh, plane: Plane) -> SectionReport:
    """用有理平面剖切四面体网格，纯只读，不修改、不分裂原网格。

    负向单元先经 :func:`repair_orientation` 的同一规则处理（只交换最后
    两个顶点）；存在零体积单元时抛 :class:`ZeroVolumeError`，不产生部分
    结果。
    """

    if not isinstance(mesh, TetMesh):
        raise TypeError("section 需要 TetMesh 实例（请先用 TetMesh.build 构造）")
    if not isinstance(plane, Plane):
        raise TypeError("section 需要 Plane 实例（请先用 Plane.build 构造）")

    # 原修正规则：负向单元交换最后两顶点；零体积单元拒绝（抛 ZeroVolumeError）。
    work = repair_orientation(mesh)
    assert all(a is b for a, b in zip(work.vertices, mesh.vertices, strict=True))

    n = plane.normal
    d = plane.offset
    nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2]
    verts = work.vertices

    # ---- 全局：每个原顶点相对平面的符号只算一次 ----
    signed = [
        (n[0] * v[0] + n[1] * v[1] + n[2] * v[2] - d) for v in verts
    ]

    def sign(x: Fraction) -> int:
        if x > 0:
            return 1
        if x < 0:
            return -1
        return 0

    signs = [sign(s) for s in signed]

    # ---- 全局点表：身份 -> 坐标（交点以原顶点对为身份）----
    coords: dict[PointId, Coord] = {}

    def vertex_point(i: int) -> PointId:
        pid = (_VERTEX, i)
        coords.setdefault(pid, verts[i])
        return pid

    def edge_point(i: int, j: int) -> PointId:
        a, b = (i, j) if i < j else (j, i)
        pid = (_EDGE, (a, b))
        if pid not in coords:
            sa, sb = signed[a], signed[b]
            # 精确凸插值：p = (sb·v_a − sa·v_b)/(sb − sa)。
            den = sb - sa
            coords[pid] = tuple(
                (sb * verts[a][k] - sa * verts[b][k]) / den for k in range(3)
            )
        return pid

    polygons: list[SectionPolygon] = []
    coplanar_sources: dict[tuple[int, int, int], list[int]] = {}
    contact_pt_src: dict[PointId, list[int]] = {}
    contact_ed_src: dict[tuple[PointId, PointId], list[int]] = {}

    for ti, tet in enumerate(work.tets):
        sg = [signs[i] for i in tet]
        nzero = sg.count(0)
        npos = sg.count(1)
        nneg = sg.count(-1)

        # ---- 完全共面的原面：登记到全局去重表 ----
        if nzero == 3:
            face = tuple(sorted(tet[k] for k in range(4) if sg[k] == 0))
            coplanar_sources.setdefault(face, []).append(ti)
            continue

        # ---- 其余情形：由四个面与平面的交线段构造截面边图 ----
        edge_set: set[frozenset[PointId]] = set()
        incident: set[PointId] = set()
        for local in FACE_LOCAL:
            fidx = (tet[local[0]], tet[local[1]], tet[local[2]])
            fpts: list[PointId] = []
            seen_local: set[PointId] = set()
            for vi in fidx:
                if signs[vi] == 0:
                    pid = vertex_point(vi)
                    if pid not in seen_local:
                        fpts.append(pid)
                        seen_local.add(pid)
            for a, b in ((fidx[0], fidx[1]), (fidx[1], fidx[2]), (fidx[2], fidx[0])):
                if signs[a] * signs[b] < 0:
                    pid = edge_point(a, b)
                    if pid not in seen_local:
                        fpts.append(pid)
                        seen_local.add(pid)
            for pid in fpts:
                incident.add(pid)
            if len(fpts) == 2:
                edge_set.add(frozenset(fpts))

        # ---- 零面积接触：退化成一点或一条棱 ----
        if not incident:
            # 四个顶点严格同侧：该单元与平面完全不相交。
            continue
        if len(incident) == 1 and not edge_set:
            pid = next(iter(incident))
            contact_pt_src.setdefault(pid, []).append(ti)
            continue
        if len(incident) == 2 and len(edge_set) == 1:
            p1, p2 = sorted(incident, key=id_key)
            contact_ed_src.setdefault((p1, p2), []).append(ti)
            continue

        # ---- 面积片：边图必须是一个 3 或 4 点的环 ----
        adj: dict[PointId, set[PointId]] = {pid: set() for pid in incident}
        for e in edge_set:
            a, b = tuple(e)
            adj[a].add(b)
            adj[b].add(a)
        cycle = _walk_cycle(adj)
        cycle = _orient_cycle(cycle, n, coords)
        tris = _triangulate(cycle)
        alpha = sum(_triangle_alpha(t, n, nn, coords) for t in tris)
        polygons.append(
            SectionPolygon(
                tet_index=ti,
                vertices=tuple(cycle),
                triangles=tuple(tris),
                area_coefficient=alpha,
            )
        )

    # ---- 共面原面：每个三点键一份，按法向统一朝向 ----
    coplanar_patches: list[CoplanarFacePatch] = []
    for key in sorted(coplanar_sources):
        cyc = [vertex_point(i) for i in key]
        cyc = _orient_cycle(cyc, n, coords)
        tris = _triangulate(cyc)
        alpha = _triangle_alpha(tris[0], n, nn, coords)
        coplanar_patches.append(
            CoplanarFacePatch(
                key=key,
                vertices=tuple(cyc),
                triangles=tuple(tris),
                tet_indices=tuple(sorted(coplanar_sources[key])),
                area_coefficient=alpha,
            )
        )

    points = tuple(
        SectionPoint(pid, coords[pid]) for pid in sorted(coords, key=id_key)
    )
    contact_points = tuple(
        ContactPoint(pid, tuple(sorted(src)))
        for pid, src in sorted(contact_pt_src.items(), key=lambda kv: id_key(kv[0]))
    )
    contact_edges = tuple(
        ContactEdge(key, tuple(sorted(src)))
        for key, src in sorted(
            contact_ed_src.items(), key=lambda kv: tuple(id_key(x) for x in kv[0])
        )
    )
    total_alpha = sum(
        [poly.area_coefficient for poly in polygons]
        + [cp.area_coefficient for cp in coplanar_patches],
        Fraction(0),
    )

    return SectionReport(
        plane=plane,
        points=points,
        polygons=tuple(polygons),
        coplanar_faces=tuple(coplanar_patches),
        contact_points=contact_points,
        contact_edges=contact_edges,
        total_area_coefficient=total_alpha,
    )


# ---------------------------------------------------------------- 几何工具


def id_key(pid: PointId) -> tuple[int, int, int]:
    """全局点身份的全序：原顶点优先，其次按棱顶点对。"""

    if pid[0] == _VERTEX:
        return (0, pid[1], 0)
    return (1, pid[1][0], pid[1][1])


def _walk_cycle(
    adj: dict[PointId, set[PointId]],
) -> list[PointId]:
    """从边图按确定性顺序走出环（每点度数必须为 2）。"""

    start = min(adj, key=id_key)
    cycle = [start]
    prev = start
    cur = min(adj[start], key=id_key)
    while cur != start:
        cycle.append(cur)
        nxt_candidates = [x for x in adj[cur] if x != prev]
        if len(nxt_candidates) != 1:
            raise AssertionError(  # 组合数学上不应发生
                f"截面边图不是环：单元交点邻接异常于 {cur}: {adj[cur]}"
            )
        prev, cur = cur, nxt_candidates[0]
    if len(cycle) != len(adj) or len(cycle) not in (3, 4):
        raise AssertionError(
            f"截面交点环大小非法：{len(cycle)} 个顶点（应为 3 或 4）"
        )
    return cycle


def _cross(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def _dot(u, v) -> Fraction:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _sub(u, v):
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def _orient_cycle(
    cycle: list[PointId], n, coords: dict[PointId, Coord]
) -> list[PointId]:
    """把绕序统一为 n·(b−a)×(c−a) > 0（点全在平面上，必不为零）。"""

    p0, p1, p2 = coords[cycle[0]], coords[cycle[1]], coords[cycle[2]]
    oriented_sign = _dot(n, _cross(_sub(p1, p0), _sub(p2, p0)))
    if oriented_sign == 0:
        raise AssertionError("截面片面积为零却进入面积片分支（组合分类错误）")
    if oriented_sign < 0:
        return [cycle[0], *reversed(cycle[1:])]
    return cycle


def _triangulate(
    cycle: list[PointId],
) -> list[tuple[PointId, PointId, PointId]]:
    """三角形直接返回；四边形以全局 ID 最小顶点与其对点的对角线拆分。"""

    if len(cycle) == 3:
        return [(cycle[0], cycle[1], cycle[2])]
    # 确定性对角线：全局 ID 最小的顶点 -> 环上的对点。
    k = min(range(len(cycle)), key=lambda i: id_key(cycle[i]))
    a, b, c, d = (cycle[k], cycle[(k + 1) % 4], cycle[(k + 2) % 4],
                  cycle[(k + 3) % 4])
    return [(a, b, c), (a, c, d)]


def _triangle_alpha(tri, n, nn: Fraction, coords: dict[PointId, Coord]) -> Fraction:
    """面积系数 α = n·(q−p)×(r−p)/(2 n·n)，定向正确时恒正。"""

    p, q, r = coords[tri[0]], coords[tri[1]], coords[tri[2]]
    return _dot(n, _cross(_sub(q, p), _sub(r, p))) / (2 * nn)
