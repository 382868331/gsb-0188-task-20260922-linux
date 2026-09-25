"""平面截面：有理平面与四面体网格的精确有理相交。

契约摘要（详见 README）：

- 平面由有理法向 ``normal`` 与偏置 ``offset`` 给出：``normal·x = offset``；
- 每个顶点按精确符号 ``sign(normal·v - offset)`` 统一分类
  （正 / 负 / 恰好在平面上），全程 :class:`fractions.Fraction`，无容差；
- 截面点身份由来源决定：原顶点给 ``("vertex", i)``，棱交点给
  ``("edge", i, j)``（``i < j`` 的原顶点对）。共享棱的交点天然同身份，
  不做任何浮点坐标近似焊接；
- 每个非退化四面体贡献三角形或四边形截面；四边形按**全局交点 ID**
  选择确定对角线（含最小全局 ID 的那条）剖成两个三角形；
- 完全共面的原面只保留一份（``kind == "coplanar_face"``），
  ``source_tets`` 记录所有来源单元；
- 只接触点 / 棱的单元列入零面积接触（``point_contacts`` /
  ``edge_contacts``），不混进面片；
- 所有面片朝向统一与给定法向一致（片法向与 ``normal`` 点积为正）；
- 负向输入单元先按原修正规则（只交换最后两个顶点，与
  :func:`orient_audit.audit.repair_orientation` 相同）在内部修正后再截面；
  零体积（退化）单元不参与截面，列入 ``degenerate_tets``；
- 输入网格不被修改，也不做任何边分裂（截面点是新建的有理点，
  原顶点对象恒等保留）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .errors import MeshInputError
from .exact import parse_coord, signed_volume6
from .mesh import TetMesh, Vertex

#: 截面片类别：普通三角截面 / 四边形截面（已剖两三角）/ 完全共面的原面。
PATCH_KINDS = ("triangle", "quad", "coplanar_face")


def _sign(x: Fraction) -> int:
    return (x > 0) - (x < 0)


def _sub(a: Vertex, b: Vertex) -> Vertex:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(u: Vertex, v: Vertex) -> Vertex:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


@dataclass(frozen=True)
class SectionPoint:
    """截面点表中的一个点（精确有理坐标）。

    属性:
        pid: 全局交点 ID（在本报告的 ``points`` 中的下标）。
        coords: 精确 :class:`Fraction` 三元组。
        origin: 点身份——``("vertex", i)`` 表示就是原网格顶点 ``i``；
                ``("edge", i, j)`` 表示原网格棱 ``(i, j)``（``i < j``）
                与平面的交点。身份只依赖原顶点对，与坐标数值无关。
    """

    pid: int
    coords: Vertex
    origin: tuple


@dataclass(frozen=True)
class SectionPatch:
    """一张面积大于零的有理多边形片。

    属性:
        kind: ``"triangle"`` / ``"quad"`` / ``"coplanar_face"``。
        polygon: 有序边界点的全局 ID（朝向与平面法向一致）。
        triangles: 三角剖分（全局点 ID 三元组）；四边形给出两个三角形，
                   对角线为含最小全局交点 ID 的那一条。
        source_tets: 来源单元下标（共面面可多于一个，升序）。
        area2: 面积的精确平方（``Fraction``）。有理点的多边形面积
               一般不是有理数，故精确出口给平方；``area()`` 是浮点便利。
    """

    kind: str
    polygon: tuple[int, ...]
    triangles: tuple[tuple[int, int, int], ...]
    source_tets: tuple[int, ...]
    area2: Fraction

    def area(self) -> float:
        """面积的一次浮点求值（对精确 ``area2`` 开平方，非精确出口）。"""

        return math.sqrt(self.area2)


@dataclass(frozen=True)
class PointContact:
    """零面积点接触：某单元只在原顶点处碰到平面。"""

    vertex: int
    point_id: int
    source_tets: tuple[int, ...]


@dataclass(frozen=True)
class EdgeContact:
    """零面积棱接触：某单元只沿原网格的一条棱贴在平面上。"""

    edge: tuple[int, int]
    point_ids: tuple[int, int]
    source_tets: tuple[int, ...]


@dataclass(frozen=True)
class SectionReport:
    """一次平面截面的完整结果（不可变）。

    属性:
        normal: 平面法向（精确 Fraction 三元组）。
        offset: 平面偏置，平面方程为 ``normal·x = offset``。
        points: 全部截面点（按全局 ID 排序）。
        patches: 面积大于零的多边形片。
        point_contacts / edge_contacts: 零面积接触，与面片分开列出。
        degenerate_tets: 零体积单元下标（不参与截面）。
        repaired_tets: 截面前按原规则（交换最后两顶点）内部修正的
                       负向单元下标；输入网格本身未被修改。
    """

    normal: Vertex
    offset: Fraction
    points: tuple[SectionPoint, ...]
    patches: tuple[SectionPatch, ...]
    point_contacts: tuple[PointContact, ...]
    edge_contacts: tuple[EdgeContact, ...]
    degenerate_tets: tuple[int, ...]
    repaired_tets: tuple[int, ...]

    def total_area(self) -> float:
        """全部面片面积之和的浮点值（逐片对精确 area2 开方后求和）。"""

        return math.fsum(p.area() for p in self.patches)

    def plane_residuals(self) -> tuple[Fraction, ...]:
        """每个截面点的精确平面残差 ``normal·p - offset``（应全为 0）。"""

        return tuple(
            sum(self.normal[k] * p.coords[k] for k in range(3)) - self.offset
            for p in self.points
        )


def _parse_plane(normal: Any, offset: Any) -> tuple[Vertex, Fraction]:
    """解析并校验平面参数；所有错误带 ``plane.*`` 定位。"""

    if isinstance(normal, str) or not isinstance(normal, (list, tuple)):
        raise MeshInputError("平面法向必须是长度为 3 的有理序列", "plane.normal")
    if len(normal) != 3:
        raise MeshInputError(
            f"平面法向必须恰好 3 个分量，得到 {len(normal)} 个", "plane.normal"
        )
    n = tuple(
        parse_coord(c, f"plane.normal[{k}]") for k, c in enumerate(normal)
    )
    if all(c == 0 for c in n):
        raise MeshInputError("平面法向不能是零向量", "plane.normal")
    d = parse_coord(offset, "plane.offset")
    return n, d  # type: ignore[return-value]


def section(mesh: TetMesh, normal: Any, offset: Any) -> SectionReport:
    """计算四面体网格与有理平面的截面，纯只读，不修改输入。

    参数:
        mesh: 已校验的 :class:`TetMesh`。
        normal: 平面法向，3 个有理分量（int / Fraction / 有理字符串，
                拒绝 float 与零向量）。
        offset: 平面偏置（同上有理类型），平面为 ``normal·x = offset``。

    返回:
        :class:`SectionReport`；语义见模块 docstring 与 README。
    """

    if not isinstance(mesh, TetMesh):
        raise TypeError("section 需要 TetMesh 实例（请先用 TetMesh.build 构造）")
    n, d = _parse_plane(normal, offset)

    # 每个顶点的精确带符号值 s(v) = n·v - d，全网格统一分类一次。
    svals = [sum(n[k] * v[k] for k in range(3)) - d for v in mesh.vertices]

    # 点身份表：("vertex", i) 或 ("edge", i, j) -> 精确坐标。
    point_coords: dict[tuple, Vertex] = {}

    def vertex_point(i: int) -> tuple:
        key = ("vertex", i)
        point_coords.setdefault(key, mesh.vertices[i])
        return key

    def edge_point(i: int, j: int) -> tuple:
        # 身份 = 原顶点对（升序）；坐标用规范顺序计算，保证同一棱
        # 无论从哪个单元进入都得到完全相同的 Fraction。
        a, b = (i, j) if i < j else (j, i)
        key = ("edge", a, b)
        if key not in point_coords:
            sa, sb = svals[a], svals[b]
            t = sa / (sa - sb)  # 两端异号 => 分母非零，t ∈ (0,1) 精确
            va, vb = mesh.vertices[a], mesh.vertices[b]
            point_coords[key] = tuple(
                va[k] + t * (vb[k] - va[k]) for k in range(3)
            )
        return key

    # 中间产物（仍以身份键表示，最后统一换成全局 ID）。
    raw_patches: list[tuple[str, tuple[tuple, ...], int]] = []
    coplanar: dict[tuple[int, int, int], list[int]] = {}
    point_contacts_raw: dict[int, list[int]] = {}
    edge_contacts_raw: dict[tuple[int, int], list[int]] = {}
    degenerate: list[int] = []
    repaired: list[int] = []

    for ti, tet in enumerate(mesh.tets):
        pts = tuple(mesh.vertices[k] for k in tet)
        v6 = signed_volume6(pts[0], pts[1], pts[2], pts[3])
        if v6 == 0:
            degenerate.append(ti)  # 退化单元不参与截面
            continue
        if v6 < 0:
            # 与 repair_orientation 相同的规则：只交换最后两个顶点。
            tet = (tet[0], tet[1], tet[3], tet[2])
            repaired.append(ti)

        signs = [_sign(svals[k]) for k in tet]
        zeros = [m for m in range(4) if signs[m] == 0]
        pos = [m for m in range(4) if signs[m] > 0]
        neg = [m for m in range(4) if signs[m] < 0]
        z = len(zeros)

        if z == 0:
            if not pos or not neg:
                continue  # 严格同侧，无交
            if len(pos) == 1 or len(pos) == 3:
                # 1-3 分布：孤立顶点与另外三顶点的三条棱各交一点。
                solo = pos[0] if len(pos) == 1 else neg[0]
                others = neg if len(pos) == 1 else pos
                poly = tuple(edge_point(tet[solo], tet[o]) for o in others)
                raw_patches.append(("triangle", poly, ti))
            else:
                # 2-2 分布：四条跨侧棱给出四边形，环序由面结构确定。
                up = [tet[m] for m in pos]
                dn = [tet[m] for m in neg]
                poly = (
                    edge_point(up[0], dn[0]),
                    edge_point(up[1], dn[0]),
                    edge_point(up[1], dn[1]),
                    edge_point(up[0], dn[1]),
                )
                raw_patches.append(("quad", poly, ti))
        elif z == 1:
            w = zeros[0]
            if not pos or not neg:
                # 其余三点严格同侧：仅在原顶点处接触。
                point_contacts_raw.setdefault(tet[w], []).append(ti)
                vertex_point(tet[w])
            else:
                # 过原顶点且切两条棱：三角形。
                solo = pos[0] if len(pos) == 1 else neg[0]
                others = neg if len(pos) == 1 else pos
                poly = (
                    vertex_point(tet[w]),
                    edge_point(tet[solo], tet[others[0]]),
                    edge_point(tet[solo], tet[others[1]]),
                )
                raw_patches.append(("triangle", poly, ti))
        elif z == 2:
            w0, w1 = zeros
            if not pos or not neg:
                # 一条原棱贴在平面上，另两点严格同侧：棱接触。
                a, b = sorted((tet[w0], tet[w1]))
                edge_contacts_raw.setdefault((a, b), []).append(ti)
                vertex_point(a)
                vertex_point(b)
            else:
                # 两个在平面上 + 一正一负：三角形（面上棱 + 一个棱交点）。
                poly = (
                    vertex_point(tet[w0]),
                    vertex_point(tet[w1]),
                    edge_point(tet[pos[0]], tet[neg[0]]),
                )
                raw_patches.append(("triangle", poly, ti))
        elif z == 3:
            # 一个原面完全共面：按键去重，记录所有来源单元。
            face = tuple(sorted(tet[m] for m in zeros))
            coplanar.setdefault(face, []).append(ti)
            for m in zeros:
                vertex_point(tet[m])
        else:  # z == 4 意味着零体积，已在上面拦截
            raise AssertionError("非退化四面体不可能四顶点共面")

    # ---- 全局交点 ID：按身份键规范排序，与遍历顺序无关 ----
    keys_sorted = sorted(
        point_coords, key=lambda k: (k[0] != "vertex", k[1:])
    )
    pid_of = {key: pid for pid, key in enumerate(keys_sorted)}
    points = tuple(
        SectionPoint(pid, point_coords[key], key)
        for pid, key in enumerate(keys_sorted)
    )

    def orient(poly_keys: tuple[tuple, ...]) -> tuple[tuple, ...]:
        """统一朝向：片法向与给定法向点积必须为正（精确判定）。"""

        c = [point_coords[k] for k in poly_keys]
        normal_vec = _cross(_sub(c[1], c[0]), _sub(c[2], c[0]))
        dot = sum(normal_vec[k] * n[k] for k in range(3))
        # 非退化截面片的法向必与平面法向平行，dot 不可能为 0。
        assert dot != 0
        return poly_keys if dot > 0 else tuple(reversed(poly_keys))

    def triangulate(poly_pids: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
        if len(poly_pids) == 3:
            return (poly_pids,)
        # 四边形：选含最小全局交点 ID 的对角线，规则全局确定。
        q = poly_pids
        m = min(q)
        if q[0] == m or q[2] == m:
            return ((q[0], q[1], q[2]), (q[0], q[2], q[3]))
        return ((q[1], q[2], q[3]), (q[1], q[3], q[0]))

    def area2_of(triangles: tuple[tuple[int, int, int], ...]) -> Fraction:
        # 同向共面三角形的叉积和 = 2 倍面积向量，取模方 / 4，精确。
        acc = [Fraction(0), Fraction(0), Fraction(0)]
        for tri in triangles:
            c = [points[p].coords for p in tri]
            cr = _cross(_sub(c[1], c[0]), _sub(c[2], c[0]))
            for k in range(3):
                acc[k] += cr[k]
        return sum(acc[k] * acc[k] for k in range(3)) / 4

    def make_patch(
        kind: str, poly_keys: tuple[tuple, ...], source_tets: tuple[int, ...]
    ) -> SectionPatch:
        poly_pids = tuple(pid_of[k] for k in orient(poly_keys))
        triangles = triangulate(poly_pids)
        return SectionPatch(
            kind, poly_pids, triangles, source_tets, area2_of(triangles)
        )

    patches: list[SectionPatch] = [
        make_patch(kind, poly, (ti,)) for kind, poly, ti in raw_patches
    ]
    for face, tis in sorted(coplanar.items()):
        poly = tuple(("vertex", i) for i in face)
        patches.append(
            make_patch("coplanar_face", poly, tuple(sorted(tis)))
        )

    point_contacts = tuple(
        PointContact(v, pid_of[("vertex", v)], tuple(sorted(tis)))
        for v, tis in sorted(point_contacts_raw.items())
    )
    edge_contacts = tuple(
        EdgeContact(
            e,
            (pid_of[("vertex", e[0])], pid_of[("vertex", e[1])]),
            tuple(sorted(tis)),
        )
        for e, tis in sorted(edge_contacts_raw.items())
    )

    return SectionReport(
        normal=n,
        offset=d,
        points=points,
        patches=tuple(patches),
        point_contacts=point_contacts,
        edge_contacts=edge_contacts,
        degenerate_tets=tuple(degenerate),
        repaired_tets=tuple(repaired),
    )
