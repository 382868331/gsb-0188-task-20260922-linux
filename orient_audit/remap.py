"""守恒重映射：把源四面体网格上的单元常量密度精确转到目标网格。

权重来自源/目标四面体对的真实交叠体积：对每个源-目标四面体对，
用目标四面体的四个半空间精确裁剪源四面体（凸多面体裁剪，全部走
fractions.Fraction），得到交叠凸多面体并精确求体积。

只处理四面体（凸）单元；不实现一般非凸单元。输入网格对象只读。

输出：

- 稀疏传递权重 weights：(目标下标, 源下标, 交叠体积)，仅含正体积项；
- 每个目标单元的覆盖体积、未覆盖体积（缺口保留，不归一）；
- 每个目标单元的质量与覆盖部分平均密度；
- 总量误差：目标总质量减去源网格被覆盖部分的质量（精确为 0 才守恒）。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .errors import MeshInputError
from .exact import parse_coord, signed_volume6
from .mesh import FACE_LOCAL, TetMesh, Vertex

# 一个凸多面体：顶点表 + 面表（面为顶点下标环，朝向一致朝外）。
_Polytope = tuple[list[Vertex], list[list[int]]]


def _dot3(a: Vertex, b: Vertex) -> Fraction:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub3(a: Vertex, b: Vertex) -> Vertex:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross3(a: Vertex, b: Vertex) -> Vertex:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _tet_polytope(vertices, tet) -> _Polytope:
    """把四面体转成朝外定向的凸多面体表示。"""

    pts = [vertices[i] for i in tet]
    faces: list[list[int]] = []
    for local in FACE_LOCAL:
        a, b, c = (pts[local[0]], pts[local[1]], pts[local[2]])
        # 面法向应指向单元外：与 (对面顶点 - a) 点积 < 0。
        opposite = next(k for k in range(4) if k not in local)
        outward = _dot3(_cross3(_sub3(b, a), _sub3(c, a)), _sub3(pts[opposite], a)) < 0
        faces.append([local[0], local[1], local[2]] if outward else [local[0], local[2], local[1]])
    return pts, faces


def _clip_by_halfspace(poly: _Polytope, normal: Vertex, offset: Fraction) -> _Polytope:
    """用半空间 normal·x <= offset 裁剪凸多面体（精确有理运算）。

    返回裁剪后的多面体；若被裁空则顶点表为空。
    """

    verts, faces = poly
    values = [_dot3(normal, v) - offset for v in verts]
    if all(v > 0 for v in values):
        return [], []
    new_verts: list[Vertex] = list(verts)
    new_faces: list[list[int]] = []
    cap_edges: list[tuple[int, int]] = []
    n_old = len(verts)
    cache: dict[Vertex, int] = {}
    on_plane: dict[int, bool] = {i: v == 0 for i, v in enumerate(values)}
    for face in faces:
        poly_face: list[int] = []
        m = len(face)
        for i in range(m):
            i0, i1 = face[i], face[(i + 1) % m]
            d0, d1 = values[i0], values[i1]
            if d0 <= 0:
                poly_face.append(i0)
            # 仅统计严格穿越（两端异号且都不在平面上）；
            # 恰好落在平面上的顶点保留，不产生新点。
            if (d0 < 0 < d1) or (d1 < 0 < d0):
                t = d0 / (d0 - d1)
                p = tuple(
                    verts[i0][k] + t * (verts[i1][k] - verts[i0][k]) for k in range(3)
                )
                # 同一条棱被相邻两个面共享，交点按坐标去重，
                # 否则截面线段无法按端点链接成环。
                idx = cache.get(p)
                if idx is None:
                    idx = len(new_verts)
                    new_verts.append(p)
                    cache[p] = idx
                on_plane[idx] = True
                poly_face.append(idx)
        if len(poly_face) >= 3:
            new_faces.append(poly_face)
        # 被保留多边形中两端都落在平面上的边是截面与该面的交线。
        # 截面边若是新产生的（穿过面内部）只出现一次；若落在多面体
        # 已有棱上（含被裁成线段的退化面）会被两个相邻面各贡献一次，
        # 方向相反——按无向边去重即可。
        m = len(poly_face)
        if m == 2:
            a, b = poly_face
            if on_plane[a] and on_plane[b]:
                cap_edges.append((a, b))
        elif m >= 3:
            for i in range(m):
                a, b = poly_face[i], poly_face[(i + 1) % m]
                if on_plane[a] and on_plane[b]:
                    cap_edges.append((a, b))
    # 只有真正裁掉了顶点（存在严格外侧点）才需要补截面；
    # 全部顶点都在半空间内（含恰好在平面上）时多面体原样封闭。
    if cap_edges and any(v > 0 for v in values):
        undirected = {tuple(sorted(e)) for e in cap_edges}
        adj: dict[int, list[int]] = {}
        for a, b in undirected:
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
        # 截面是单个凸多边形：沿无向边链接成环。
        start = next(iter(adj))
        loop = [start]
        prev = -1
        for _ in range(len(undirected)):
            nxts = [w for w in adj.get(loop[-1], ()) if w != prev]
            if not nxts:
                break
            prev, cur = loop[-1], nxts[0]
            if cur == start:
                break
            loop.append(cur)
        # 截面外法向必须指向被裁掉的一侧（+normal）；用 Newell 公式
        # 求环法向并据此定向。零面积截面（退化）直接丢弃。
        nx = ny = nz = Fraction(0)
        for i in range(len(loop)):
            p = new_verts[loop[i]]
            q = new_verts[loop[(i + 1) % len(loop)]]
            nx += (p[1] - q[1]) * (p[2] + q[2])
            ny += (p[2] - q[2]) * (p[0] + q[0])
            nz += (p[0] - q[0]) * (p[1] + q[1])
        dot = nx * normal[0] + ny * normal[1] + nz * normal[2]
        if dot < 0:
            loop.reverse()
        if dot != 0 and len(loop) >= 3:
            new_faces.append(loop)
    return new_verts, new_faces


def _polytope_volume(poly: _Polytope) -> Fraction:
    """精确体积：对一致朝外的封闭面做扇形三角化，取行列式和的绝对值。"""

    verts, faces = poly
    if not verts:
        return Fraction(0)
    total = Fraction(0)
    for face in faces:
        v0 = verts[face[0]]
        for i in range(1, len(face) - 1):
            total += _dot3(v0, _cross3(verts[face[i]], verts[face[i + 1]]))
    return abs(total) / 6


def tet_tet_intersection_volume(vertices_a, tet_a, vertices_b, tet_b) -> Fraction:
    """两个四面体交叠体积的精确值（有理数）。

    用 tet_b 的四个面半空间依次裁剪 tet_a。
    """

    poly = _tet_polytope(vertices_a, tet_a)
    pts_b = [vertices_b[i] for i in tet_b]
    for local in FACE_LOCAL:
        if not poly[0]:
            return Fraction(0)
        a, b, c = (pts_b[local[0]], pts_b[local[1]], pts_b[local[2]])
        normal = _cross3(_sub3(b, a), _sub3(c, a))
        opposite = next(k for k in range(4) if k not in local)
        if _dot3(normal, _sub3(pts_b[opposite], a)) > 0:
            normal = (-normal[0], -normal[1], -normal[2])
        poly = _clip_by_halfspace(poly, normal, _dot3(normal, a))
    return _polytope_volume(poly)


def _bbox(vertices, tet):
    pts = [vertices[i] for i in tet]
    return tuple(
        (min(p[k] for p in pts), max(p[k] for p in pts)) for k in range(3)
    )


def _bbox_overlap(box_a, box_b) -> bool:
    return all(box_a[k][0] <= box_b[k][1] and box_b[k][0] <= box_a[k][1] for k in range(3))


def _check_mesh(mesh, name: str) -> None:
    if not isinstance(mesh, TetMesh):
        raise MeshInputError(
            name + " 必须是 orient_audit.TetMesh（用 TetMesh.build 构造）", name
        )
    zero = [
        ti
        for ti, tet in enumerate(mesh.tets)
        if signed_volume6(*(mesh.vertices[i] for i in tet)) == 0
    ]
    if zero:
        raise MeshInputError(
            name + " 含零体积退化单元 " + repr(zero) + "，无法定义单元平均量", name
        )


@dataclass(frozen=True)
class RemapReport:
    """守恒重映射结果（全部为精确 Fraction）。

    属性:
        weights: 稀疏传递权重 (目标下标, 源下标, 交叠体积)，仅正体积项，
            按 (目标, 源) 升序。
        target_volumes: 每个目标单元的体积。
        covered_volumes: 每个目标单元被源网格覆盖的体积。
        uncovered_volumes: 每个目标单元未被覆盖的体积（缺口，不归一）。
        target_masses: 每个目标单元接收到的质量 sum(weights * rho_source)。
        target_averages: 每个目标单元覆盖部分的平均密度
            （target_masses / covered_volumes；完全未覆盖时为 None）。
        source_mass: 源网格总质量。
        covered_source_mass: 源网格被目标网格覆盖部分的质量。
        target_mass: 全部目标单元质量之和。
        mass_error: target_mass - covered_source_mass（守恒时精确为 0）。
    """

    weights: tuple[tuple[int, int, Fraction], ...]
    target_volumes: tuple[Fraction, ...]
    covered_volumes: tuple[Fraction, ...]
    uncovered_volumes: tuple[Fraction, ...]
    target_masses: tuple[Fraction, ...]
    target_averages: tuple
    source_mass: Fraction
    covered_source_mass: Fraction
    target_mass: Fraction
    mass_error: Fraction


def remap(source: TetMesh, source_density, target: TetMesh) -> RemapReport:
    """把源网格上的单元常量密度守恒地重映射到目标网格。

    参数:
        source: 源 TetMesh（只读）。
        source_density: 每个源单元的常量密度，有理数序列
            （int / Fraction / 有理字符串，拒绝 float 与 bool）。
        target: 目标 TetMesh（只读）。

    目标未被源完全覆盖时，缺口计入 uncovered_volumes，
    不做任何归一化；守恒口径是“目标总质量 == 源被覆盖部分质量”。
    """

    _check_mesh(source, "source")
    _check_mesh(target, "target")
    if isinstance(source_density, str) or not isinstance(source_density, (list, tuple)):
        raise MeshInputError("source_density 必须是 list/tuple 序列", "source_density")
    if len(source_density) != len(source.tets):
        raise MeshInputError(
            "source_density 长度 " + str(len(source_density))
            + " 与源单元数 " + str(len(source.tets)) + " 不一致",
            "source_density",
        )
    rho = tuple(
        parse_coord(v, "source_density[" + str(i) + "]")
        for i, v in enumerate(source_density)
    )

    src_vol = [
        abs(signed_volume6(*(source.vertices[i] for i in tet))) / 6
        for tet in source.tets
    ]
    tgt_vol = [
        abs(signed_volume6(*(target.vertices[i] for i in tet))) / 6
        for tet in target.tets
    ]
    src_box = [_bbox(source.vertices, tet) for tet in source.tets]
    tgt_box = [_bbox(target.vertices, tet) for tet in target.tets]

    weights: list[tuple[int, int, Fraction]] = []
    covered = [Fraction(0) for _ in target.tets]
    src_covered = [Fraction(0) for _ in source.tets]
    for ti, ttet in enumerate(target.tets):
        for si, stet in enumerate(source.tets):
            if not _bbox_overlap(tgt_box[ti], src_box[si]):
                continue
            vol = tet_tet_intersection_volume(
                target.vertices, ttet, source.vertices, stet
            )
            if vol > 0:
                weights.append((ti, si, vol))
                covered[ti] += vol
                src_covered[si] += vol

    masses = [Fraction(0) for _ in target.tets]
    for ti, si, vol in weights:
        masses[ti] += vol * rho[si]
    averages = tuple(
        (masses[i] / covered[i] if covered[i] > 0 else None)
        for i in range(len(target.tets))
    )
    source_mass = sum(v * r for v, r in zip(src_vol, rho))
    covered_source_mass = sum(c * r for c, r in zip(src_covered, rho))
    target_mass = sum(masses)
    return RemapReport(
        weights=tuple(weights),
        target_volumes=tuple(tgt_vol),
        covered_volumes=tuple(covered),
        uncovered_volumes=tuple(t - c for t, c in zip(tgt_vol, covered)),
        target_masses=tuple(masses),
        target_averages=averages,
        source_mass=source_mass,
        covered_source_mass=covered_source_mass,
        target_mass=target_mass,
        mass_error=target_mass - covered_source_mass,
    )
