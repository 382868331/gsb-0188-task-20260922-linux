"""守恒重映射：两张四面体网格之间的精确交叠权重。

输入两张覆盖（部分）同一区域的四面体网格与源单元常量密度，用
半空间交（两四面体共 8 个半空间）精确求每对源/目标四面体的
交叠体积，得到稀疏传递权重、每个目标单元的质量与平均密度、
未覆盖体积（缺口保留，不归一）以及总量守恒误差。

全部几何运算使用 :class:`fractions.Fraction`，无容差、无浮点判定。
输入网格对象只读；仅支持（凸的）四面体单元，不实现一般非凸单元。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key
from itertools import combinations
from typing import Any

from .errors import MeshInputError
from .exact import parse_coord
from .mesh import TetMesh, Vertex

__all__ = [
    "TransferWeight",
    "RemapResult",
    "intersection_volume",
    "tet_volume",
    "remap",
]


def _sub(a: Vertex, b: Vertex) -> Vertex:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Vertex, b: Vertex) -> Vertex:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vertex, b: Vertex) -> Fraction:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _det3(a: Vertex, b: Vertex, c: Vertex) -> Fraction:
    return _dot(a, _cross(b, c))


def _as_fraction_tet(tet) -> tuple:
    """把 4 个顶点的坐标统一转成 Fraction（接受 int/Fraction 输入）。"""

    return tuple(tuple(Fraction(c) for c in p) for p in tet)


def tet_volume(verts) -> Fraction:
    """四面体的精确（非负）体积，与顶点定向无关。"""

    verts = _as_fraction_tet(verts)
    v6 = _det3(
        _sub(verts[1], verts[0]),
        _sub(verts[2], verts[0]),
        _sub(verts[3], verts[0]),
    )
    return abs(v6) / 6


# 半空间表示：n·x <= c，n 为外法向。
_FACE_LOCAL = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (2, 1, 0))


def _tet_halfspaces(tet):
    """四面体的 4 个外法向半空间（与单元顶点定向无关）。"""

    halfspaces = []
    for k, local in enumerate(_FACE_LOCAL):
        p0, p1, p2 = tet[local[0]], tet[local[1]], tet[local[2]]
        normal = _cross(_sub(p1, p0), _sub(p2, p0))
        # 依据对顶点定向，使 normal 指向单元外侧。
        if _dot(normal, _sub(tet[k], p0)) > 0:
            normal = (-normal[0], -normal[1], -normal[2])
        halfspaces.append((normal, _dot(normal, p0)))
    return tuple(halfspaces)


def _solve3(planes):
    """精确解三张平面交点（Cramer 法则）；共面/平行时返回 None。"""

    (n0, c0), (n1, c1), (n2, c2) = planes
    d = _det3(n0, n1, n2)
    if d == 0:
        return None
    x = _det3((c0, n0[1], n0[2]), (c1, n1[1], n1[2]), (c2, n2[1], n2[2])) / d
    y = _det3((n0[0], c0, n0[2]), (n1[0], c1, n1[2]), (n2[0], c2, n2[2])) / d
    z = _det3((n0[0], n0[1], c0), (n1[0], n1[1], c1), (n2[0], n2[1], c2)) / d
    return (x, y, z)


def _plane_key(normal, c):
    """平面的规范键：法向整体缩放与符号无关（用于去重相同平面）。"""

    for comp in normal:
        if comp != 0:
            # 只归一化缩放、保留朝向：同平面异侧的半空间不是重复。
            inv = 1 / abs(comp)
            return (
                normal[0] * inv,
                normal[1] * inv,
                normal[2] * inv,
                c * inv,
            )
    raise AssertionError("零法向不构成平面")


def _unique_planes(halfspaces):
    """按规范键去重完全相同的半空间（同平面且同侧才视为重复）。"""

    seen = {}
    for normal, c in halfspaces:
        seen.setdefault(_plane_key(normal, c), (normal, c))
    return tuple(seen.values())


def _inside_all(p, halfspaces) -> bool:
    return all(_dot(n, p) <= c for n, c in halfspaces)


def _sort_on_plane(points, normal):
    """把共面点按绕外法向的正向（右手）排序，全部比较精确。"""

    if len(points) < 3:
        return points
    count = len(points)
    centroid = (
        sum(p[0] for p in points) / count,
        sum(p[1] for p in points) / count,
        sum(p[2] for p in points) / count,
    )
    ref = _sub(points[0], centroid)

    def half_of(p) -> int:
        v = _sub(p, centroid)
        cross_s = _dot(_cross(ref, v), normal)
        if cross_s > 0 or (cross_s == 0 and _dot(ref, v) > 0):
            return 0
        return 1

    def cmp(a, b) -> int:
        ha, hb = half_of(a), half_of(b)
        if ha != hb:
            return ha - hb
        cross_s = _dot(_cross(_sub(a, centroid), _sub(b, centroid)), normal)
        if cross_s > 0:
            return -1
        if cross_s < 0:
            return 1
        return 0

    return sorted(points, key=cmp_to_key(cmp))


def intersection_volume(tet_a, tet_b) -> Fraction:
    """两个（凸）四面体交叠区域的精确体积。

    方法：取两四面体共 8 个外法向半空间，枚举全部平面三元组求交点，
    保留落在所有半空间内的点作为凸交体顶点；再对每张平面上的顶点
    按外法向排序、扇形三角化，用散度定理公式
    ``V = (1/6) * Σ det(p, q, r)`` 精确求和。退化（零体积）交集自然
    得到 0，不需要任何容差。
    """

    tet_a = _as_fraction_tet(tet_a)
    tet_b = _as_fraction_tet(tet_b)
    halfspaces = _unique_planes(_tet_halfspaces(tet_a) + _tet_halfspaces(tet_b))
    vertices = set()
    for planes in combinations(halfspaces, 3):
        p = _solve3(planes)
        if p is not None and _inside_all(p, halfspaces):
            vertices.add(p)
    if len(vertices) < 4:
        return Fraction(0)

    volume6 = Fraction(0)
    for normal, c in halfspaces:
        on_plane = [p for p in vertices if _dot(normal, p) == c]
        if len(on_plane) < 3:
            continue
        ordered = _sort_on_plane(on_plane, normal)
        p0 = ordered[0]
        for k in range(1, len(ordered) - 1):
            volume6 += _det3(p0, ordered[k], ordered[k + 1])
    # 外法向约定下体积非负；退化时各项抵消为 0。
    return max(volume6, Fraction(0)) / 6


@dataclass(frozen=True)
class TransferWeight:
    """一条稀疏传递权重：源单元对目标单元的交叠体积（只存正体积项）。"""

    source_index: int
    target_index: int
    volume: Fraction


@dataclass(frozen=True)
class RemapResult:
    """守恒重映射结果（全部为精确 :class:`Fraction`）。

    属性:
        weights: 稀疏传递权重（仅交叠体积 > 0 的源/目标对）。
        target_mass: 每个目标单元接收的总质量 ``Σ w·ρ``。
        target_average_density: 目标单元平均密度，按**整个**目标单元
            体积平均（未覆盖部分按 0 贡献，不归一到满覆盖）。
        target_covered_volume: 每个目标单元被源网格覆盖的体积。
        target_uncovered_volume: 每个目标单元未被覆盖的缺口体积。
        total_target_mass: 全部目标质量之和。
        covered_source_mass: 源网格被目标覆盖部分的源质量
            （按源侧独立汇总 ``Σ ρ_i · (i 被覆盖的体积)``）。
        mass_error: ``total_target_mass - covered_source_mass``，
            精确算术中恒为 0，作为守恒性核验输出。
    """

    weights: tuple
    target_mass: tuple
    target_average_density: tuple
    target_covered_volume: tuple
    target_uncovered_volume: tuple
    total_target_mass: Fraction
    covered_source_mass: Fraction
    mass_error: Fraction


def _parse_density(raw: Any, n_source: int) -> tuple:
    if isinstance(raw, str) or not isinstance(raw, (list, tuple)):
        raise MeshInputError("source_density 必须是 list/tuple 序列", "source_density")
    if len(raw) != n_source:
        raise MeshInputError(
            f"source_density 长度 {len(raw)} 与源单元数 {n_source} 不一致",
            "source_density",
        )
    return tuple(
        parse_coord(v, f"source_density[{k}]") for k, v in enumerate(raw)
    )


def _mesh_tets(mesh: TetMesh, role: str) -> tuple:
    tets = []
    for ti, tet in enumerate(mesh.tets):
        verts = tuple(mesh.vertices[k] for k in tet)
        if tet_volume(verts) == 0:
            raise MeshInputError(
                f"{role}网格含零体积（共面退化）单元，无法参与重映射",
                f"{role}.tets[{ti}]",
            )
        tets.append(verts)
    return tuple(tets)


def _bbox(tet):
    lo = tuple(min(p[d] for p in tet) for d in range(3))
    hi = tuple(max(p[d] for p in tet) for d in range(3))
    return lo, hi


def _bbox_overlap(a, b) -> bool:
    return all(a[0][d] <= b[1][d] and b[0][d] <= a[1][d] for d in range(3))


def remap(source: TetMesh, target: TetMesh, source_density: Any) -> RemapResult:
    """把源网格上的单元常量密度守恒地重映射到目标网格。

    参数:
        source: 源 :class:`TetMesh`（只读，不被修改）。
        target: 目标 :class:`TetMesh`（只读，不被修改）。
        source_density: 每个源单元的常量密度，有理数序列
            （``int`` / ``Fraction`` / 有理字符串；拒绝 ``float``）。

    返回:
        :class:`RemapResult`。权重来自真实的四面体交叠体积，
        目标未覆盖部分保留为缺口（``target_uncovered_volume``），
        平均密度按整个目标单元体积计算，不做满覆盖归一。

    异常:
        :class:`MeshInputError`: 密度长度不符、密度非有理数、
        或任一网格含零体积单元（错误带定位）。
    """

    if not isinstance(source, TetMesh) or not isinstance(target, TetMesh):
        raise MeshInputError(
            "source/target 必须是 orient_audit.TetMesh 实例", "source/target"
        )
    density = _parse_density(source_density, len(source.tets))
    src_tets = _mesh_tets(source, "source")
    tgt_tets = _mesh_tets(target, "target")
    src_bbox = [_bbox(t) for t in src_tets]
    tgt_bbox = [_bbox(t) for t in tgt_tets]
    tgt_volume = [tet_volume(t) for t in tgt_tets]

    weights = []
    covered_by_target = [Fraction(0)] * len(tgt_tets)
    covered_by_source = [Fraction(0)] * len(src_tets)
    for ti, tgt in enumerate(tgt_tets):
        for si, src in enumerate(src_tets):
            if not _bbox_overlap(tgt_bbox[ti], src_bbox[si]):
                continue
            overlap = intersection_volume(tgt, src)
            if overlap > 0:
                weights.append(TransferWeight(si, ti, overlap))
                covered_by_target[ti] += overlap
                covered_by_source[si] += overlap

    target_mass = tuple(
        sum(
            (w.volume * density[w.source_index] for w in weights if w.target_index == ti),
            Fraction(0),
        )
        for ti in range(len(tgt_tets))
    )
    target_average = tuple(
        target_mass[ti] / tgt_volume[ti] for ti in range(len(tgt_tets))
    )
    uncovered = tuple(
        tgt_volume[ti] - covered_by_target[ti] for ti in range(len(tgt_tets))
    )
    total_target_mass = sum(target_mass, Fraction(0))
    covered_source_mass = sum(
        (density[si] * covered_by_source[si] for si in range(len(src_tets))),
        Fraction(0),
    )
    return RemapResult(
        weights=tuple(weights),
        target_mass=target_mass,
        target_average_density=target_average,
        target_covered_volume=tuple(covered_by_target),
        target_uncovered_volume=uncovered,
        total_target_mass=total_target_mass,
        covered_source_mass=covered_source_mass,
        mass_error=total_target_mass - covered_source_mass,
    )
