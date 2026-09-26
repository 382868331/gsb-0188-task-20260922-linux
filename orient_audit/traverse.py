"""线段穿越采样与积分：:func:`traverse_segment`。

调用方为每个顶点给一个标量（有理数），本模块计算线段 ``P(t) = p0 + t*(p1-p0)``
（``t ∈ [0, 1]``）穿过四面体网格的完整信息：

- 线段依次落入哪些四面体、每段的参数区间 ``[t0, t1]``；
- 每段端点处由所在单元重心插值得到的场值（精确 :class:`Fraction`）；
- 每段积分与总积分（场在单元内线性，故每段是精确梯形）；
- 网格外的间隙单独列为 ``gaps``，**不**自动补 0 混进结果。

实现要点（契约见 README）：

- **邻接行走**：先用面关联表建立 ``(单元, 面) -> 相邻单元`` 映射，
  正常跨面时直接用邻接表走到下一个单元，不重新扫描全网格；
  只有"起点定位"、"离开网格后的再进入"、"恰好穿过共享棱/顶点的退化事件"
  才退化为一次全网格定位扫描；
- **确定去重**：线恰过共享棱/顶点时，候选下一单元取**最小单元下标**；
  仅与线相切（区间长度为零）的单元直接丢弃，同一长度不会计两次；
- **沿边界走的归属**：线段在某段区间内同时位于多个单元的公共面/棱上时，
  该区间固定归属下标最小的候选单元，每个参数值只属于一个分段；
- **无效网格拒绝**：负体积、零体积、非流形面、定向冲突内面任一存在时
  抛 :class:`InvalidTraversalMeshError`，不输出看似完整的穿越序列；
  负定向单元请先按原审计规则用 ``repair_orientation`` 修复后再调用。

所有几何判定与插值、积分都在 :class:`fractions.Fraction` 上精确完成。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .audit import audit
from .errors import InvalidTraversalMeshError, MeshInputError
from .exact import parse_coord, signed_volume6
from .mesh import FACE_LOCAL, TetMesh, build_face_incidence

_ZERO = Fraction(0)
_ONE = Fraction(1)


@dataclass(frozen=True)
class SegmentPiece:
    """线段落在单个四面体内的一段。

    属性:
        tet_index: 所属单元下标。
        t0, t1: 参数区间端点（``0 <= t0 < t1 <= 1``，Fraction）。
        value0, value1: 两端点处由该单元重心坐标插值得到的场值。
        integral: 该段积分 ``(t1 - t0) * (value0 + value1) / 2``（精确）。
    """

    tet_index: int
    t0: Fraction
    t1: Fraction
    value0: Fraction
    value1: Fraction
    integral: Fraction


@dataclass(frozen=True)
class SegmentGap:
    """线段处于网格外的一段参数区间（不计入积分，不补 0）。"""

    t0: Fraction
    t1: Fraction


@dataclass(frozen=True)
class SegmentTraversal:
    """一次线段穿越的完整结果（不可变）。

    属性:
        pieces: 有序的分段（按参数升序，互不重叠）。
        gaps: 有序的网格外间隙（与 pieces 互补，合并后覆盖 [0, 1]）。
        total_integral: 全部分段积分之和（精确 Fraction）。
    """

    pieces: tuple[SegmentPiece, ...]
    gaps: tuple[SegmentGap, ...]
    total_integral: Fraction

    def total_integral_float(self) -> float:
        """总积分的一次正确舍入 double 转换（误差 ≤ 0.5 ULP）。"""

        return float(self.total_integral)

    @property
    def tet_sequence(self) -> tuple[int, ...]:
        """依次穿过的单元下标。"""

        return tuple(p.tet_index for p in self.pieces)


@dataclass(frozen=True)
class _TetLineData:
    """单个四面体与线段的精确关系（缓存）。

    ``coeffs[i] = (a, b)`` 表示把单元第 i 个顶点换成 P(t) 后的
    6 倍有向体积 ``f_i(t) = a + b*t``；点在单元内当且仅当四个
    ``f_i(t) >= 0``（单元已验证为正体积）。``lo/hi`` 为裁剪到
    [0, 1] 后的参数区间；区间为空时整个对象为 ``None``（不缓存本对象）。
    """

    vol6: Fraction
    coeffs: tuple[tuple[Fraction, Fraction], ...]
    lo: Fraction
    hi: Fraction


def _parse_point(raw: Any, name: str) -> tuple[Fraction, Fraction, Fraction]:
    if isinstance(raw, str) or not isinstance(raw, (list, tuple)):
        raise MeshInputError("线段端点必须是长度为 3 的坐标序列", name)
    coords = list(raw)
    if len(coords) != 3:
        raise MeshInputError(
            f"线段端点必须恰好有 3 个坐标分量，得到 {len(coords)} 个", name
        )
    return tuple(  # type: ignore[return-value]
        parse_coord(c, f"{name}[{k}]") for k, c in enumerate(coords)
    )


def _parse_values(values: Any, n_verts: int) -> tuple[Fraction, ...]:
    if isinstance(values, str) or not isinstance(values, (list, tuple)):
        raise MeshInputError("values 必须是 list/tuple 形式的逐顶点标量序列", "values")
    if len(values) != n_verts:
        raise MeshInputError(
            f"values 长度 {len(values)} 与顶点数 {n_verts} 不一致", "values"
        )
    return tuple(parse_coord(v, f"values[{i}]") for i, v in enumerate(values))


def _tet_line_data(mesh: TetMesh, ti: int, p0, p1) -> _TetLineData | None:
    """计算单元与线段（裁剪到 t∈[0,1]）的精确交集区间；空则 None。"""

    tet = mesh.tets[ti]
    pts = [mesh.vertices[k] for k in tet]
    vol6 = signed_volume6(pts[0], pts[1], pts[2], pts[3])
    # 调用前已拒绝零/负体积单元，此处是不变量断言。
    assert vol6 > 0

    lo, hi = _ZERO, _ONE
    coeffs: list[tuple[Fraction, Fraction]] = []
    for i in range(4):
        rep0 = list(pts)
        rep0[i] = p0
        rep1 = list(pts)
        rep1[i] = p1
        a = signed_volume6(rep0[0], rep0[1], rep0[2], rep0[3])
        b = signed_volume6(rep1[0], rep1[1], rep1[2], rep1[3]) - a
        coeffs.append((a, b))
        # 约束 f_i(t) = a + b*t >= 0
        if b > 0:
            bound = -a / b
            if bound > lo:
                lo = bound
        elif b < 0:
            bound = -a / b
            if bound < hi:
                hi = bound
        elif a < 0:
            return None
        if lo > hi:
            return None
    return _TetLineData(vol6, tuple(coeffs), lo, hi)


def _build_adjacency(mesh: TetMesh) -> dict[tuple[int, tuple[int, int, int]], int]:
    """(单元下标, 无序三点面键) -> 共享该面的另一个单元。"""

    adj: dict[tuple[int, tuple[int, int, int]], int] = {}
    for face in build_face_incidence(mesh).values():
        if face.multiplicity == 2:
            o1, o2 = face.occurrences
            adj[(o1.tet_index, face.key)] = o2.tet_index
            adj[(o2.tet_index, face.key)] = o1.tet_index
    return adj


def traverse_segment(
    mesh: TetMesh, values: Any, start: Any, end: Any
) -> SegmentTraversal:
    """线段穿越采样与积分入口。

    参数:
        mesh: 已校验的 :class:`TetMesh`；必须正体积、流形且定向一致，
              否则抛 :class:`InvalidTraversalMeshError`（负定向请先
              ``repair_orientation``）。
        values: 逐顶点标量（int/Fraction/有理字符串），长度等于顶点数。
        start, end: 线段两端点，各为 3 个有理坐标。

    返回:
        :class:`SegmentTraversal`，含有序分段、网格外间隙与总积分。
    """

    if not isinstance(mesh, TetMesh):
        raise TypeError("traverse_segment 需要 TetMesh 实例（请先用 TetMesh.build 构造）")

    report = audit(mesh)
    if (
        report.negative_tets
        or report.zero_tets
        or report.non_manifold_faces
        or report.inconsistent_faces
    ):
        raise InvalidTraversalMeshError(
            negative=report.negative_tets,
            zero=report.zero_tets,
            non_manifold=tuple(f.key for f in report.non_manifold_faces),
            inconsistent=tuple(f.key for f in report.inconsistent_faces),
        )

    vals = _parse_values(values, len(mesh.vertices))
    p0 = _parse_point(start, "start")
    p1 = _parse_point(end, "end")

    adjacency = _build_adjacency(mesh)
    face_keys = tuple(
        tuple(
            tuple(sorted((tet[f[0]], tet[f[1]], tet[f[2]])))
            for f in FACE_LOCAL
        )
        for tet in mesh.tets
    )

    cache: dict[int, _TetLineData | None] = {}

    def data(ti: int) -> _TetLineData | None:
        if ti not in cache:
            cache[ti] = _tet_line_data(mesh, ti, p0, p1)
        return cache[ti]

    def locate(t: Fraction) -> int | None:
        """全网格定位：找 lo <= t < hi 的最小下标单元（仅定位事件调用）。"""

        for ti in range(len(mesh.tets)):
            d = data(ti)
            if d is not None and d.lo <= t < d.hi:
                return ti
        return None

    def next_entry(t: Fraction) -> tuple[Fraction, int] | None:
        """t 之后最近的网格入口（间隙结束时调用一次）。"""

        best: tuple[Fraction, int] | None = None
        for ti in range(len(mesh.tets)):
            d = data(ti)
            if d is None or d.hi <= d.lo:
                continue
            if d.lo > t and (best is None or (d.lo, ti) < best):
                best = (d.lo, ti)
        return best

    def value_at(ti: int, d: _TetLineData, t: Fraction) -> Fraction:
        tet = mesh.tets[ti]
        total = _ZERO
        for i in range(4):
            a, b = d.coeffs[i]
            total += (a + b * t) * vals[tet[i]]
        return total / d.vol6

    pieces: list[SegmentPiece] = []
    gaps: list[SegmentGap] = []
    total = _ZERO

    t = _ZERO
    cur = locate(t)
    while t < _ONE:
        if cur is None:
            nxt = next_entry(t)
            if nxt is None:
                gaps.append(SegmentGap(t, _ONE))
                break
            t_next, cur = nxt
            gaps.append(SegmentGap(t, t_next))
            t = t_next
            continue

        d = data(cur)
        assert d is not None
        t_enter = t if t > d.lo else d.lo
        t_exit = d.hi
        v0 = value_at(cur, d, t_enter)
        v1 = value_at(cur, d, t_exit)
        integ = (t_exit - t_enter) * (v0 + v1) / 2
        pieces.append(SegmentPiece(cur, t_enter, t_exit, v0, v1, integ))
        total += integ
        t = t_exit
        if t >= _ONE:
            break

        # 跨面行走：在 t 处 f_i(t)==0 的面即出口面，查邻接表得候选。
        candidates: set[int] = set()
        for i in range(4):
            a, b = d.coeffs[i]
            if a + b * t == 0:
                nb = adjacency.get((cur, face_keys[cur][i]))
                if nb is None:
                    continue
                nd = data(nb)
                if nd is not None and nd.lo <= t < nd.hi:
                    candidates.add(nb)
        if candidates:
            # 确定规则：多个候选（线恰过共享棱/顶点）取最小下标。
            cur = min(candidates)
        else:
            # 退化出口（恰过棱/顶点但邻接面不连续）或离开网格：定位一次。
            cur = locate(t)

    return SegmentTraversal(tuple(pieces), tuple(gaps), total)
