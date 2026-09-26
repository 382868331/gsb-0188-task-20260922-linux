"""线段穿越采样与积分：sample_segment。

调用方为每个顶点给一个标量（分段线性场），给定线段起点/终点，本模块
精确（Fraction，无浮点、无容差）给出：

- 线段依次落入哪些四面体（SegmentHit：单元下标、参数区间 [t0, t1]、
  线性插值端点值 value0/value1、该段积分）；
- 网格外间隙（SegmentGap），单独列出，不自动补 0 混进结果；
- 总积分（各命中段积分之和，对分段线性场是精确值）。

行走策略：从当前单元经**面邻接**跨到相邻单元，不在每一步重扫全网格。
全网格扫描只发生两次：初始定位（t=0 或首次进入）与进入间隙后的
下一次进入定位。跨面本身只查面邻接表。

退化情形的确定规则（精确有理数下无容差歧义）：

- 线恰过共享棱/顶点（同一时刻 t 有多个面同时到期）：候选下一单元取
  所有到期面的邻居中"区间覆盖 t 之后"者，按下标最小者唯一选定；
  邻居表无合法候选时（线穿顶点进入仅共点的单元），回退到一次性
  预计算的区间表取同样条件下下标最小者；
- 长度为零的接触（擦过顶点/棱）不产生命中段，不计长度；
- 沿边界/共享面行走的区间归属唯一：每个参数 t 至多属于一个命中段，
  由"下标最小合法单元"规则确定归属，同一长度不会计两次。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .audit import audit
from .errors import MeshInputError, TraversalError
from .exact import parse_coord, signed_volume6
from .mesh import TetMesh, build_face_incidence

Vec3 = tuple[Fraction, Fraction, Fraction]


@dataclass(frozen=True)
class SegmentHit:
    """线段在一个四面体内的一段。

    属性:
        tet_index: 单元下标。
        t0, t1: 参数区间（Fraction，0 <= t0 < t1 <= 1，半开归属见模块文档）。
        value0, value1: 分段线性场在 t0 / t1 处的精确插值。
    """

    tet_index: int
    t0: Fraction
    t1: Fraction
    value0: Fraction
    value1: Fraction

    @property
    def integral(self) -> Fraction:
        """该段积分（线性场梯形公式即精确值）。"""

        return (self.t1 - self.t0) * (self.value0 + self.value1) / 2


@dataclass(frozen=True)
class SegmentGap:
    """线段在网格外的一段参数区间（不补 0，不混入积分）。"""

    t0: Fraction
    t1: Fraction


@dataclass(frozen=True)
class SegmentTraversal:
    """一次线段穿越的完整结果（不可变）。

    hits 按参数 t 升序、互不重叠；gaps 与 hits 互补地覆盖 [0, 1]。
    total_integral 为各命中段积分之和（Fraction 精确值）。
    """

    hits: tuple[SegmentHit, ...]
    gaps: tuple[SegmentGap, ...]
    total_integral: Fraction

    def total_integral_float(self) -> float:
        """精确总积分的一次正确舍入 double 转换（误差 <= 0.5 ULP）。"""

        return float(self.total_integral)


class _TetLine:
    """一个四面体与直线的精确关系：lambda_i(t) = (n_i + t*m_i) / vol6。"""

    __slots__ = ("vol6", "n", "m", "lo", "hi")

    def __init__(self, vol6: Fraction, n: list[Fraction], m: list[Fraction]) -> None:
        self.vol6 = vol6
        self.n = n
        self.m = m
        lo: Fraction | None = None
        hi: Fraction | None = None
        for ni, mi in zip(n, m, strict=True):
            if mi == 0:
                if ni < 0:
                    lo, hi = Fraction(1), Fraction(0)  # 空区间
                    break
                continue
            root = -ni / mi
            if mi > 0:
                if lo is None or root > lo:
                    lo = root
            else:
                if hi is None or root < hi:
                    hi = root
        self.lo = lo if lo is not None else Fraction(-10**30)
        self.hi = hi if hi is not None else Fraction(10**30)

    @property
    def nonempty(self) -> bool:
        return self.lo < self.hi

    def value_at(self, t: Fraction, field: tuple[Fraction, ...], tet) -> Fraction:
        total = Fraction(0)
        for i in range(4):
            lam = (self.n[i] + t * self.m[i]) / self.vol6
            total += lam * field[tet[i]]
        return total


def _parse_point(raw, location: str) -> Vec3:
    if isinstance(raw, str) or not isinstance(raw, (list, tuple)):
        raise MeshInputError("线段端点必须是长度为 3 的坐标序列", location)
    if len(raw) != 3:
        raise MeshInputError(
            f"线段端点必须恰好有 3 个坐标分量，得到 {len(raw)} 个", location
        )
    return tuple(  # type: ignore[return-value]
        parse_coord(c, f"{location}[{k}]") for k, c in enumerate(raw)
    )


def _parse_values(values, n_verts: int) -> tuple[Fraction, ...]:
    if isinstance(values, str) or not isinstance(values, (list, tuple)):
        raise MeshInputError("values 必须是每顶点一个标量的序列", "values")
    if len(values) != n_verts:
        raise MeshInputError(
            f"values 长度 {len(values)} 与顶点数 {n_verts} 不一致", "values"
        )
    return tuple(parse_coord(v, f"values[{k}]") for k, v in enumerate(values))


def _check_traversable(mesh: TetMesh) -> None:
    report = audit(mesh)
    problems: list[str] = []
    if report.zero_tets:
        problems.append(
            "零体积单元 " + ", ".join(f"tets[{i}]" for i in report.zero_tets)
        )
    if report.negative_tets:
        problems.append(
            "负定向单元 "
            + ", ".join(f"tets[{i}]" for i in report.negative_tets)
            + "（请先 repair_orientation）"
        )
    if report.non_manifold_faces:
        problems.append(
            "非流形面 " + ", ".join(str(f.key) for f in report.non_manifold_faces)
        )
    if problems:
        raise TraversalError(tuple(problems))


def sample_segment(mesh: TetMesh, values, start, end) -> SegmentTraversal:
    """沿线段 start->end 采样每顶点标量场并积分（全部精确有理）。

    参数:
        mesh: TetMesh（须通过可行走性检查，否则抛 TraversalError）。
        values: 每顶点一个标量（int / Fraction / 有理字符串；拒绝 float/bool）。
        start, end: 线段端点，各为 3 个有理坐标；t=0 在 start，t=1 在 end。

    返回 SegmentTraversal；网格外区间在 gaps 中单独列出。
    """

    if not isinstance(mesh, TetMesh):
        raise TypeError("sample_segment 需要 TetMesh 实例（请先用 TetMesh.build 构造）")
    _check_traversable(mesh)

    field = _parse_values(values, len(mesh.vertices))
    a = _parse_point(start, "start")
    b = _parse_point(end, "end")
    d = tuple(b[k] - a[k] for k in range(3))
    if d == (Fraction(0), Fraction(0), Fraction(0)):
        raise MeshInputError("线段两端点重合，无法定义参数化", "end")

    # ---- 一次性预计算：每单元的直线区间与重心分子（之后跨面不再重扫） ----
    others = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (2, 1, 0))
    lines: list[_TetLine] = []
    for tet in mesh.tets:
        p = [mesh.vertices[k] for k in tet]
        vol6 = signed_volume6(p[0], p[1], p[2], p[3])
        # vol6 > 0 已由 _check_traversable 保证（无零体积、无负向）。
        n_at_a: list[Fraction] = []
        n_at_b: list[Fraction] = []
        for o in others:
            q = (p[o[0]], p[o[1]], p[o[2]])
            n_at_a.append(signed_volume6(a, q[0], q[1], q[2]))
            n_at_b.append(signed_volume6(b, q[0], q[1], q[2]))
        lines.append(
            _TetLine(
                vol6,
                n_at_a,
                [nb - na for na, nb in zip(n_at_a, n_at_b, strict=True)],
            )
        )

    active = [ti for ti, ln in enumerate(lines) if ln.nonempty]
    active_set = set(active)

    # ---- 面邻接表（仅用于跨面行走） ----
    incidence = build_face_incidence(mesh)
    neighbor: dict[tuple[int, tuple[int, int, int]], int] = {}
    for face in incidence.values():
        if face.multiplicity == 2:
            (o1, o2) = face.occurrences
            neighbor[(o1.tet_index, face.key)] = o2.tet_index
            neighbor[(o2.tet_index, face.key)] = o1.tet_index

    def face_key_excluding(tet, i: int) -> tuple[int, int, int]:
        return tuple(sorted(tet[j] for j in range(4) if j != i))  # type: ignore[return-value]

    def covers(ti: int, t: Fraction) -> bool:
        ln = lines[ti]
        return ln.lo <= t < ln.hi

    def pick(candidates, t: Fraction) -> int | None:
        valid = [ti for ti in candidates if ti in active_set and covers(ti, t)]
        return min(valid) if valid else None

    def next_entry(t: Fraction) -> int | None:
        """t 之后（含 t）最早进入的单元（用一次性预计算的区间表）。"""

        best: int | None = None
        best_lo: Fraction | None = None
        for ti in active:
            ln = lines[ti]
            enter = ln.lo if ln.lo > t else t
            if enter < ln.hi and (best_lo is None or enter < best_lo):
                best, best_lo = ti, enter
        return best

    hits: list[SegmentHit] = []
    gaps: list[SegmentGap] = []
    one = Fraction(1)
    zero = Fraction(0)

    t = zero
    current = pick(active, zero)
    if current is None:
        current = next_entry(zero)

    while t < one:
        if current is None:
            break
        ln = lines[current]
        enter = ln.lo if ln.lo > t else t
        if enter >= one:
            break
        if enter > t:
            gaps.append(SegmentGap(t, min(enter, one)))
            t = enter
        exit_t = min(ln.hi, one)
        tet = mesh.tets[current]
        hits.append(
            SegmentHit(
                tet_index=current,
                t0=t,
                t1=exit_t,
                value0=ln.value_at(t, field, tet),
                value1=ln.value_at(exit_t, field, tet),
            )
        )
        if exit_t >= one:
            t = one
            break
        # ---- 跨面行走：找出同时到期的面（线恰过棱/顶点时不止一个） ----
        due: list[int] = []
        for i in range(4):
            mi, ni = ln.m[i], ln.n[i]
            if mi < 0 and -ni / mi == exit_t:
                due.append(i)
        candidates = []
        for i in due:
            key = face_key_excluding(tet, i)
            nb = neighbor.get((current, key))
            if nb is not None:
                candidates.append(nb)
        nxt = pick(candidates, exit_t)
        if nxt is None:
            # 线恰穿仅共点/共棱的单元：回退到预计算区间表（仍非逐步重扫）。
            nxt = pick((ti for ti in active if ti != current), exit_t)
        if nxt is None:
            # 离开网格：定位下一次进入。
            follow = next_entry(exit_t)
            if follow is None or lines[follow].lo >= one:
                gaps.append(SegmentGap(exit_t, one))
                t = one
                current = None
                break
            t = exit_t
            current = follow
        else:
            t = exit_t
            current = nxt

    if t < one and (current is None or not covers(current, t)):
        # 起点在网格外且之后无进入：整段为间隙。
        if not gaps or gaps[-1].t1 < one:
            gaps.append(SegmentGap(t, one))

    total = sum((h.integral for h in hits), Fraction(0))
    return SegmentTraversal(tuple(hits), tuple(gaps), total)

