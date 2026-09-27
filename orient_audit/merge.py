"""两份四面体网格的显式点合并。

核心语义（见 README 与 :func:`merge_meshes`）：

- **只合并调用方明确列出的重合点对**；除此之外，哪怕两点坐标完全相同
  也保持为两个独立顶点，绝不按距离做任何自动吸附；
- 指定点对的坐标必须**精确相等**（有理比较，无容差），否则视为调用方
  错误并抛出 :class:`~orient_audit.errors.MergePointError`；
- 合并后重新做一遍完整审计：重复单元抛
  :class:`~orient_audit.errors.DuplicateTetrahedraError`、零体积抛
  :class:`~orient_audit.errors.ZeroVolumeError`，非流形 / 定向冲突 /
  连通分量等几何状态由 :class:`~orient_audit.audit.AuditReport` 报告；
- 任何一步失败都不返回半成品，两份输入网格保持原样（不可变对象 +
  新对象全部在成功后才交付）。
"""

from __future__ import annotations

from dataclasses import dataclass

from .audit import AuditReport, audit
from .errors import (
    DuplicateTetrahedraError,
    MergePointError,
    MeshInputError,
    ZeroVolumeError,
)
from .exact import MAX_TETS, signed_volume6
from .mesh import TetMesh


@dataclass(frozen=True)
class MergeResult:
    """成功合并的交付物（不可变）。

    属性:
        mesh: 合并后的新 :class:`TetMesh`（与两份输入均不同的新对象）。
        vertex_map_a: ``vertex_map_a[i]`` 为网格 A 旧顶点 ``i`` 的新编号。
        vertex_map_b: ``vertex_map_b[j]`` 为网格 B 旧顶点 ``j`` 的新编号；
                      被合并的点与对应 A 点编号相同。
        merged_pairs: 实际合并的点对 ``(a_index, b_index)``，按输入顺序。
        n_vertices_a: 网格 A 的顶点数（编号分界，便于反查）。
        report: 对合并后网格重新执行的完整 :class:`AuditReport`。
    """

    mesh: TetMesh
    vertex_map_a: tuple[int, ...]
    vertex_map_b: tuple[int, ...]
    merged_pairs: tuple[tuple[int, int], ...]
    n_vertices_a: int
    report: AuditReport

    def origin_of(self, new_index: int) -> tuple[str, int]:
        """反查新编号来源：``("a", 旧编号)`` 或 ``("b", 旧编号)``。"""

        if new_index < self.n_vertices_a:
            return ("a", new_index)
        b_old = self.vertex_map_b.index(new_index)
        return ("b", b_old)


def _check_pair_index(value, location: str, n_verts: int) -> int:
    """校验点对中的单个下标：必须是真 int（拒绝 bool）且不越界。"""

    if isinstance(value, bool) or not isinstance(value, int):
        raise MergePointError(
            f"可合并点下标必须是 int，得到 {type(value).__name__}", location
        )
    if not 0 <= value < n_verts:
        raise MergePointError(
            f"可合并点下标 {value} 越界（合法区间 [0, {n_verts - 1}]）",
            location,
        )
    return value


def merge_meshes(
    mesh_a: TetMesh,
    mesh_b: TetMesh,
    merge_pairs=(),
) -> MergeResult:
    """按调用方显式指定的点对，把两份四面体网格拼成一份。

    参数:
        mesh_a / mesh_b: 两份已由 ``TetMesh.build`` 构造的网格。
        merge_pairs: ``(a_index, b_index)`` 点对序列。**只有列在这里的
            点会被视为同一个点**；未列出的点即使坐标完全重合也各自保留，
            本函数不会因距离近而自动吸附。每对两点的坐标必须精确相等。

    返回: :class:`MergeResult`，含新网格、新旧编号关系和重新审计报告。

    异常（失败时两份输入网格均不被修改，也不留下半成品）：

        TypeError: 参数不是 :class:`TetMesh` / 点对序列结构非法时抛
            :class:`MergePointError`（后者带定位）；
        MergePointError: 点对结构非法、同一下标重复指定、或两点坐标不相等；
        DuplicateTetrahedraError: 合并后出现四点集合相同的重复单元；
        ZeroVolumeError: 合并后存在精确体积为零的单元；
        MeshInputError: 合并后单元数超过 :data:`~orient_audit.exact.MAX_TETS`。

    非流形面、定向冲突内面、连通分量等**不**作为异常，由返回的
    ``MergeResult.report`` 如实给出——它们是合并结果的拓扑状态，
    与单网格 :func:`~orient_audit.audit.audit` 的宽容语义一致。
    """

    if not isinstance(mesh_a, TetMesh):
        raise TypeError("mesh_a 必须是 TetMesh 实例（请先用 TetMesh.build 构造）")
    if not isinstance(mesh_b, TetMesh):
        raise TypeError("mesh_b 必须是 TetMesh 实例（请先用 TetMesh.build 构造）")
    if isinstance(merge_pairs, str) or not isinstance(merge_pairs, (list, tuple)):
        raise MergePointError(
            "merge_pairs 必须是 (a_index, b_index) 点对的 list/tuple",
            "merge_pairs",
        )

    na = len(mesh_a.vertices)
    nb = len(mesh_b.vertices)
    used_a: dict[int, int] = {}
    used_b: dict[int, int] = {}
    pairs: list[tuple[int, int]] = []
    for k, raw_pair in enumerate(merge_pairs):
        loc = f"merge_pairs[{k}]"
        if isinstance(raw_pair, str) or not isinstance(raw_pair, (list, tuple)):
            raise MergePointError("点对必须是两个下标的序列 (a, b)", loc)
        if len(raw_pair) != 2:
            raise MergePointError(
                f"点对必须恰好含 2 个下标，得到 {len(raw_pair)} 个", loc
            )
        ia = _check_pair_index(raw_pair[0], f"{loc}[0]", na)
        ib = _check_pair_index(raw_pair[1], f"{loc}[1]", nb)
        if ia in used_a:
            raise MergePointError(
                f"网格 A 的顶点 {ia} 已在 merge_pairs[{used_a[ia]}] 中指定，"
                "不能再与另一个 B 顶点合并（每个点至多出现一次）",
                f"{loc}[0]",
            )
        if ib in used_b:
            raise MergePointError(
                f"网格 B 的顶点 {ib} 已在 merge_pairs[{used_b[ib]}] 中指定，"
                "不能再与另一个 A 顶点合并（每个点至多出现一次）",
                f"{loc}[1]",
            )
        # 精确有理相等比较；绝无容差，更不会把“离得近”的点吸在一起。
        if mesh_a.vertices[ia] != mesh_b.vertices[ib]:
            pa = "(" + ", ".join(map(str, mesh_a.vertices[ia])) + ")"
            pb = "(" + ", ".join(map(str, mesh_b.vertices[ib])) + ")"
            raise MergePointError(
                f"指定合并的两点坐标不相等：A[{ia}]={pa} B[{ib}]={pb}；"
                "本库不按距离自动吸附，请只对位置完全相同的点建立点对",
                loc,
            )
        used_a[ia] = k
        used_b[ib] = k
        pairs.append((ia, ib))

    # ---- 编号映射：A 点保持原编号，追加 B 中未被合并的点 ----
    map_a = tuple(range(na))
    b_to_new: list[int | None] = [None] * nb
    new_vertices = list(mesh_a.vertices)  # A 顶点坐标对象原样保留
    for ia, ib in pairs:
        b_to_new[ib] = ia
    for j in range(nb):
        if b_to_new[j] is None:
            b_to_new[j] = len(new_vertices)
            new_vertices.append(mesh_b.vertices[j])
    map_b = tuple(b_to_new)  # type: ignore[arg-type]

    merged_tets = list(mesh_a.tets) + [
        tuple(map_b[k] for k in tet) for tet in mesh_b.tets
    ]

    # ---- 重建即复检：单元数上限、重复单元（含跨网格重复）。
    # 两份输入已是通过 TetMesh.build 校验的不可变对象，映射又保持“四点
    # 互异”，因此这里直接组装数据类（顶点坐标对象原样保留，不重新解析），
    # 仅把契约要求的两项整体复检显式做一遍。
    if len(merged_tets) > MAX_TETS:
        raise MeshInputError(
            f"合并后四面体数量 {len(merged_tets)} 超过上限 {MAX_TETS}", "tets"
        )
    seen: dict[tuple[int, ...], list[int]] = {}
    for ti, tet in enumerate(merged_tets):
        seen.setdefault(tuple(sorted(tet)), []).append(ti)
    dups = tuple((key, tuple(pos)) for key, pos in seen.items() if len(pos) > 1)
    if dups:
        raise DuplicateTetrahedraError(dups)
    merged = TetMesh(tuple(new_vertices), tuple(merged_tets))

    # ---- 零体积复检（合并本身不引入退化，但拼接结果必须重新证明）----
    zero: list[int] = []
    for ti, tet in enumerate(merged.tets):
        p0, p1, p2, p3 = (merged.vertices[k] for k in tet)
        if signed_volume6(p0, p1, p2, p3) == 0:
            zero.append(ti)
    if zero:
        raise ZeroVolumeError(tuple(zero))

    # ---- 合并后重新审计：非流形、内面朝向、连通分量在此重新计算 ----
    report = audit(merged)
    assert not report.has_zero_volume  # 上面已拦截

    return MergeResult(
        mesh=merged,
        vertex_map_a=map_a,
        vertex_map_b=map_b,
        merged_pairs=tuple(pairs),
        n_vertices_a=na,
        report=report,
    )
