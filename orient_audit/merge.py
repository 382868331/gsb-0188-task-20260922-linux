"""两份四面体网格的显式重合点合并。

合并是**显式 opt-in** 的：调用方必须通过 merge_pairs 逐一给出
（第一份网格顶点下标, 第二份网格顶点下标）。只有被明确配对、且坐标
**精确相等**的点才会合成一个点；其他顶点即使坐标完全相同也保持各自独立，
绝不因距离近而自动吸附。

合并会用新编号重新构造 TetMesh，因此重复单元、零体积单元与非流形面
都会在合并后的编号体系下重新检查。任一检查失败都抛异常，不返回半份
结果；两个输入网格本身不可变，失败后仍保持原状。
"""

from __future__ import annotations

from dataclasses import dataclass

from .audit import AuditReport, audit
from .errors import MeshInputError, NonManifoldError, ZeroVolumeError
from .exact import MAX_TETS
from .mesh import TetMesh

#: 来源标记：(0, 旧下标) 表示第一份网格，(1, 旧下标) 表示第二份网格。
VertexOrigin = tuple[int, int]


@dataclass(frozen=True)
class MergeMapping:
    """合并后编号与两份旧网格编号之间的双向关系。

    属性:
        first_vertex_map: 第一份网格旧顶点下标 -> 新顶点下标。
        second_vertex_map: 第二份网格旧顶点下标 -> 新顶点下标。
        vertex_origins: 新顶点下标 -> 一个或两个 (来源, 旧下标)；
                        来源 0/1 分别表示第一/第二份网格；被配对的点
                        恰有两个来源，未配对点只有一个。
        merged_vertices: 每项 (新下标, 第一份旧下标, 第二份旧下标)，
                         列出所有实际发生合并的顶点。
        first_tet_map: 第一份网格旧单元下标 -> 新单元下标。
        second_tet_map: 第二份网格旧单元下标 -> 新单元下标。
        tet_origins: 新单元下标 -> (来源, 旧单元下标)。
    """

    first_vertex_map: tuple[int, ...]
    second_vertex_map: tuple[int, ...]
    vertex_origins: tuple[tuple[VertexOrigin, ...], ...]
    merged_vertices: tuple[tuple[int, int, int], ...]
    first_tet_map: tuple[int, ...]
    second_tet_map: tuple[int, ...]
    tet_origins: tuple[VertexOrigin, ...]


@dataclass(frozen=True)
class MergeResult:
    """成功合并的结果。

    属性:
        mesh: 合并后的新 TetMesh（失败时不会产生）。
        mapping: 新旧顶点 / 单元编号关系。
        report: 对合并后网格执行的完整 AuditReport。
    """

    mesh: TetMesh
    mapping: MergeMapping
    report: AuditReport


def merge_meshes(first: TetMesh, second: TetMesh, merge_pairs=()) -> MergeResult:
    """按调用方明确给出的顶点配对合并两份网格。

    参数:
        first: 第一份网格（其顶点编号在合并结果中保持在前）。
        second: 第二份网格。
        merge_pairs: (a_index, b_index) 序列。两点必须坐标精确相等；
            同一个旧顶点至多出现一次（一对一合并）。未列出的重合点
            不会被合并，本函数不存在距离容差或自动吸附。

    返回:
        MergeResult，含合并网格、编号映射与合并后审计报告。

    失败语义（均不产生部分结果，两个原网格保持不变）:

    - 配对非法（类型、越界、重复、坐标不精确相等）->
      MeshInputError，错误位置如 merge_pairs[2][0]；
    - 合并后单元数超过 MAX_TETS -> MeshInputError；
    - 合并后出现重复四面体 -> DuplicateTetrahedraError；
    - 合并后存在零体积四面体 -> ZeroVolumeError；
    - 合并后存在引用数 > 2 的非流形面 -> NonManifoldError。

    单元负向或两单元共享面朝向不一致不在此拒绝：它们是可由 audit
    报告、repair_orientation 修复的定向状态。
    """

    if not isinstance(first, TetMesh) or not isinstance(second, TetMesh):
        raise TypeError(
            "merge_meshes 需要两个 TetMesh 实例（请先用 TetMesh.build 构造）"
        )
    if isinstance(merge_pairs, str) or not isinstance(merge_pairs, (list, tuple)):
        raise MeshInputError(
            "merge_pairs 必须是 (第一份下标, 第二份下标) 的 list/tuple",
            "merge_pairs",
        )
    if len(first.tets) + len(second.tets) > MAX_TETS:
        raise MeshInputError(
            f"合并后四面体数量 {len(first.tets) + len(second.tets)} "
            f"超过上限 {MAX_TETS}",
            "tets",
        )

    paired_b_to_a = _validate_pairs(first, second, merge_pairs)

    # ---- 组装新顶点：第一份原样在前；第二份配对点复用，其余追加 ----
    new_vertices = list(first.vertices)
    origins: list[list[VertexOrigin]] = [[(0, i)] for i in range(len(first.vertices))]
    second_to_new = [-1] * len(second.vertices)
    merged_vertices: list[tuple[int, int, int]] = []
    for b_index, coord in enumerate(second.vertices):
        if b_index in paired_b_to_a:
            a_index = paired_b_to_a[b_index]
            second_to_new[b_index] = a_index
            origins[a_index].append((1, b_index))
            merged_vertices.append((a_index, a_index, b_index))
        else:
            new_index = len(new_vertices)
            second_to_new[b_index] = new_index
            new_vertices.append(coord)
            origins.append([(1, b_index)])

    new_tets = [tuple(tet) for tet in first.tets]
    for tet in second.tets:
        new_tets.append(tuple(second_to_new[i] for i in tet))

    # 在新编号体系下重新构造：结构校验与重复单元复查都在这里发生。
    merged = TetMesh.build(new_vertices, new_tets)
    report = audit(merged)

    if report.zero_tets:
        raise ZeroVolumeError(report.zero_tets)
    if report.non_manifold_faces:
        raise NonManifoldError(
            tuple((face.key, face.tet_indices) for face in report.non_manifold_faces)
        )

    n_first_tets = len(first.tets)
    mapping = MergeMapping(
        first_vertex_map=tuple(range(len(first.vertices))),
        second_vertex_map=tuple(second_to_new),
        vertex_origins=tuple(tuple(items) for items in origins),
        merged_vertices=tuple(merged_vertices),
        first_tet_map=tuple(range(n_first_tets)),
        second_tet_map=tuple(range(n_first_tets, n_first_tets + len(second.tets))),
        tet_origins=tuple(
            [(0, i) for i in range(n_first_tets)]
            + [(1, i) for i in range(len(second.tets))]
        ),
    )
    return MergeResult(merged, mapping, report)


def _validate_pairs(
    first: TetMesh, second: TetMesh, merge_pairs
) -> dict[int, int]:
    """校验显式配对，返回 第二份顶点下标 -> 第一份顶点下标 表。"""

    b_to_a: dict[int, int] = {}
    used_a: set[int] = set()

    for k, pair in enumerate(merge_pairs):
        loc = f"merge_pairs[{k}]"
        if isinstance(pair, str) or not isinstance(pair, (list, tuple)):
            raise MeshInputError("每个配对必须是 (a_index, b_index) 序列", loc)
        if len(pair) != 2:
            raise MeshInputError(
                f"每个配对必须恰好含 2 个顶点下标，得到 {len(pair)} 个", loc
            )
        a_index = _check_pair_index(pair[0], f"{loc}[0]", first, "第一份网格")
        b_index = _check_pair_index(pair[1], f"{loc}[1]", second, "第二份网格")
        if a_index in used_a:
            raise MeshInputError(
                f"第一份网格顶点 {a_index} 已参与其他配对；合并必须一对一", loc
            )
        if b_index in b_to_a:
            raise MeshInputError(
                f"第二份网格顶点 {b_index} 已参与其他配对；合并必须一对一", loc
            )
        if first.vertices[a_index] != second.vertices[b_index]:
            raise MeshInputError(
                f"配对顶点坐标不精确相等：第一份 verts[{a_index}]="
                f"{first.vertices[a_index]}，第二份 verts[{b_index}]="
                f"{second.vertices[b_index]}；不会按距离自动吸附",
                loc,
            )
        used_a.add(a_index)
        b_to_a[b_index] = a_index

    return b_to_a


def _check_pair_index(value, location: str, mesh: TetMesh, mesh_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise MeshInputError(
            f"顶点下标必须是 int，得到 {type(value).__name__}", location
        )
    count = len(mesh.vertices)
    if not 0 <= value < count:
        raise MeshInputError(
            f"{mesh_name}顶点下标 {value} 越界（顶点数 {count}，"
            f"合法区间 [0, {count - 1}]）",
            location,
        )
    return value

