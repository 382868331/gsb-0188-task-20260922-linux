"""定向与边界审计核心：:func:`audit` 与 :func:`repair_orientation`。

契约摘要（详见 README）：

- 体积极性用精确有理行列式判定；负向单元修复时**只交换最后两个顶点**
  （下标位 2 与 3），正体积单元不动，零体积单元拒绝修复；
- 面按无序三点键关联：出现 1 次为边界面，恰 2 次且诱导方向相反为
  合法内面，恰 2 次同向为定向冲突内面，>2 次为非流形面；
- 修复不移动任何顶点（返回新网格，顶点对象恒等保留）。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .errors import ZeroVolumeError
from .exact import signed_volume6
from .mesh import (
    FaceIncidence,
    TetMesh,
    build_face_incidence,
)


@dataclass(frozen=True)
class BoundaryFace:
    """边界面（仅被一个单元引用）。

    属性:
        ordered: 该单元诱导朝向的三个顶点下标；对修复后的网格，
                 此即外法向朝向。
        tet_index: 所属单元下标。
    """

    ordered: tuple[int, int, int]
    tet_index: int


@dataclass(frozen=True)
class InteriorFace:
    """恰被两个单元引用的内面（无论方向是否一致）。"""

    key: tuple[int, int, int]
    tet_indices: tuple[int, int]
    opposite: bool


@dataclass(frozen=True)
class NonManifoldFace:
    """被超过两个单元引用的面。"""

    key: tuple[int, int, int]
    tet_indices: tuple[int, ...]


@dataclass(frozen=True)
class AuditReport:
    """一次审计的完整结果（不可变）。

    属性:
        n_tets: 单元数。
        negative_tets: 精确行列式为负的单元下标（按升序）。
        zero_tets: 精确体积为零（共面退化）的单元下标（按升序）。
        total_signed_volume6: 全部单元 6 倍有向体积的精确和（Fraction）。
        boundary_faces: 边界面，带单元诱导朝向。
        interior_faces: 所有恰被两单元引用的面；``opposite`` 标明朝向是否相反。
        non_manifold_faces: 引用数 > 2 的面。
        components: 按共面邻接（共享任一三点面即相连）得到的单元连通分量，
                    每个分量为单元下标元组，按最小下标排序。
    """

    n_tets: int
    negative_tets: tuple[int, ...]
    zero_tets: tuple[int, ...]
    total_signed_volume6: Fraction
    boundary_faces: tuple[BoundaryFace, ...]
    interior_faces: tuple[InteriorFace, ...]
    non_manifold_faces: tuple[NonManifoldFace, ...]
    components: tuple[tuple[int, ...], ...]

    # ---- 派生判定 ----
    @property
    def total_signed_volume(self) -> Fraction:
        """全部单元有向体积的精确和（``total_signed_volume6 / 6``）。"""

        return self.total_signed_volume6 / 6

    def total_signed_volume_float(self) -> float:
        """精确有向体积和的一次正确舍入 double 转换。

        误差口径见 README：误差不超过 0.5 ULP，量级正常的数相对误差
        <= 2**-53（约 1.11e-16）。该舍入误差与几何/拓扑判定无关。
        """

        return float(self.total_signed_volume)

    @property
    def inconsistent_faces(self) -> tuple[InteriorFace, ...]:
        """两单元共享但诱导方向相同的内面（定向冲突，可定位）。"""

        return tuple(f for f in self.interior_faces if not f.opposite)

    @property
    def has_zero_volume(self) -> bool:
        return bool(self.zero_tets)

    @property
    def has_non_manifold(self) -> bool:
        return bool(self.non_manifold_faces)

    @property
    def consistently_oriented(self) -> bool:
        """所有单元正体积且每个内面两单元诱导方向相反。"""

        return not self.negative_tets and not self.inconsistent_faces

    @property
    def is_manifold(self) -> bool:
        """没有引用数 >2 的面（不含定向判定）。"""

        return not self.non_manifold_faces


def audit(mesh: TetMesh) -> AuditReport:
    """对网格执行定向与边界审计，纯只读，不修改输入。"""

    if not isinstance(mesh, TetMesh):
        raise TypeError("audit 需要 TetMesh 实例（请先用 TetMesh.build 构造）")

    vol6_by_tet: list[Fraction] = []
    negative: list[int] = []
    zero: list[int] = []
    total6 = Fraction(0)
    for ti, tet in enumerate(mesh.tets):
        pts = tuple(mesh.vertices[k] for k in tet)
        v6 = signed_volume6(pts[0], pts[1], pts[2], pts[3])
        vol6_by_tet.append(v6)
        total6 += v6
        if v6 < 0:
            negative.append(ti)
        elif v6 == 0:
            zero.append(ti)

    incidence = build_face_incidence(mesh)
    boundary: list[BoundaryFace] = []
    interior: list[InteriorFace] = []
    nonmanifold: list[NonManifoldFace] = []
    for face in incidence.values():
        if face.is_boundary:
            occ = face.occurrences[0]
            boundary.append(BoundaryFace(occ.ordered, occ.tet_index))
        elif face.multiplicity == 2:
            o1, o2 = face.occurrences
            interior.append(
                InteriorFace(
                    face.key,
                    (o1.tet_index, o2.tet_index),
                    opposite=face.is_manifold_interior,
                )
            )
        else:  # > 2
            nonmanifold.append(
                NonManifoldFace(
                    face.key,
                    tuple(o.tet_index for o in face.occurrences),
                )
            )

    components = _connected_components(mesh.tets, incidence)

    return AuditReport(
        n_tets=len(mesh.tets),
        negative_tets=tuple(negative),
        zero_tets=tuple(zero),
        total_signed_volume6=total6,
        boundary_faces=tuple(boundary),
        interior_faces=tuple(interior),
        non_manifold_faces=tuple(nonmanifold),
        components=components,
    )


def _connected_components(
    tets: tuple[tuple[int, int, int, int], ...],
    incidence: dict[tuple[int, int, int], FaceIncidence],
) -> tuple[tuple[int, ...], ...]:
    """并查集：共享任一三点面的单元归为同一连通分量。"""

    parent = list(range(len(tets)))

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != x:
            parent[x], x = root, parent[x]
        return root

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for face in incidence.values():
        if face.multiplicity >= 2:
            first = face.occurrences[0].tet_index
            for occ in face.occurrences[1:]:
                union(first, occ.tet_index)

    groups: dict[int, list[int]] = {}
    for i in range(len(tets)):
        groups.setdefault(find(i), []).append(i)
    return tuple(tuple(sorted(members)) for members in sorted(groups.values()))


def repair_orientation(mesh: TetMesh) -> TetMesh:
    """修复单元定向：负体积单元只交换最后两个顶点（位 2 与位 3）。

    - 正体积单元保持不变；
    - 零体积（共面退化）单元无法通过交换顶点修复，抛出
      :class:`ZeroVolumeError`，且不产生任何部分修改；
    - 顶点坐标对象恒等保留，修复不移动顶点；
    - 非流形 / 定向冲突是拓扑状态，由 :func:`audit` 报告，不阻止本修复
      （本函数的职责仅是逐单元体积极性）。
    """

    if not isinstance(mesh, TetMesh):
        raise TypeError(
            "repair_orientation 需要 TetMesh 实例（请先用 TetMesh.build 构造）"
        )

    zero: list[int] = []
    repaired: list[tuple[int, int, int, int]] = []
    for ti, tet in enumerate(mesh.tets):
        pts = tuple(mesh.vertices[k] for k in tet)
        v6 = signed_volume6(pts[0], pts[1], pts[2], pts[3])
        if v6 == 0:
            zero.append(ti)
            repaired.append(tet)
        elif v6 < 0:
            # 奇置换：交换最后两个顶点，行列式变号。
            repaired.append((tet[0], tet[1], tet[3], tet[2]))
        else:
            repaired.append(tet)

    if zero:
        # 原子语义：拒绝时不返回部分结果，输入网格保持不变。
        raise ZeroVolumeError(tuple(zero))

    new_mesh = mesh.with_repaired_tets(tuple(repaired))
    # 显式不变量：顶点一个都没动（对象恒等）。
    assert all(a is b for a, b in zip(new_mesh.vertices, mesh.vertices, strict=True))
    return new_mesh
