"""定向审计、面关联、非流形检测与修复规则。

所有符号判定均为精确有理运算（见 :mod:`tetmesh_audit.geometry`），
不使用任何浮点容差。
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Optional, Tuple

from .errors import DegenerateTetrahedronError, Issue
from .geometry import (
    FACE_REMOVAL_ORDER,
    chain_triangle,
    induced_orientations_opposite,
    outward_triangle,
    signed_six_volume,
)
from .mesh import TetMesh


class BoundaryTriangle:
    """边界三角面（仅被一个单元引用的面）。

    属性：
        vertices: 按该单元外法向定向的 ``(p, q, r)``（零体积单元的面
            方向无定义，此时取边界链顺序且 ``orientation_valid`` 为 False）。
        tet_index: 所属单元下标。
        local_face: 在单元中的局部面号 0..3（0 对顶点 0，依此类推）。
        component: 所属四面体连通分量编号。
        orientation_valid: 定向是否有定义（所属单元非零体积时为 True）。
    """

    __slots__ = ("vertices", "tet_index", "local_face", "component", "orientation_valid")

    def __init__(self, vertices, tet_index, local_face, component, orientation_valid):
        self.vertices = tuple(vertices)
        self.tet_index = tet_index
        self.local_face = local_face
        self.component = component
        self.orientation_valid = orientation_valid

    def __repr__(self):  # pragma: no cover
        return (
            f"BoundaryTriangle(vertices={self.vertices}, tet={self.tet_index}, "
            f"face={self.local_face}, component={self.component})"
        )


class InteriorFace:
    """恰被两个单元共享的内面。

    属性：
        vertices: 排序后的三点索引（无序面键）。
        tet_indices: 两侧单元下标。
        consistent: 两单元诱导方向是否相反（True 为定向一致）。
        induced: 两个单元各自诱导的有序三元组。
    """

    __slots__ = ("vertices", "tet_indices", "consistent", "induced")

    def __init__(self, vertices, tet_indices, consistent, induced):
        self.vertices = tuple(vertices)
        self.tet_indices = tuple(tet_indices)
        self.consistent = consistent
        self.induced = tuple(induced)

    def __repr__(self):  # pragma: no cover
        return f"InteriorFace(vertices={self.vertices}, tets={self.tet_indices}, consistent={self.consistent})"


class Defect:
    """修复规则无法消除的结构性缺陷（位置精确定位）。

    kind 取值：
        ``"non_manifold_face"``: 同一面被超过两个单元引用；
        ``"inconsistent_interior_face"``: 内面两侧诱导方向相同；
        ``"zero_volume_tet"``: 单元四点共面，体积为零。
    """

    __slots__ = ("kind", "face", "tet_indices", "detail")

    def __init__(self, kind, face, tet_indices, detail):
        self.kind = kind
        self.face = tuple(face) if face is not None else None
        self.tet_indices = tuple(tet_indices)
        self.detail = detail

    def __repr__(self):  # pragma: no cover
        return f"Defect(kind={self.kind!r}, face={self.face}, tets={self.tet_indices})"


class Component:
    """通过共享面连接的四面体连通分量。

    属性：
        index: 分量编号（按分量中最小单元下标升序分配）。
        tet_indices: 分量包含的单元下标（升序）。
        signed_six_volume: 分量内非零体积单元的 6 倍有向体积之和（精确）。
        num_boundary_triangles: 该分量拥有的边界三角面数量。
    """

    __slots__ = ("index", "tet_indices", "signed_six_volume", "num_boundary_triangles")

    def __init__(self, index, tet_indices, signed_six_volume, num_boundary_triangles):
        self.index = index
        self.tet_indices = tuple(tet_indices)
        self.signed_six_volume = signed_six_volume
        self.num_boundary_triangles = num_boundary_triangles

    def __repr__(self):  # pragma: no cover
        return f"Component(index={self.index}, tets={len(self.tet_indices)})"


class AuditResult:
    """完整审计结果。

    关键精确量（整数坐标为 ``int``，否则为 ``Fraction``）：
        * ``six_volumes[i]``: 单元 i 的 6 倍有向体积；
        * ``total_signed_six_volume``: 非零体积单元之和；
        * ``total_signed_volume``: 总有向体积（Fraction，= 前者 / 6）；
        * ``boundary_det_sum``: 边界定向三角面 det(p,q,r) 之和。
          当且仅当 ``is_consistent`` 为真（闭流形且全部定向一致）时，
          恒等式 ``boundary_det_sum == total_signed_six_volume`` 成立
          （散度定理对恒等向量场的离散形式，可作独立核对）。
    """

    __slots__ = (
        "mesh",
        "six_volumes",
        "issues",
        "negative_tets",
        "zero_volume_tets",
        "boundary_triangles",
        "interior_faces",
        "defects",
        "components",
        "total_signed_six_volume",
        "total_signed_volume",
        "boundary_det_sum",
        "is_consistent",
    )

    def __init__(self, **kw):
        for slot in self.__slots__:
            setattr(self, slot, kw[slot])

    def summary(self) -> str:
        """单行人类可读摘要。"""
        flags = []
        if self.negative_tets:
            flags.append(f"{len(self.negative_tets)} negative")
        if self.zero_volume_tets:
            flags.append(f"{len(self.zero_volume_tets)} zero-volume")
        nm = sum(1 for d in self.defects if d.kind == "non_manifold_face")
        ic = sum(1 for d in self.defects if d.kind == "inconsistent_interior_face")
        if nm:
            flags.append(f"{nm} non-manifold")
        if ic:
            flags.append(f"{ic} inconsistent-interface")
        state = "consistent" if self.is_consistent else ", ".join(flags) or "defective"
        return (
            f"tets={self.mesh.num_tetrahedra}, boundary_triangles="
            f"{len(self.boundary_triangles)}, components={len(self.components)}, "
            f"signed_volume={self.total_signed_volume} ({state})"
        )


class RepairResult:
    """修复结果（函数式：原网格不被修改）。"""

    __slots__ = ("source", "repaired_mesh", "flipped_tets", "audit_before", "audit_after")

    def __init__(self, source, repaired_mesh, flipped_tets, audit_before, audit_after):
        self.source = source
        self.repaired_mesh = repaired_mesh
        self.flipped_tets = tuple(flipped_tets)
        self.audit_before = audit_before
        self.audit_after = audit_after


def audit(mesh: TetMesh) -> AuditResult:
    """对网格执行定向与边界审计，返回 :class:`AuditResult`。

    审计是纯函数，不修改输入。检测是穷举的（不存在搜索预算问题）：
    每个面的关联单元数与诱导方向都被精确枚举。
    """
    verts = mesh.vertices
    tets = mesh.tetrahedra
    n = len(tets)

    # 1) 逐单元精确体积。
    six = []
    issues: List[Issue] = []
    negative: List[int] = []
    zero: List[int] = []
    for i, t in enumerate(tets):
        a, b, c, d = (verts[t[0]], verts[t[1]], verts[t[2]], verts[t[3]])
        s = signed_six_volume(a, b, c, d)
        six.append(s)
        if s < 0:
            negative.append(i)
            issues.append(Issue(i, "negative", s))
        elif s == 0:
            zero.append(i)
            issues.append(Issue(i, "zero_volume", s))

    # 2) 面关联：无序三点键 -> [(tet, local_face, chain_triple, outward_triple)]。
    #    chain_triple：边界链诱导方向（判两单元定向一致性，与体积符号无关）；
    #    outward_triple：几何外法向方向（边界三角面输出，负体积时取反）。
    incidence: Dict[tuple, List[Tuple[int, int, tuple, tuple]]] = {}
    for i, t in enumerate(tets):
        for lf, removal in enumerate(FACE_REMOVAL_ORDER):
            triple_idx = (t[removal[0]], t[removal[1]], t[removal[2]])
            key = tuple(sorted(triple_idx))
            incidence.setdefault(key, []).append(
                (i, lf, chain_triangle(t, lf), outward_triangle(t, lf, six[i]))
            )

    # 3) 并查集：任意共享面（含非流形面）都连接对应单元。
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)

    for key, recs in incidence.items():
        owners = [r[0] for r in recs]
        for o in owners[1:]:
            union(owners[0], o)

    comp_of: Dict[int, int] = {}
    comp_members: Dict[int, List[int]] = {}
    for i in range(n):
        root = find(i)
        if root not in comp_of:
            comp_of[root] = len(comp_of)
            comp_members[root] = []
        comp_members[root].append(i)
    tet_component = {i: comp_of[find(i)] for i in range(n)}

    # 4) 分类面：边界 / 一致内面 / 缺陷。
    boundary: List[BoundaryTriangle] = []
    interior: List[InteriorFace] = []
    defects: List[Defect] = []

    for key in sorted(incidence):
        recs = incidence[key]
        if len(recs) == 1:
            (ti, lf, _chain, out_tri), = recs
            boundary.append(
                BoundaryTriangle(
                    out_tri, ti, lf, tet_component[ti], six[ti] != 0
                )
            )
        elif len(recs) == 2:
            (t0, _, ch0, _o0), (t1, _, ch1, _o1) = recs
            opposite = induced_orientations_opposite(key, ch0, ch1)
            if six[t0] != 0 and six[t1] != 0 and not opposite:
                consistent = False
                defects.append(
                    Defect(
                        "inconsistent_interior_face",
                        key,
                        (t0, t1),
                        f"face {key} shared by tets {t0},{t1} has equal boundary "
                        f"chain orientation (induced {ch0} and {ch1}); one "
                        f"tetrahedron must be flipped",
                    )
                )
            else:
                # 零体积单元一侧方向无定义，不能称为一致。
                consistent = opposite and six[t0] != 0 and six[t1] != 0
            interior.append(
                InteriorFace(key, (t0, t1), consistent, (ch0, ch1))
            )
        else:
            owners = tuple(r[0] for r in recs)
            defects.append(
                Defect(
                    "non_manifold_face",
                    key,
                    owners,
                    f"face {key} is incident to {len(recs)} tetrahedra {owners}; "
                    f"manifold interior faces must have exactly 2",
                )
            )

    boundary.sort(key=lambda bt: (bt.tet_index, bt.local_face))
    for z in zero:
        defects.append(
            Defect(
                "zero_volume_tet",
                None,
                (z,),
                f"tetrahedron {z} has zero signed volume (four coplanar/coincident "
                f"vertices); orientation is undefined and repair is refused",
            )
        )
    defects.sort(key=lambda d: (d.kind, d.face if d.face is not None else (), d.tet_indices))

    # 5) 分量汇总。
    components: List[Component] = []
    bd_count: Dict[int, int] = {}
    for bt in boundary:
        bd_count[bt.component] = bd_count.get(bt.component, 0) + 1
    for root in sorted(comp_members, key=lambda r: comp_members[r][0]):
        cid = comp_of[root]
        members = sorted(comp_members[root])
        vol6 = sum((six[i] for i in members if six[i] != 0), 0)
        components.append(Component(cid, members, vol6, bd_count.get(cid, 0)))

    total_six = sum((s for s in six if s != 0), 0)
    boundary_det = 0
    for bt in boundary:
        if not bt.orientation_valid:
            continue
        p, q, r = verts[bt.vertices[0]], verts[bt.vertices[1]], verts[bt.vertices[2]]
        boundary_det += signed_six_volume((0, 0, 0), p, q, r)

    is_consistent = not negative and not zero and not defects

    return AuditResult(
        mesh=mesh,
        six_volumes=tuple(six),
        issues=tuple(issues),
        negative_tets=tuple(negative),
        zero_volume_tets=tuple(zero),
        boundary_triangles=tuple(boundary),
        interior_faces=tuple(interior),
        defects=tuple(defects),
        components=tuple(components),
        total_signed_six_volume=total_six,
        total_signed_volume=Fraction(total_six, 6),
        boundary_det_sum=boundary_det,
        is_consistent=is_consistent,
    )


def repair_orientation(mesh: TetMesh) -> RepairResult:
    """按固定规则修复负向单元。

    规则（题面契约）：

    * 精确行列式判号；负向单元**只交换最后两个顶点**（位置 2、3）；
    * 零体积单元方向无定义：只要存在就拒绝整个修复，抛
      :class:`~tetmesh_audit.errors.DegenerateTetrahedronError`
      （原子失败：先全部检查，再构造新网格，不存在部分修复）；
    * 顶点永不移动；重复单元在网格构造阶段已被单独拒绝。

    返回 :class:`RepairResult`（含修复后的新网格与修复前后审计）。
    非流形等结构缺陷不阻止修复：修复只处理单元符号，结构缺陷仍会出现在
    修复后审计中。
    """
    before = audit(mesh)
    if before.zero_volume_tets:
        raise DegenerateTetrahedronError(list(before.zero_volume_tets))

    flipped = list(before.negative_tets)
    new_tets: List[Tuple[int, int, int, int]] = []
    for i, t in enumerate(mesh.tetrahedra):
        if six_negative(before, i):
            new_tets.append((t[0], t[1], t[3], t[2]))
        else:
            new_tets.append(t)

    repaired = TetMesh(mesh.vertices, new_tets)
    # 不变量：修复不移动任何顶点。
    assert repaired.vertices == mesh.vertices
    after = audit(repaired)
    return RepairResult(mesh, repaired, flipped, before, after)


def six_negative(result: AuditResult, i: int) -> bool:
    return result.six_volumes[i] < 0
