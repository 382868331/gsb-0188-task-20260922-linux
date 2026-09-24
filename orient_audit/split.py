"""内部边分裂：在指定有理位置插点并重建相邻四面体。

语义（详见 README）：

- :func:`split_internal_edge(mesh, edge, t)` 在内部边
  ``edge = (a, b)`` 的精确有理位置 ``c = (1-t)·a + t·b``（``0 < t < 1``）
  插入新顶点；该边星形（star，即所有含该边的四面体）中的每个单元，
  在其顶点序的*槽位*上把边端点替换成新点，得到两个子单元（体积是被
  替换顶点的仿射函数，在另一端点处取值为零）::

      A（近 a）：父序中 b 所在槽位换成 c → 有向体积 = t       · 父单元
      B（近 b）：父序中 a 所在槽位换成 c → 有向体积 = (1 - t) · 父单元

  因此两个子单元严格继承父单元的定向符号与材料标签，且逐单元体积守恒
  （全精确 :class:`fractions.Fraction`，无浮点、无容差）。
- 只允许**内部边**：边上的每个侧面（含边两端点的面）必须恰被两个单元
  引用；边界面、非流形侧面、星形邻接环不闭合、星内零体积单元、插入点
  与现有顶点重合等情况一律拒绝。
- **保护面**：边不得落在任何保护面（``mesh.protected_faces``）上，
  否则拒绝；分裂后保护面的三点集合与引用重数保持不变。
- 操作是**原子**的：:class:`~orient_audit.mesh.TetMesh` 不可变，
  分裂在新对象上构造；构造后再做体积、朝向、邻接、边界与保护面守恒
  自检，任一不通过即抛 :class:`EdgeSplitError`，原网格不会有任何
  半更新（函数要么返回完整 :class:`SplitResult` 要么抛异常）。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .errors import EdgeSplitError, MeshInputError
from .exact import MAX_TETS, parse_index, signed_volume6
from .mesh import FACE_LOCAL, TetMesh, build_face_incidence


@dataclass(frozen=True)
class SplitResult:
    """一次内部边分裂的完整结果（不可变）。

    属性:
        mesh: 分裂后的新网格（旧网格不被修改；旧顶点坐标对象恒等保留）。
        edge: 被分裂的边 ``(a, b)``。
        parameter: 精确有理参数 ``t``，新点为 ``(1-t)·a + t·b``。
        new_vertex_index: 新顶点在 ``mesh.vertices`` 中的下标。
        star_tet_indices: 星形单元（含该边的旧单元）在**原网格**中的下标。
        child_tet_indices: 与 ``star_tet_indices`` 对齐的
            ``(A 下标, B 下标)``，二者都是**新网格**中的单元下标。
    """

    mesh: TetMesh
    edge: tuple[int, int]
    parameter: Fraction
    new_vertex_index: int
    star_tet_indices: tuple[int, ...]
    child_tet_indices: tuple[tuple[int, int], ...]


def _parse_parameter(value: object) -> Fraction:
    """分裂参数 t：int / Fraction / 有理字符串，拒绝 float、bool。"""

    if isinstance(value, bool) or not isinstance(value, int | Fraction | str):
        raise MeshInputError(
            "分裂参数 t 必须是 int / fractions.Fraction / 有理字符串，"
            f"得到 {type(value).__name__}（不接受 float，请用 Fraction）",
            "t",
        )
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    try:
        return Fraction(value)
    except (ZeroDivisionError, ValueError, ArithmeticError):
        raise MeshInputError(
            f"分裂参数 {value!r} 不是有限有理数（拒绝 NaN/Infinity）", "t"
        )


def _parse_edge(edge: object, n_verts: int) -> tuple[int, int]:
    """校验边：恰 2 个互不相同的合法 int 下标（拒绝 bool）。"""

    if isinstance(edge, str) or not isinstance(edge, (list, tuple)):
        raise MeshInputError(
            "edge 必须是两个顶点下标组成的序列 (a, b)", "edge"
        )
    if len(edge) != 2:
        raise MeshInputError(
            f"edge 必须恰好含 2 个顶点下标，得到 {len(edge)} 个", "edge"
        )
    a = parse_index(edge[0], "edge[0]", n_verts)
    b = parse_index(edge[1], "edge[1]", n_verts)
    if a == b:
        raise MeshInputError(f"边的两个端点必须不同，得到 ({a}, {b})", "edge")
    return a, b


def split_internal_edge(mesh: TetMesh, edge, t) -> SplitResult:
    """在内部边 ``edge=(a, b)`` 的有理位置 ``t`` 插点并分裂相邻四面体。

    参数:
        mesh: 已构造的 :class:`TetMesh`（不会被修改）。
        edge: 两个不同顶点下标组成的序列。
        t: 精确有理参数，``0 < t < 1``；新点
           ``c = (1-t)·verts[a] + t·verts[b]``。接受 ``int``、
           :class:`~fractions.Fraction` 与有理字符串，拒绝 ``float``。

    返回:
        :class:`SplitResult`。

    异常:
        MeshInputError: 参数形状 / 类型 / 范围非法（可定位到
            ``edge[...]`` 或 ``t``）。
        EdgeSplitError: 操作会破坏拓扑或保护边界、星内有零体积单元、
            新点与旧顶点重合、结果超 :data:`~orient_audit.exact.MAX_TETS`
            或分裂后自检失败；此时原网格保持不变。
    """

    if not isinstance(mesh, TetMesh):
        raise TypeError(
            "split_internal_edge 需要 TetMesh 实例（请先用 TetMesh.build 构造）"
        )

    a, b = _parse_edge(edge, len(mesh.vertices))
    param = _parse_parameter(t)
    if not 0 < param < 1:
        raise MeshInputError(
            f"分裂参数 t 必须严格位于开区间 (0, 1)，得到 {param}", "t"
        )

    tets = mesh.tets
    verts = mesh.vertices
    incidence = build_face_incidence(mesh)
    multiplicity = {key: f.multiplicity for key, f in incidence.items()}
    boundary_keys = {key for key, m in multiplicity.items() if m == 1}
    protected = set(mesh.protected_faces)

    # ---- 星形：所有同时含 a、b 的单元 ----
    star = tuple(ti for ti, tet in enumerate(tets) if a in tet and b in tet)
    if not star:
        raise EdgeSplitError(
            "edge_not_found",
            f"边 ({a}, {b}) 不属于任何四面体，无法分裂",
            edge=(a, b),
        )

    # 每个星单元含边两端点的两个侧面（link faces）。
    link_keys: set[tuple[int, int, int]] = set()
    for ti in star:
        tet = tets[ti]
        for local in FACE_LOCAL:
            face = (tet[local[0]], tet[local[1]], tet[local[2]])
            key = tuple(sorted(face))
            if a in key and b in key:
                link_keys.add(key)  # type: ignore[arg-type]

    # ---- 保护面：边不得落在保护面上 ----
    touched_protected = tuple(sorted(k for k in link_keys if k in protected))
    if touched_protected:
        raise EdgeSplitError(
            "edge_on_protected_face",
            f"边 ({a}, {b}) 落在保护面 {list(touched_protected)} 上，"
            "分裂会破坏保护边界，已拒绝（网格未改动）",
            edge=(a, b),
            faces=touched_protected,
        )

    # ---- 内部边：任何侧面都不能是边界面 ----
    boundary_links = tuple(sorted(k for k in link_keys if k in boundary_keys))
    if boundary_links:
        raise EdgeSplitError(
            "edge_on_boundary",
            f"边 ({a}, {b}) 的侧面 {list(boundary_links)} 是边界面，"
            "只能分裂内部边，已拒绝（网格未改动）",
            edge=(a, b),
            faces=boundary_links,
        )

    # ---- 流形：每个侧面必须恰被两个单元引用 ----
    bad_links = tuple(
        sorted(k for k in link_keys if multiplicity.get(k, 0) != 2)
    )
    if bad_links:
        raise EdgeSplitError(
            "non_manifold_edge",
            f"边 ({a}, {b}) 的侧面 {list(bad_links)} 引用数 != 2，"
            "星形非流形，已拒绝（网格未改动）",
            edge=(a, b),
            faces=bad_links,
        )

    # ---- 邻接环：星单元经侧面构成的图必须是单一闭合环 ----
    # 每个侧面恰关联两个星单元；每个星单元有两个侧面 → 各节点度数 2。
    neighbors: dict[int, set[int]] = {ti: set() for ti in star}
    for key in link_keys:
        owners = tuple(
            o.tet_index for o in incidence[key].occurrences
            if o.tet_index in neighbors
        )
        if len(owners) == 2:
            x, y = owners
            neighbors[x].add(y)
            neighbors[y].add(x)
    # 度数 2 且连通 ⇔ 单一环；不连通则是多个小环（几何上不可能由一条
    # 边的星构成，出现即说明输入拓扑异常）。
    if any(len(ns) != 2 for ns in neighbors.values()) or not _ring_connected(
        star, neighbors
    ):
        raise EdgeSplitError(
            "ring_disconnected",
            f"边 ({a}, {b}) 的 {len(star)} 个邻接单元不构成单一闭合邻接环，"
            "已拒绝（网格未改动）",
            edge=(a, b),
        )

    # ---- 星内不得有零体积单元（交换/分裂都无法消除零行列式）----
    zero_star = tuple(
        ti
        for ti in star
        if signed_volume6(*(verts[k] for k in tets[ti])) == 0
    )
    if zero_star:
        raise EdgeSplitError(
            "zero_volume_star",
            f"边 ({a}, {b}) 的星形中存在零体积单元 {list(zero_star)}，"
            "分裂会产生零体积子单元，已拒绝（网格未改动）",
            edge=(a, b),
        )

    # ---- 新点坐标（精确）与重合检查 ----
    pa, pb = verts[a], verts[b]
    new_point = tuple(
        (Fraction(1) - param) * pa[k] + param * pb[k] for k in range(3)
    )
    colliding = tuple(
        vi for vi, p in enumerate(verts) if p == new_point
    )
    if colliding:
        raise EdgeSplitError(
            "point_collides",
            f"t = {param} 处插入点与现有顶点 {list(colliding)} 坐标重合，"
            "已拒绝（网格未改动）",
            edge=(a, b),
        )

    # ---- 上限：n_tets + |star| 不得超过 MAX_TETS ----
    predicted = len(tets) + len(star)
    if predicted > MAX_TETS:
        raise EdgeSplitError(
            "too_many_tets",
            f"分裂将产生 {predicted} 个四面体，超过上限 {MAX_TETS}"
            f"（当前 {len(tets)} 个，星形 {len(star)} 个单元各一分为二）",
            edge=(a, b),
        )

    # ---- 构造新顶点 / 新单元（先在普通列表里完成，不触碰旧网格）----
    #
    # 子单元顶点序在*父单元的槽位*上替换得到（体积是被替换顶点的仿射
    # 函数，另一端点处取值为零）：
    #   A（近 a）：父序中槽位 j 的 b 换成 c → 有向体积 = t 份；
    #   B（近 b）：父序中槽位 i 的 a 换成 c → 有向体积 = (1-t) 份。
    # 因此无论 i,j 落在哪两个槽，子单元都严格继承父单元的定向符号。
    new_index = len(verts)
    new_vertices = (*verts, new_point)
    new_tets: list[tuple[int, int, int, int]] = []
    expected_materials: list[int] = []
    expected_vol6: list[Fraction] = []
    star_set = set(star)
    child_indices: dict[int, tuple[int, int]] = {}
    for ti, tet in enumerate(tets):
        parent_vol = signed_volume6(*(verts[k] for k in tet))
        if ti not in star_set:
            new_tets.append(tet)
            expected_materials.append(mesh.materials[ti])
            expected_vol6.append(parent_vol)
            continue
        i = tet.index(a)
        j = tet.index(b)
        near_a = tuple(new_index if q == j else tet[q] for q in range(4))
        near_b = tuple(new_index if q == i else tet[q] for q in range(4))
        idx_a = len(new_tets)
        new_tets.append(near_a)  # type: ignore[arg-type]
        expected_materials.append(mesh.materials[ti])
        expected_vol6.append(param * parent_vol)
        idx_b = len(new_tets)
        new_tets.append(near_b)  # type: ignore[arg-type]
        expected_materials.append(mesh.materials[ti])
        expected_vol6.append((Fraction(1) - param) * parent_vol)
        child_indices[ti] = (idx_a, idx_b)

    # build 做完整结构校验（下标、重复体、标签、保护面）；它会把坐标
    # 重新解析成新 Fraction，因此最终对象沿用其校验后的 tets/标签，而
    # 顶点直接挂回旧对象（坐标与解析结果逐一相等），保证旧顶点恒等。
    validated = TetMesh.build(
        new_vertices,
        new_tets,
        materials=tuple(expected_materials),
        protected_faces=mesh.protected_faces,
    )
    new_mesh = TetMesh(
        (*verts, new_point),
        validated.tets,
        validated.materials,
        validated.protected_faces,
    )

    # ---- 分裂后自检：体积 / 朝向 / 邻接 / 边界 / 保护面 / 标签守恒 ----
    _verify_split(
        mesh,
        new_mesh,
        a,
        b,
        new_index,
        star,
        child_indices,
        param,
        link_keys,
        tuple(expected_vol6),
        tuple(expected_materials),
    )

    # 旧顶点坐标对象恒等保留（不移动顶点）。
    assert all(
        x is y for x, y in zip(new_mesh.vertices[:new_index], verts, strict=True)
    )

    return SplitResult(
        mesh=new_mesh,
        edge=(a, b),
        parameter=param,
        new_vertex_index=new_index,
        star_tet_indices=star,
        child_tet_indices=tuple(child_indices[ti] for ti in star),
    )


def _ring_connected(
    star: tuple[int, ...], neighbors: dict[int, set[int]]
) -> bool:
    """星单元经侧面邻接是否连通（各节点度数已为 2 时即单一环）。"""

    seen = {star[0]}
    stack = [star[0]]
    while stack:
        cur = stack.pop()
        for nxt in neighbors[cur]:
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return len(seen) == len(star)


def _verify_split(
    old: TetMesh,
    new: TetMesh,
    a: int,
    b: int,
    c: int,
    star: tuple[int, ...],
    child_indices: dict[int, tuple[int, int]],
    param: Fraction,
    link_keys: set[tuple[int, int, int]],
    expected_vol6: tuple[Fraction, ...],
    expected_materials: tuple[int, ...],
) -> None:
    """构造后守恒自检；任何不通过都按原子语义拒绝（旧网格未被改动）。

    期望值 ``expected_vol6`` / ``expected_materials`` 由分裂方按仿射
    关系与继承规则预算；这里重新从结果网格直接计算行列式与面关联来
    比对，不信任构造过程本身。
    """

    def fail(detail: str) -> None:
        raise EdgeSplitError(
            "verification_failed",
            f"边 ({a}, {b}) 分裂后自检未通过：{detail}（已丢弃新网格，"
            "原网格保持不变）",
            edge=(a, b),
        )

    # 1) 从结果网格直接重算每个单元的精确 6 倍有向体积，与预算逐单元比对。
    new_vol = [
        signed_volume6(*(new.vertices[k] for k in tet)) for tet in new.tets
    ]
    if len(new_vol) != len(expected_vol6):
        fail("结果单元数与预期不符")
    for qi, (got, want) in enumerate(zip(new_vol, expected_vol6, strict=True)):
        if got != want:
            fail(f"新单元 {qi} 的 6 倍有向体积 {got} != 预期 {want}")
    old_total = sum(
        (signed_volume6(*(old.vertices[k] for k in tet)) for tet in old.tets),
        Fraction(0),
    )
    if sum(new_vol, Fraction(0)) != old_total:
        fail("分裂前后总 6 倍有向体积不守恒")

    # 2) 朝向继承：每个子单元体积符号必须与父单元一致且严格非零。
    for ti in star:
        parent_vol = signed_volume6(*(old.vertices[k] for k in old.tets[ti]))
        ia, ib = child_indices[ti]
        for qi in (ia, ib):
            if new_vol[qi] == 0 or (new_vol[qi] > 0) != (parent_vol > 0):
                fail(f"单元 {ti} 的子单元 {qi} 未继承父单元朝向")

    # 3) 面关联：旧侧面消失且各替换为两个重数 2 的新内面；每个父单元
    #    新增的 (c,k,l) 恰被其两个子单元共享；其余面（含边界面、盖面、
    #    保护面）键与重数完全不变。
    old_inc = build_face_incidence(old)
    new_inc = build_face_incidence(new)
    old_mult = {k: f.multiplicity for k, f in old_inc.items()}
    new_mult = {k: f.multiplicity for k, f in new_inc.items()}

    expected_mult: dict[tuple[int, int, int], int] = dict(old_mult)
    for key in link_keys:
        if expected_mult.pop(key, None) is None:
            fail(f"旧侧面 {key} 缺失")
        x = next(v for v in key if v not in (a, b))
        ka = tuple(sorted((a, c, x)))
        kb = tuple(sorted((b, c, x)))
        expected_mult[ka] = expected_mult.get(ka, 0) + 2
        expected_mult[kb] = expected_mult.get(kb, 0) + 2
    for ti in star:
        others = tuple(v for v in old.tets[ti] if v not in (a, b))
        key = tuple(sorted((c, others[0], others[1])))
        expected_mult[key] = expected_mult.get(key, 0) + 2
    if expected_mult != new_mult:
        fail("分裂后面关联（面键与引用重数）与预期不符")

    # 4) 朝向一致性保持：父单元间共享的 (c,k,l) 是全新内面，必须诱导
    #    相反；旧侧面的两个“半边”必须继承该侧面分裂前的一致/冲突状态
    #    （分裂既不制造也不修复定向冲突）。
    for ti in star:
        others = tuple(v for v in old.tets[ti] if v not in (a, b))
        key = tuple(sorted((c, others[0], others[1])))
        if not new_inc[key].is_manifold_interior:
            fail(f"父子共享新内面 {key} 诱导朝向相同")
    for key in link_keys:
        was = old_inc[key].is_manifold_interior
        x = next(v for v in key if v not in (a, b))
        for half in (tuple(sorted((a, c, x))), tuple(sorted((b, c, x)))):
            if new_inc[half].is_manifold_interior != was:
                fail(f"侧面 {key} 的半边 {half} 未继承定向一致性状态")

    # 5) 边界有序面多重集不变（内部边分裂不得改动边界几何与朝向）。
    old_boundary = sorted(
        f.occurrences[0].ordered
        for f in old_inc.values()
        if f.multiplicity == 1
    )
    new_boundary = sorted(
        f.occurrences[0].ordered
        for f in new_inc.values()
        if f.multiplicity == 1
    )
    if old_boundary != new_boundary:
        fail("边界面有序三重集合在分裂后发生变化")

    # 6) 保护面：键与引用重数逐一不变（新点不允许落在保护面内）。
    for key in old.protected_faces:
        if new_mult.get(key) != old_mult.get(key):
            fail(f"保护面 {key} 的引用关系发生变化")

    # 7) 材料标签：新网格逐单元标签必须等于继承/保留预算。
    if tuple(new.materials) != expected_materials:
        fail("材料标签与继承预期不符")
