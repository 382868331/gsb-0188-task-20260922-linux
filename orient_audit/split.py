"""内部边分裂：在有理参数位置插点并重建相邻四面体。

几何（全部精确 :class:`fractions.Fraction`）：

- 边 ``(a, b)`` 的两端点为 A、B，插点 ``P = A + t·(B - A)``，``0 < t < 1``；
- 每个含该边的四面体（“星内”单元）写成 ``(a, b, u, w)``（``u,w`` 为
  另外两个顶点，保持原序），沿切面 ``{P, u, w}`` 分裂成两个子单元::

      T_a = (P, b, u, w)   # a 槽换成 P，体积 = (1 - t) · V(T)
      T_b = (a, P, u, w)   # b 槽换成 P，体积 = t · V(T)

  （手算：det(P-a, b-a, u-a, w-a) 对 P=a+t(b-a) 展开即得上述系数。）
  两者与原单元体积同号（t、1-t 均为正），故朝向自然保持；
  材料标签由两个子单元分别继承。
- 与分裂边无关的单元原样保留。

保护语义（见 :class:`~orient_audit.errors.ProtectedFaceError`）：

- 全部边界面默认受保护——边界边 / 边界面上的边一律拒绝分裂
  （否则会把边界三角面一分为二、改变保护边界）；
- ``mesh.protected_faces`` 中的显式保护面（如材料分界面）同样拦截。

拒绝（原子回滚）：下列任一情况都只抛异常、不产生新网格，
输入 :class:`TetMesh` 是不可变对象，调用方引用保持原状：

- 输入类型 / 形状非法（:class:`~orient_audit.errors.MeshInputError`）；
- 边不存在、``t`` 不在开区间 ``(0, 1)``、星为空（不是网格的边）；
- 边落在保护面（边界或显式保护面）上；
- 边是非流形边（含该边的面被 >2 个单元引用），或边的链环是多个
  分离环（一条边穿过多张互不接触的单元片，夹伤非流形）；
- 星内存在零体积单元（``t=0/1`` 之外无法保证子体积严格非零）；
- 分裂后单元数超过 :data:`~orient_audit.exact.MAX_TETS`。

成功时除上述重建外，还用与构造相互独立的不变量（总 6 倍有向体积
精确守恒、边界面键集合逐点不变、每个切面恰为两子单元反向共享）
做结果自检，自检不过即 :class:`AssertionError`（表示实现缺陷而非
输入问题），同样不返回半成品。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .errors import EdgeSplitError, MeshInputError, ProtectedFaceError
from .exact import MAX_TETS, signed_volume6
from .mesh import FACE_LOCAL, FaceKey, TetMesh, build_face_incidence


@dataclass(frozen=True)
class EdgeSplitResult:
    """内部边分裂的结果（不可变）。

    属性:
        mesh: 分裂后的新网格（新对象；旧网格不被修改）。
        edge: 实际分裂的边 ``(a, b)``，保持入参顺序（不排序）。
        new_vertex_index: 新插点 P 的顶点下标（等于原顶点数）。
        parameter: 实际使用的参数 ``t``（精确 :class:`Fraction`）。
        star_tet_indices: 被分裂的原单元下标（升序）。
        child_tet_pairs: 每个星内单元对应的两个新单元下标
                         ``(原下标, (T_a 下标, T_b 下标))``，按原下标升序。
        untouched_tet_indices: 未参与分裂、原样保留的单元下标（升序）。
    """

    mesh: TetMesh
    edge: tuple[int, int]
    new_vertex_index: int
    parameter: Fraction
    star_tet_indices: tuple[int, ...]
    child_tet_pairs: tuple[tuple[int, tuple[int, int]], ...]
    untouched_tet_indices: tuple[int, ...]


def _parse_t(value: Any) -> Fraction:
    """解析分裂参数：int / Fraction / 有理字符串，拒绝 float、bool、非有限数。"""

    if isinstance(value, bool) or not isinstance(value, int | Fraction | str):
        raise MeshInputError(
            "分裂参数 t 必须是 int / fractions.Fraction / 有理字符串"
            f"，得到 {type(value).__name__}（不接受 float，请用 Fraction）",
            "t",
        )
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    try:
        return Fraction(value)
    except (ZeroDivisionError, ValueError, ArithmeticError):
        raise MeshInputError(f"分裂参数 {value!r} 不是有限有理数", "t")


def _edge_endpoints(value: Any, name: str, n_verts: int) -> int:
    """校验边端点：真正的 int（拒绝 bool）且不越界。"""

    if isinstance(value, bool) or not isinstance(value, int):
        raise MeshInputError(
            f"边端点 {name} 必须是 int，得到 {type(value).__name__}", name
        )
    if not 0 <= value < n_verts:
        raise MeshInputError(
            f"边端点 {name}={value} 越界（当前顶点数 {n_verts}，"
            f"合法区间 [0, {n_verts - 1}]）",
            name,
        )
    return value


def _face_keys_of_tet(tet: tuple[int, int, int, int]) -> tuple[FaceKey, ...]:
    return tuple(
        tuple(sorted((tet[local[0]], tet[local[1]], tet[local[2]])))  # type: ignore[misc]
        for local in FACE_LOCAL
    )


def split_interior_edge(
    mesh: TetMesh,
    a: int,
    b: int,
    t: int | Fraction | str = Fraction(1, 2),
) -> EdgeSplitResult:
    """在内部边 ``(a, b)`` 的参数位置 ``t`` 插点，重建相邻四面体。

    详见模块 docstring。成功返回 :class:`EdgeSplitResult`；
    任何前提不满足都抛异常且不修改 ``mesh``（不可变对象，原子语义）。
    """

    # ---------------------------------------------------------------- 输入校验
    if not isinstance(mesh, TetMesh):
        raise TypeError(
            "split_interior_edge 需要 TetMesh 实例（请先用 TetMesh.build 构造）"
        )
    n_verts = len(mesh.vertices)
    a = _edge_endpoints(a, "a", n_verts)
    b = _edge_endpoints(b, "b", n_verts)
    if a == b:
        raise MeshInputError(f"边的两个端点必须不同，得到 a == b == {a}", "b")
    tt = _parse_t(t)
    if not 0 < tt < 1:
        raise EdgeSplitError(
            f"分裂参数 t 必须在开区间 (0, 1) 内，得到 t = {tt}（端点 0/1 不插新点）",
            edge=(a, b),
            reason="t_out_of_range",
        )

    # -------------------------------- 星：收集含边 (a,b) 的单元与 (槽位, u, w)
    # star 项：(单元下标, 槽位 a, 槽位 b, 其余两顶点（按原序）)
    star: list[tuple[int, int, int, tuple[int, int]]] = []
    for ti, tet in enumerate(mesh.tets):
        slots = {v: k for k, v in enumerate(tet)}
        if a in slots and b in slots:
            ka, kb = slots[a], slots[b]
            others = tuple(v for k, v in enumerate(tet) if k not in (ka, kb))
            star.append((ti, ka, kb, others))  # type: ignore[arg-type]
    if not star:
        raise EdgeSplitError(
            f"边 ({a}, {b}) 不是网格中的边：没有任何四面体同时包含这两个顶点",
            edge=(a, b),
            reason="edge_not_found",
        )
    star.sort()
    star_indices = tuple(ti for ti, *_ in star)

    # --------------------------------------------- 星的拓扑（边的链环）
    # 每个星内单元 {a,b,u,w} 贡献：
    #   - 两个含整条分裂边的“链环面” {a,b,u}、{a,b,w}；
    #   - 两个只含单端点的“帽面” {a,u,w}、{b,u,w}（分裂后分别由
    #     T_b / T_a 原样继承，重数不变，无需额外检查）；
    #   - 一条链环边 (u, w)。
    # 统计全网格中含 a,b 的面的重数（引用该面的单元必含边 ab，故都在星内）。
    link_face_count: dict[FaceKey, int] = {}
    link_edges: list[tuple[int, int]] = []
    link_vertices: set[int] = set()
    for ti, _ka, _kb, others in star:
        u, w = others
        link_edges.append((u, w))
        link_vertices.update((u, w))
        for v in (u, w):
            key = tuple(sorted((a, b, v)))
            link_face_count[key] = link_face_count.get(key, 0) + 1  # type: ignore[index]

    # --------------------------------------------- 保护面检查（边界默认保护）
    explicit_faces = set(mesh.protected_faces)
    boundary = mesh.boundary_face_keys()
    hit_boundary = {key for key in link_face_count if key in boundary}
    hit_explicit = {key for key in link_face_count if key in explicit_faces}
    if hit_boundary or hit_explicit:
        raise ProtectedFaceError(
            edge=(a, b),
            boundary_faces=tuple(hit_boundary),
            explicit_faces=tuple(hit_explicit),
        )

    # --------------------------------------------- 非流形边：链环面重数 > 2
    over = tuple(sorted(k for k, n in link_face_count.items() if n > 2))
    if over:
        raise EdgeSplitError(
            f"边 ({a}, {b}) 是非流形边：含该边的面 {list(over)} 被超过 2 个"
            "单元引用，无法保证分裂后邻接一致，已拒绝",
            edge=(a, b),
            tet_indices=star_indices,
            reason="non_manifold_edge",
        )

    # 到这里链环面重数只能是 1 或 2；重数 1 即边界面含 ab，
    # 已在上面作为保护面拦截。故每个链环顶点度数恰为 2，链环是若干个环。
    # 还必须恰为一个环：多个环意味着同一条边穿过多张互不接触的片
    # （夹伤非流形），分裂会把该缺陷原样复制到两条新边，明确拒绝。
    degrees: dict[int, int] = {v: 0 for v in link_vertices}
    adj: dict[int, set[int]] = {v: set() for v in link_vertices}
    for u, w in link_edges:
        degrees[u] += 1
        degrees[w] += 1
        adj[u].add(w)
        adj[w].add(u)
    bad_deg = tuple(sorted(v for v, d in degrees.items() if d != 2))
    if bad_deg:
        # 理论上不该到达（重数 1 已被保护检查拦），作为防御性检查保留。
        raise EdgeSplitError(
            f"边 ({a}, {b}) 的星在顶点 {bad_deg} 处不闭合（链环不是环），"
            "该边位于开放边界或非流形处，已拒绝",
            edge=(a, b),
            tet_indices=star_indices,
            reason="open_link",
        )
    seen: set[int] = set()
    n_cycles = 0
    for start in link_vertices:
        if start in seen:
            continue
        n_cycles += 1
        stack = [start]
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v)
            stack.extend(adj[v] - seen)
    if seen != link_vertices or n_cycles != 1:
        raise EdgeSplitError(
            f"边 ({a}, {b}) 穿过多个互不接触的单元片（链环是 {n_cycles} 个"
            "分离环而非单环，夹伤非流形），分裂无法唯一定义邻接，已拒绝",
            edge=(a, b),
            tet_indices=star_indices,
            reason="disconnected_link",
        )

    # --------------------------------------------- 星内零体积（含退化星）
    zero_star: list[int] = []
    for ti, *_ in star:
        p0, p1, p2, p3 = (mesh.vertices[k] for k in mesh.tets[ti])
        if signed_volume6(p0, p1, p2, p3) == 0:
            zero_star.append(ti)
    if zero_star:
        raise EdgeSplitError(
            f"边 ({a}, {b}) 的相邻四面体中存在零体积单元 {tuple(zero_star)}，"
            "分裂会产生退化子单元，已拒绝",
            edge=(a, b),
            tet_indices=tuple(zero_star),
            reason="zero_volume_star",
        )

    # --------------------------------------------- 规模上限
    new_count = len(mesh.tets) + len(star)
    if new_count > MAX_TETS:
        raise EdgeSplitError(
            f"分裂后四面体数 {new_count}（{len(mesh.tets)} + "
            f"{len(star)} 个星内单元）超过上限 {MAX_TETS}，已拒绝",
            edge=(a, b),
            tet_indices=star_indices,
            reason="too_many_tets",
        )

    # ---------------------------------------------------------------- 构造新点
    A = mesh.vertices[a]
    B = mesh.vertices[b]
    P = tuple(A[k] + tt * (B[k] - A[k]) for k in range(3))
    p_index = n_verts
    new_vertices = mesh.vertices + (P,)  # type: ignore[assignment]

    # 槽位 / 其余两点按单元下标查表。
    star_by_tet = {ti: (ka, kb, others) for ti, ka, kb, others in star}
    star_pos = set(star_indices)
    new_tets: list[tuple[int, int, int, int]] = []
    new_labels: list[Any] = []
    untouched: list[int] = []
    # (父单元, T_a 新下标, T_b 新下标, 切面键)
    splits: list[tuple[int, int, int, FaceKey]] = []

    for ti, old_tet in enumerate(mesh.tets):
        if ti not in star_pos:
            untouched.append(len(new_tets))
            new_tets.append(old_tet)
            new_labels.append(mesh.tet_labels[ti])
            continue
        ka, kb, others = star_by_tet[ti]
        u, w = others
        label = mesh.tet_labels[ti]
        # 槽位替换：T_a 把 a 槽换成 P；T_b 把 b 槽换成 P。
        ta = tuple(p_index if k == ka else v for k, v in enumerate(old_tet))
        tb = tuple(p_index if k == kb else v for k, v in enumerate(old_tet))
        ia = len(new_tets)
        new_tets.append(ta)  # type: ignore[arg-type]
        new_labels.append(label)
        ib = len(new_tets)
        new_tets.append(tb)  # type: ignore[arg-type]
        new_labels.append(label)
        splits.append((ti, ia, ib, tuple(sorted((p_index, u, w)))))  # type: ignore[arg-type]

    # 直接构造（不经过 build 的重复检测）：星内原单元严格非零体积、
    # 0<t<1，两子单元四点集合互不相同、也不可能与其他单元重复。
    # 输入是不可变对象，在新对象构造成功前完全未被触碰。
    new_mesh = TetMesh(
        new_vertices,
        tuple(new_tets),
        tuple(new_labels),
        mesh.protected_faces,
    )

    # ---------------------------------------------------------------- 结果自检
    _verify_split_invariants(mesh, new_mesh, a, b, p_index, tt, splits)

    child_pairs = tuple((parent, (ia, ib)) for parent, ia, ib, _ in splits)
    return EdgeSplitResult(
        mesh=new_mesh,
        edge=(a, b),
        new_vertex_index=p_index,
        parameter=tt,
        star_tet_indices=star_indices,
        child_tet_pairs=child_pairs,
        untouched_tet_indices=tuple(untouched),
    )


def _verify_split_invariants(
    old: TetMesh,
    new: TetMesh,
    a: int,
    b: int,
    p_index: int,
    tt: Fraction,
    splits: list[tuple[int, int, int, FaceKey]],
) -> None:
    """构造后独立自检；任一不过即 AssertionError（不返回半成品）。"""

    # 1) 单元数：每个星内单元 1 -> 2，其余不变。
    assert len(new.tets) == len(old.tets) + len(splits), (
        f"单元数不符: {len(new.tets)} != {len(old.tets)} + {len(splits)}"
    )

    def total6(m: TetMesh) -> Fraction:
        total = Fraction(0)
        for tet in m.tets:
            p0, p1, p2, p3 = (m.vertices[k] for k in tet)
            total += signed_volume6(p0, p1, p2, p3)
        return total

    # 2) 总 6 倍有向体积精确守恒（直接用精确行列式独立重算）。
    assert total6(new) == total6(old), "分裂前后总 6 倍有向体积不守恒"

    # 3) 每个星内单元：子体积精确为 (1-t)·V（T_a）与 t·V（T_b），
    #    严格非零且与父同号；两子单元顶点恰为“父单元把对应槽位换成 P”。
    for parent, ia, ib, cut_key in splits:
        old_tet = old.tets[parent]
        V = signed_volume6(*(old.vertices[k] for k in old_tet))
        assert V != 0, f"父单元 {parent} 零体积不应进入分裂"
        va = signed_volume6(*(new.vertices[k] for k in new.tets[ia]))
        vb = signed_volume6(*(new.vertices[k] for k in new.tets[ib]))
        assert va == (1 - tt) * V, (
            f"单元 {parent} 的 T_a 体积 {va} 不等于 (1-t)·V = {(1 - tt) * V}"
        )
        assert vb == tt * V, (
            f"单元 {parent} 的 T_b 体积 {vb} 不等于 t·V = {tt * V}"
        )
        assert (V > 0) == (va > 0) == (vb > 0), (
            f"单元 {parent} 子单元朝向与父单元不一致"
        )
        # 两子单元四点集合 = 父四点集合把 a（或 b）换成 P。
        sa, sb, so = set(new.tets[ia]), set(new.tets[ib]), set(old_tet)
        assert sa == (so - {a}) | {p_index}, f"T_a 顶点集合错误: {sa}"
        assert sb == (so - {b}) | {p_index}, f"T_b 顶点集合错误: {sb}"
        # 切面 {P,u,w} 必须是两子单元的公共面。
        assert cut_key in _face_key_set(new.tets[ia]), f"T_a 不含切面 {cut_key}"
        assert cut_key in _face_key_set(new.tets[ib]), f"T_b 不含切面 {cut_key}"

    # 4) 切面在新网格中恰被两个子单元以相反诱导朝向共享（合法内面）。
    incidence = build_face_incidence(new)
    for _parent, ia, ib, cut_key in splits:
        face = incidence.get(cut_key)
        assert face is not None, f"切面 {cut_key} 在新网格中缺失"
        assert sorted(o.tet_index for o in face.occurrences) == sorted([ia, ib]), (
            f"切面 {cut_key} 引用单元 {[o.tet_index for o in face.occurrences]}，"
            f"应为 {(ia, ib)}"
        )
        assert face.is_manifold_interior, f"切面 {cut_key} 两子单元诱导朝向未相反"

    # 5) 边界面三点键集合逐点不变（保护边界在几何上确实未被改动）。
    assert new.boundary_face_keys() == old.boundary_face_keys(), (
        "分裂改变了边界面集合，违反保护边界不变量"
    )

    # 5b) 不产生新的非流形面：>2 重面键集合分裂前后逐点相同
    #     （非流形边已在前置检查拒绝；帽面被原样继承）。
    old_inc_all = build_face_incidence(old)
    old_nm = {k for k, f in old_inc_all.items() if f.multiplicity > 2}
    new_nm = {k for k, f in incidence.items() if f.multiplicity > 2}
    assert new_nm == old_nm, (
        f"分裂引入了新的非流形面: 新 {sorted(new_nm - old_nm)} / "
        f"消失 {sorted(old_nm - new_nm)}"
    )

    # 6) 旧顶点坐标对象恒等保留（只追加一个新顶点）。
    assert len(new.vertices) == len(old.vertices) + 1
    assert all(x is y for x, y in zip(new.vertices, old.vertices)), (
        "分裂移动了既有顶点"
    )


def _face_key_set(tet: tuple[int, int, int, int]) -> set[FaceKey]:
    return set(_face_keys_of_tet(tet))
