#!/usr/bin/env python3
"""固定输入演示：四面体定向与边界审计 + 内部边分裂。

真实计算，无 sleep、无预录结果、无网络。三部分：

A. 正常结果：单位立方体 5 四面体剖分中故意翻转两个单元 ->
   audit 定位负向单元 -> repair_orientation 修复（只交换最后两顶点）->
   再审计：精确总体积 1、12 个外法向边界面、单一连通分量。
   总体积另用*独立*的面通量公式 V = (1/6)·Σ det(p,q,r) 核对
   （该公式不经过任何单元体积，也不调用本库核心）。

B. 实际触发的边界：三单元共享一个面 -> 报告非流形位置；
   加入共面退化单元 -> repair_orientation 抛 ZeroVolumeError 拒绝修复；
   重复单元 -> TetMesh.build 抛 DuplicateTetrahedraError。

C. 内部边分裂：外四面体重心剖分（4 子体，总体积解析值 32/3）->
   在内部边 (4,0) 的 t=1/3 处插点，3 个邻接体各一分为二（4->7 单元）->
   新点坐标、子体积 (1-t)·V / t·V、材料标签继承、4 个边界面逐点不变，
   全部用解析关系与独立面通量公式核对；再演示中点/反向边边界，
   以及边界保护面拒绝、非流形边拒绝、零体积星拒绝、t 越界与 500 上限拒绝。

运行：python demo.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import (
    DuplicateTetrahedraError,
    EdgeSplitError,
    MeshInputError,
    ProtectedFaceError,
    TetMesh,
    ZeroVolumeError,
    audit,
    repair_orientation,
    split_interior_edge,
)

# ---------------------------------------------------------------- 固定数据
VERTICES = (
    (0, 0, 0),  # 0
    (1, 0, 0),  # 1
    (1, 1, 0),  # 2
    (0, 1, 0),  # 3
    (0, 0, 1),  # 4
    (1, 0, 1),  # 5
    (1, 1, 1),  # 6
    (0, 1, 1),  # 7
)

# 标准 5 四面体剖分：中心正四面体体积 1/3，四个角各 1/6，合计 1。
TETS = (
    (0, 1, 2, 5),
    (0, 2, 3, 7),
    (0, 4, 5, 7),
    (2, 5, 6, 7),
    (0, 2, 7, 5),
)


def flip(tet):
    """奇置换：交换最后两个顶点。"""

    return (tet[0], tet[1], tet[3], tet[2])


def independent_flux_volume(vertices, boundary_triples) -> Fraction:
    """散度定理面通量公式（独立于被测核心，只用点坐标）。

    对外法向闭合三角网格：V = (1/6) * Σ det(p, q, r)。
    """

    def det3(p, q, r):
        return (
            p[0] * (q[1] * r[2] - q[2] * r[1])
            - p[1] * (q[0] * r[2] - q[2] * r[0])
            + p[2] * (q[0] * r[1] - q[1] * r[0])
        )

    total = Fraction(0)
    for i, j, k in boundary_triples:
        total += det3(vertices[i], vertices[j], vertices[k])
    return total / 6


def main() -> None:
    t0 = time.perf_counter()

    # ============ A. 正常结果：审计 -> 修复 -> 复核 ============
    print("=" * 68)
    print("A. 单位立方体（5 四面体剖分）：故意翻转单元 1 和 3")
    print("=" * 68)
    broken = [list(t) for t in TETS]
    broken[1] = list(flip(TETS[1]))
    broken[3] = list(flip(TETS[3]))
    mesh = TetMesh.build(VERTICES, broken)

    before = audit(mesh)
    print(f"单元数: {before.n_tets}")
    print(f"负向单元（精确行列式<0）: {list(before.negative_tets)}")
    print(f"零体积单元: {list(before.zero_tets)}")
    print(f"非流形面数: {len(before.non_manifold_faces)}")
    print(f"修复前总精确有向体积: {before.total_signed_volume} (= 1 - 2·1/6)")
    print(
        "定向冲突内面: "
        + ", ".join(str(f.key) for f in before.inconsistent_faces)
    )

    fixed = repair_orientation(mesh)
    after = audit(fixed)
    flux = independent_flux_volume(
        fixed.vertices, [f.ordered for f in after.boundary_faces]
    )
    print("-" * 68)
    print("repair_orientation 后：")
    print(f"  tets[1] = {fixed.tets[1]}  （只交换了最后两顶点）")
    print(f"  负向单元: {list(after.negative_tets)}")
    print(f"  全部 {len(after.interior_faces)} 个内面诱导方向相反: "
          f"{all(f.opposite for f in after.interior_faces)}")
    print(f"  边界面数: {len(after.boundary_faces)}（立方体 6 面 × 2 三角）")
    print(f"  连通分量: {after.components}")
    print(f"  总精确有向体积（行列式和/6）: {after.total_signed_volume}")
    print(f"  总体积（独立面通量公式核对）: {flux}")
    assert after.total_signed_volume == Fraction(1)
    assert flux == Fraction(1)
    print("  两种独立算法一致，且等于解析体积 1。")

    # ============ B. 实际触发的拒绝/异常边界 ============
    print()
    print("=" * 68)
    print("B. 缺陷边界：非流形报告 / 零体积拒绝 / 重复体拒绝")
    print("=" * 68)

    # B1. 三个四面体共享同一个三点面 (0,1,2) -> 非流形。
    nm_verts = (
        (0, 0, 0), (1, 1, 0), (1, 0, 1),
        (1, 0, 0), (0, 1, 1), (0, 1, 0),
    )
    nm_tets = [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5)]
    nm_mesh = TetMesh.build(nm_verts, nm_tets)
    nm_report = audit(nm_mesh)
    print("[B1] 三单元共享面：")
    for f in nm_report.non_manifold_faces:
        print(f"  非流形面 key={f.key}，被单元 {list(f.tet_indices)} 引用"
              f"（{len(f.tet_indices)} > 2）")
    assert not nm_report.is_manifold

    # B2. 加入一个共面（零体积）单元后，修复必须被拒绝。
    deg_verts = (*VERTICES, (1, 1, 0))  # 下标 8 与 0,1,2 共面于 z=0
    deg_mesh = TetMesh.build(deg_verts, [*TETS, (0, 1, 2, 8)])
    deg_report = audit(deg_mesh)
    print("[B2] 共面退化：")
    print(f"  audit 精确定位零体积单元: {list(deg_report.zero_tets)}（不静默纠正）")
    try:
        repair_orientation(deg_mesh)
    except ZeroVolumeError as exc:
        print(f"  repair_orientation 拒绝修复 -> {type(exc).__name__}:")
        print(f"    {exc}")

    # B3. 重复单元单独报错（排列不同也算重复）。
    print("[B3] 重复单元：")
    try:
        perm = (TETS[0][2], TETS[0][0], TETS[0][3], TETS[0][1])
        TetMesh.build(VERTICES, [*TETS, perm])
    except DuplicateTetrahedraError as exc:
        print(f"  TetMesh.build 拒绝 -> {type(exc).__name__}:")
        key, positions = exc.duplicates[0]
        print(f"    单元 {positions} 顶点集合均为 {key}")

    # B4. 非法输入可定位（顺带演示错误语义，非核心几何边界）。
    print("[B4] 非法输入：")
    for label, build in (
        ("float 坐标", lambda: TetMesh.build([(0, 0, 0.5)], [])),
        ("bool 下标", lambda: TetMesh.build(VERTICES, [(0, 1, 2, True)])),
    ):
        try:
            build()
        except MeshInputError as exc:
            print(f"  {label} -> {type(exc).__name__} 于 {exc.location}: 已拒绝")

    # ============ C. 内部边分裂：正常结果 / 边界 / 实际错误 ============
    print()
    print("=" * 68)
    print("C. 内部边分裂：重心剖分上有理插点重建")
    print("=" * 68)

    # 外四面体 0..3（底面积 8、高 4，解析体积 32/3），严格内点 4
    # 连四个面，得到 4 个等体积子体（各 8/3）。边 (4,0) 是内部边，
    # 它的星含单元 1、2、3，链环是与顶点 0 相对的三角形面。
    bverts = (
        (0, 0, 0),  # 0
        (4, 0, 0),  # 1
        (0, 4, 0),  # 2
        (0, 0, 4),  # 3
        (1, 1, 1),  # 4 —— 严格内点
    )
    btets = (
        (4, 1, 2, 3),
        (4, 0, 3, 2),
        (4, 0, 1, 3),
        (4, 0, 2, 1),
    )
    blabels = ("copper", "rubber", "rubber", "glass")
    bmesh = TetMesh.build(bverts, btets, tet_labels=blabels)
    brep = audit(bmesh)
    print(f"[C0] 重心剖分：{brep.n_tets} 单元，总精确体积 {brep.total_signed_volume}"
          "（解析值 32/3）")
    print(f"     边界面 {len(brep.boundary_faces)} 个，内面 "
          f"{len(brep.interior_faces)} 个，材料标签 {bmesh.tet_labels}")

    # C1. 正常分裂：内部边 (4,0)，t = 1/3。
    t = Fraction(1, 3)
    result = split_interior_edge(bmesh, 4, 0, t)
    smesh = result.mesh
    srep = audit(smesh)
    P = smesh.vertices[result.new_vertex_index]
    print(f"[C1] split_interior_edge(bmesh, 4, 0, t={t})：")
    print(f"  星内单元 {list(result.star_tet_indices)} 各一分为二，"
          f"单元数 {brep.n_tets} -> {srep.n_tets}，新点下标 {result.new_vertex_index}")
    print(f"  新点坐标 P = A + t(B-A) = {tuple(str(x) for x in P)}"
          "（手算 = (2/3, 2/3, 2/3)）")
    print(f"  审计：负向 {list(srep.negative_tets)}，零体积 "
          f"{list(srep.zero_tets)}，非流形面 {len(srep.non_manifold_faces)}，"
          f"内面 {len(srep.interior_faces)} 个且全部反向: "
          f"{all(f.opposite for f in srep.interior_faces)}")
    print(f"  材料标签继承：{smesh.tet_labels}")

    # 解析核对：子体积 (1-t)·V6=32/3 与 t·V6=16/3（父 V6=16）；
    # 用独立面通量公式核对分裂后边界体积仍为 32/3。
    assert srep.total_signed_volume6 == Fraction(64), "总体积必须守恒"
    flux_after = independent_flux_volume(
        smesh.vertices, [f.ordered for f in srep.boundary_faces]
    )
    assert flux_after == Fraction(32, 3)
    assert smesh.boundary_face_keys() == bmesh.boundary_face_keys()
    parent, (ia, ib) = result.child_tet_pairs[0]
    print(f"  抽查单元 {parent}：子单元 {ia} 体积 6 倍 = "
          f"{_v6(smesh, ia)} = (1-t)·16 = 32/3；"
          f"{ib} = {_v6(smesh, ib)} = t·16 = 16/3")
    print(f"  总体积 {srep.total_signed_volume}（核心行列式和）== "
          f"{flux_after}（独立面通量），边界面三点键逐点不变。")

    # C2. 边界：t=1/2 中点；反向边 (0,4) 配 t=2/3 必须得到同一个点。
    mid = split_interior_edge(bmesh, 4, 0)
    rev = split_interior_edge(bmesh, 0, 4, Fraction(2, 3))
    mid_p = tuple(str(x) for x in mid.mesh.vertices[mid.new_vertex_index])
    print(f"[C2] 边界参数：默认 t=1/2 新点 {mid_p}；"
          f"反向边 (0,4),t=2/3 新点 {tuple(str(x) for x in rev.mesh.vertices[rev.new_vertex_index])}"
          "（与 C1 同一点）")
    assert mid.parameter == Fraction(1, 2)
    assert mid.mesh.vertices[mid.new_vertex_index] == (Fraction(1, 2),) * 3
    assert rev.mesh.vertices[rev.new_vertex_index] == P

    # C3. 实际错误（每个都真的抛出，且原网格不被修改）。
    print("[C3] 拒绝（均实际触发，原网格保持原状）：")

    def show(title, call):
        try:
            call()
            print(f"  {title}: 未拒绝（异常！）")
        except (EdgeSplitError, ProtectedFaceError, MeshInputError) as exc:
            extra = ""
            if isinstance(exc, ProtectedFaceError):
                extra = f"，受冒犯面 {list(exc.boundary_faces) + list(exc.explicit_faces)}"
            print(f"  {title} -> {type(exc).__name__}"
                  f"（reason={getattr(exc, 'reason', '-')}{extra}）")

    show("分裂立方体边界边 (0,1)",
         lambda: split_interior_edge(TetMesh.build(VERTICES, TETS), 0, 1))
    # 显式保护面：把内面 (0,3,4) 声明为材料分界面，它含内部边 (4,0)。
    prot_mesh = TetMesh.build(bverts, btets, protected_faces=[(0, 3, 4)])
    show("内部边落在显式保护面（材料分界面）",
         lambda: split_interior_edge(prot_mesh, 4, 0))
    show("非流形边（5 体共边、链环含三角形）",
         lambda: split_interior_edge(
             TetMesh.build(_PINWHEEL_V, _PINWHEEL_T), 0, 1))
    show("零体积星（一个邻接体四点共面）",
         lambda: split_interior_edge(
             TetMesh.build(_DEGEN_V, _DEGEN_T), 0, 1))
    show("参数 t=0（开区间外）", lambda: split_interior_edge(bmesh, 4, 0, 0))
    show("参数 t=0.5（float 拒绝）", lambda: split_interior_edge(bmesh, 4, 0, 0.5))
    # 边不存在：4 环夹具中链环顶点 2 与 4 从不同属一个单元。
    r4v, r4t = _ring_star(4)
    show("边 (2,4) 不是网格中的边",
         lambda: split_interior_edge(TetMesh.build(r4v, r4t), 2, 4))
    # 500 上限：251 环 -> 502 单元。
    rv, rt = _ring_star(251)
    show("分裂后 502 > 500 上限",
         lambda: split_interior_edge(TetMesh.build(rv, rt), 0, 1))

    # 原子性复核：全部失败尝试之后，bmesh 审计结果与最初一致。
    assert audit(bmesh).total_signed_volume6 == Fraction(64)
    assert len(bmesh.vertices) == 5 and len(bmesh.tets) == 4
    print("  全部拒绝后 bmesh 仍是 5 顶点 / 4 单元、总体积 64/6（原子回滚）。")

    elapsed = time.perf_counter() - t0
    print()
    print(f"演示完成，实际耗时 {elapsed:.3f} 秒。")


def _v6(mesh, ti):
    """直接用精确行列式求某单元 6 倍有向体积（演示核对用，不调审计）。"""

    from orient_audit.exact import signed_volume6

    tet = mesh.tets[ti]
    p0, p1, p2, p3 = (mesh.vertices[k] for k in tet)
    return signed_volume6(p0, p1, p2, p3)


def _ring_star(k):
    """端点 0,1 的内部边，星为 k 边形环（确定性小整数坐标）。"""

    verts = [(0, 0, 0), (0, 0, 1)] + [(1, j, Fraction(1, 2)) for j in range(k)]
    tets = [(0, 1, 2 + j, 2 + ((j + 1) % k)) for j in range(k)]
    return verts, tets


_PINWHEEL_V = (
    (0, 0, 0), (0, 0, 2),
    (2, 0, 1), (1, 2, 1), (-2, 1, 1), (0, -2, 1),
)
_PINWHEEL_T = (
    (0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5),
    (0, 1, 3, 4), (0, 1, 3, 5),
)
_DEGEN_V = ((0, 0, 0), (0, 0, 4), (1, 0, 1), (-2, 0, 1), (0, 3, 1))
_DEGEN_T = ((0, 1, 2, 3), (0, 1, 3, 4), (0, 1, 4, 2))


if __name__ == "__main__":
    main()
