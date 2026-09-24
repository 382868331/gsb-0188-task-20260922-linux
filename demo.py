#!/usr/bin/env python3
"""固定输入演示：四面体定向审计、修复与内部边分裂。

真实计算，无 sleep、无预录结果、无网络。三部分：

A. 正常结果：单位立方体 5 四面体剖分中故意翻转两个单元 ->
   audit 定位负向单元 -> repair_orientation 修复（只交换最后两顶点）->
   再审计：精确总体积 1、12 个外法向边界面、单一连通分量。
   总体积另用*独立*的面通量公式 V = (1/6)·Σ det(p,q,r) 核对
   （该公式不经过任何单元体积，也不调用本库核心）。

B. 内部边分裂：Kuhn 6 剖分的体对角线 (0,6) 在 t=1/3 处插点分裂，
   6 -> 12 单元，子单元精确继承 1/3、2/3 体积份额与材料标签，
   保护面与边界不动；演示端点顺序对称、连续分裂，以及边界边 /
   保护面 / 非法参数的实际拒绝（拒绝后网格保持原状）。

C. 实际触发的边界：三单元共享一个面 -> 报告非流形位置；
   加入共面退化单元 -> repair_orientation 抛 ZeroVolumeError 拒绝修复；
   重复单元 -> TetMesh.build 抛 DuplicateTetrahedraError。

运行：python demo.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import (
    DuplicateTetrahedraError,
    EdgeSplitError,
    MeshInputError,
    TetMesh,
    ZeroVolumeError,
    audit,
    repair_orientation,
    split_internal_edge,
)
from tests.reference import ref_signed_volume6

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

# Kuhn（Freudenthal）6 四面体剖分：6 个单元共享体对角线 (0,6)，
# 每个体积 1/6；边 (0,6) 是网格中唯一的内部边（星形为单一闭合 6 环）。
KUHN_TETS = (
    (0, 6, 1, 2),
    (0, 6, 2, 3),
    (0, 6, 3, 7),
    (0, 6, 7, 4),
    (0, 6, 4, 5),
    (0, 6, 5, 1),
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

    # ============ B. 内部边分裂：正常 / 边界 / 拒绝 ============
    print()
    print("=" * 68)
    print("B. 内部边分裂（Kuhn 6 剖分的体对角线 (0,6)，t = 1/3）")
    print("=" * 68)

    materials = (1, 2, 2, 3, 3, 1)
    kuhn = TetMesh.build(
        VERTICES, KUHN_TETS,
        materials=materials,
        protected_faces=[(0, 1, 2)],   # 保护一个边界面（与内边不相交）
    )
    k_before = audit(kuhn)
    print(f"分裂前：{k_before.n_tets} 个单元，总精确体积 "
          f"{k_before.total_signed_volume}，边界面 {len(k_before.boundary_faces)} 个")

    result = split_internal_edge(kuhn, (0, 6), Fraction(1, 3))
    split_mesh = result.mesh
    k_after = audit(split_mesh)
    print(f"新顶点下标 {result.new_vertex_index}，坐标 "
          f"{split_mesh.vertices[result.new_vertex_index]}（精确 1/3 处）")
    print(f"分裂后：{k_after.n_tets} 个单元，负向 {list(k_after.negative_tets)}，"
          f"边界面 {len(k_after.boundary_faces)} 个（不变），"
          f"内面 {len(k_before.interior_faces)} -> {len(k_after.interior_faces)}")

    # 手算核验：每个父单元 6 倍体积 1 → 子单元 1/3 与 2/3（独立 Bareiss）。
    ok_shares = True
    for ti, (ia, ib) in enumerate(result.child_tet_indices):
        va = ref_signed_volume6(split_mesh.vertices, split_mesh.tets[ia])
        vb = ref_signed_volume6(split_mesh.vertices, split_mesh.tets[ib])
        ok_shares &= (va, vb) == (Fraction(1, 3), Fraction(2, 3))
    print(f"6 个父单元均精确分成 6 倍体积 1/3 与 2/3"
          f"（独立 Bareiss 逐单元核对）: {ok_shares}")
    assert ok_shares
    assert k_after.total_signed_volume == Fraction(1)
    print(f"材料标签继承: {split_mesh.materials}")
    assert split_mesh.materials == (1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1)
    print(f"保护面保留: {split_mesh.protected_faces}")
    flux2 = independent_flux_volume(
        split_mesh.vertices, [f.ordered for f in k_after.boundary_faces]
    )
    print(f"总体积（独立面通量公式核对）: {flux2}；原网格仍为 "
          f"{audit(kuhn).n_tets} 个单元（原子：旧网格未改动）")
    assert flux2 == Fraction(1)

    # B2. 边界：端点顺序与参数方向对称 —— (0,6)@3/4 与 (6,0)@1/4 同点同剖分。
    fwd = split_internal_edge(kuhn, (0, 6), Fraction(3, 4)).mesh
    rev = split_internal_edge(kuhn, (6, 0), Fraction(1, 4)).mesh
    symmetric = (
        fwd.vertices == rev.vertices
        and sorted(tuple(sorted(q)) for q in fwd.tets)
        == sorted(tuple(sorted(q)) for q in rev.tets)
    )
    print(f"[B2] 边界对称性：(0,6)@3/4 与 (6,0)@1/4 几何等价: {symmetric}")
    assert symmetric
    # 连续两次分裂（t=1/2 再 t=1/2）恰到达 1/4 处。
    half = split_internal_edge(kuhn, (0, 6), Fraction(1, 2)).mesh
    quarter = split_internal_edge(half, (0, 8), Fraction(1, 2)).mesh
    print(f"[B2] 连续分裂：单元 6 -> {audit(half).n_tets} -> "
          f"{audit(quarter).n_tets}，新点 {quarter.vertices[9]}，"
          f"总体积仍为 {audit(quarter).total_signed_volume}")
    assert quarter.vertices[9] == (Fraction(1, 4),) * 3

    # B3. 实际触发的拒绝（网格均保持原状）。
    print("[B3] 分裂拒绝：")
    plain_kuhn = TetMesh.build(VERTICES, KUHN_TETS)
    rejections = (
        ("边界边 (0,1)",
         lambda: split_internal_edge(plain_kuhn, (0, 1), Fraction(1, 2))),
        ("参数 t=0（不在开区间）",
         lambda: split_internal_edge(kuhn, (0, 6), 0)),
        ("float 参数", lambda: split_internal_edge(kuhn, (0, 6), 0.5)),
    )
    # 保护面拒绝：把内侧面 (0,1,6) 列为保护面。
    protected_kuhn = TetMesh.build(
        VERTICES, KUHN_TETS, protected_faces=[(0, 1, 6)]
    )
    rejections += (
        ("边落在保护面 (0,1,6)",
         lambda: split_internal_edge(protected_kuhn, (0, 6), Fraction(1, 2))),
    )
    for label, action in rejections:
        try:
            action()
        except EdgeSplitError as exc:
            print(f"  {label} -> EdgeSplitError[{exc.reason}]，已拒绝")
        except MeshInputError as exc:
            print(f"  {label} -> MeshInputError 于 {exc.location}，已拒绝")
    # 原子性复核：被拒后原网格不变。
    assert len(kuhn.tets) == 6 and audit(kuhn).total_signed_volume == 1
    assert kuhn.protected_faces == ((0, 1, 2),)
    print("  拒绝后原网格单元数/体积/保护面均保持不变（无半更新）。")

    # ============ C. 原有缺陷边界：非流形报告 / 零体积拒绝 / 重复体拒绝 ============
    print()
    print("=" * 68)
    print("C. 缺陷边界：非流形报告 / 零体积拒绝 / 重复体拒绝")
    print("=" * 68)

    # C1. 三个四面体共享同一个三点面 (0,1,2) -> 非流形。
    nm_verts = (
        (0, 0, 0), (1, 1, 0), (1, 0, 1),
        (1, 0, 0), (0, 1, 1), (0, 1, 0),
    )
    nm_tets = [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5)]
    nm_mesh = TetMesh.build(nm_verts, nm_tets)
    nm_report = audit(nm_mesh)
    print("[C1] 三单元共享面：")
    for f in nm_report.non_manifold_faces:
        print(f"  非流形面 key={f.key}，被单元 {list(f.tet_indices)} 引用"
              f"（{len(f.tet_indices)} > 2）")
    assert not nm_report.is_manifold

    # C2. 加入一个共面（零体积）单元后，修复必须被拒绝。
    deg_verts = (*VERTICES, (1, 1, 0))  # 下标 8 与 0,1,2 共面于 z=0
    deg_mesh = TetMesh.build(deg_verts, [*TETS, (0, 1, 2, 8)])
    deg_report = audit(deg_mesh)
    print("[C2] 共面退化：")
    print(f"  audit 精确定位零体积单元: {list(deg_report.zero_tets)}（不静默纠正）")
    try:
        repair_orientation(deg_mesh)
    except ZeroVolumeError as exc:
        print(f"  repair_orientation 拒绝修复 -> {type(exc).__name__}:")
        print(f"    {exc}")

    # C3. 重复单元单独报错（排列不同也算重复）。
    print("[C3] 重复单元：")
    try:
        perm = (TETS[0][2], TETS[0][0], TETS[0][3], TETS[0][1])
        TetMesh.build(VERTICES, [*TETS, perm])
    except DuplicateTetrahedraError as exc:
        print(f"  TetMesh.build 拒绝 -> {type(exc).__name__}:")
        key, positions = exc.duplicates[0]
        print(f"    单元 {positions} 顶点集合均为 {key}")

    # C4. 非法输入可定位（顺带演示错误语义，非核心几何边界）。
    print("[C4] 非法输入：")
    for label, build in (
        ("float 坐标", lambda: TetMesh.build([(0, 0, 0.5)], [])),
        ("bool 下标", lambda: TetMesh.build(VERTICES, [(0, 1, 2, True)])),
    ):
        try:
            build()
        except MeshInputError as exc:
            print(f"  {label} -> {type(exc).__name__} 于 {exc.location}: 已拒绝")

    elapsed = time.perf_counter() - t0
    print()
    print(f"演示完成，实际耗时 {elapsed:.3f} 秒。")


if __name__ == "__main__":
    main()
