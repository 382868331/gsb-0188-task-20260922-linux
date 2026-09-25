#!/usr/bin/env python3
"""固定输入演示：四面体定向与边界审计。

真实计算，无 sleep、无预录结果、无网络。两部分：

A. 正常结果：单位立方体 5 四面体剖分中故意翻转两个单元 ->
   audit 定位负向单元 -> repair_orientation 修复（只交换最后两顶点）->
   再审计：精确总体积 1、12 个外法向边界面、单一连通分量。
   总体积另用*独立*的面通量公式 V = (1/6)·Σ det(p,q,r) 核对
   （该公式不经过任何单元体积，也不调用本库核心）。

B. 实际触发的边界：三单元共享一个面 -> 报告非流形位置；
   加入共面退化单元 -> repair_orientation 抛 ZeroVolumeError 拒绝修复；
   重复单元 -> TetMesh.build 抛 DuplicateTetrahedraError。

C. 平面截面（新增）：有理平面 × 四面体网格 -> 有理多边形片。
   共享棱交点身份核对、解析棱锥面积核对、共面面去重、零面积接触。

运行：python demo.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import (
    DuplicateTetrahedraError,
    MeshInputError,
    TetMesh,
    ZeroVolumeError,
    audit,
    repair_orientation,
    section,
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

    # ============ C. 平面截面：有理平面 × 四面体网格 ============
    print()
    print("=" * 68)
    print("C. 平面截面：有理多边形片 / 点身份 / 解析面积核对")
    print("=" * 68)

    # 两个共享面 (0,1,2) 的四面体（顶点 3 在 z>0，顶点 4 在 z<0）。
    two_verts = (
        (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1),
    )
    two_mesh = TetMesh.build(two_verts, [(0, 1, 2, 3), (0, 2, 1, 4)])

    # C1. 平面 x = 1/3 截两个单元：共享棱交点身份必须一致（原顶点对）。
    print("[C1] 两个共享面的四面体 × 平面 x = 1/3：")
    rep1 = section(two_mesh, (1, 0, 0), Fraction(1, 3))
    for p in rep1.patches:
        print(f"  片 kind={p.kind} 来源单元={list(p.source_tets)} "
              f"多边形={list(p.polygon)} 精确面积平方={p.area2}")
    shared = set(rep1.patches[0].polygon) & set(rep1.patches[1].polygon)
    print(f"  两片共享的全局交点 ID: {sorted(shared)}，身份: "
          f"{sorted(rep1.points[i].origin for i in shared)}")
    assert all(r == 0 for r in rep1.plane_residuals())
    print("  全部截面点精确满足平面方程（残差均为 0）。")

    # C2. 解析棱锥核对：z = 1/2 截方底棱锥，截面是边长 1 的正方形。
    print("[C2] 方底棱锥 × 平面 z = 1/2（解析截面面积 = 1）：")
    pyr_mesh = TetMesh.build(
        [(-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0), (0, 0, 1)],
        [(0, 1, 2, 4), (0, 2, 3, 4)],
    )
    rep2 = section(pyr_mesh, (0, 0, 1), Fraction(1, 2))
    print(f"  面片数: {len(rep2.patches)}，"
          f"各片精确面积平方: {[str(p.area2) for p in rep2.patches]}")
    print(f"  截面总面积（浮点便利出口）: {rep2.total_area()}")
    assert all(p.area2 == Fraction(1, 4) for p in rep2.patches)
    assert rep2.total_area() == 1.0
    print("  与解析值 1 一致（每片精确面积平方 = 1/4）。")

    # C3. 平面与共享面重合：完全共面的原面只保留一份。
    print("[C3] 平面 z = 0 与共享面重合：")
    rep3 = section(two_mesh, (0, 0, 1), 0)
    cop = rep3.patches[0]
    print(f"  面片数: {len(rep3.patches)}（去重后仅一份），"
          f"kind={cop.kind}，来源单元={list(cop.source_tets)}")
    assert len(rep3.patches) == 1 and cop.source_tets == (0, 1)

    # C4. 零面积接触：平面 x = 1 只在顶点 1 处碰到两个单元。
    print("[C4] 平面 x = 1 仅触顶点（零面积结果单独列出）：")
    rep4 = section(two_mesh, (1, 0, 0), 1)
    print(f"  面片数: {len(rep4.patches)}，"
          f"点接触: {[(c.vertex, list(c.source_tets)) for c in rep4.point_contacts]}")
    assert rep4.patches == () and len(rep4.point_contacts) == 1

    # C5. 非法平面参数实际触发拒绝。
    print("[C5] 非法平面参数：")
    try:
        section(two_mesh, (0, 0, 0), 0)
    except MeshInputError as exc:
        print(f"  零法向 -> {type(exc).__name__} 于 {exc.location}: 已拒绝")

    elapsed = time.perf_counter() - t0
    print()
    print(f"演示完成，实际耗时 {elapsed:.3f} 秒。")
if __name__ == "__main__":
    main()
