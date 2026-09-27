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

运行：python demo.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import (
    DuplicateTetrahedraError,
    MergePointError,
    MeshInputError,
    TetMesh,
    ZeroVolumeError,
    audit,
    merge_meshes,
    repair_orientation,
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

    # ============ C. 两份网格显式点对合并 ============
    print()
    print("=" * 68)
    print("C. 两立方体拼接：只合并显式指定的 4 个接合面点（无自动吸附）")
    print("=" * 68)

    # B 是 A 关于 x=1 的镜像剖分（接合面两侧同对角，共形）。
    b_verts = (
        (2, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 0),
        (2, 0, 1), (1, 0, 1), (1, 1, 1), (2, 1, 1),
    )
    b_tets = (
        (0, 1, 5, 2),
        (0, 2, 7, 3),
        (0, 4, 7, 5),
        (2, 5, 7, 6),
        (0, 2, 5, 7),
    )
    mesh_a = TetMesh.build(VERTICES, TETS)
    mesh_b = TetMesh.build(b_verts, b_tets)
    # 只列出接合面 x=1 上坐标完全相同的 4 对顶点（两侧编号恰都是 1,2,5,6）。
    pairs = ((1, 1), (2, 2), (5, 5), (6, 6))
    merged = merge_meshes(mesh_a, mesh_b, pairs)
    mr = merged.report
    print("  A 旧编号 -> 新编号:", list(merged.vertex_map_a))
    print("  B 旧编号 -> 新编号:", list(merged.vertex_map_b))
    print(f"  合并顶点数: 8+8-4 = {len(merged.mesh.vertices)}"
          f"，单元数: {len(merged.mesh.tets)}")
    print(f"  总精确有向体积: {mr.total_signed_volume}（1+1）")
    print(f"  边界面数: {len(mr.boundary_faces)}（两个立方体外表面 24 - 接合面 4）")
    print(f"  非流形面: {len(mr.non_manifold_faces)}，连通分量: {len(mr.components)}")
    assert mr.total_signed_volume == Fraction(2)
    assert len(mr.boundary_faces) == 20
    assert len(mr.components) == 1
    assert not mr.has_non_manifold

    # 不列出任何点对：即使坐标重合也各自独立，绝不自动吸附。
    separate = merge_meshes(mesh_a, mesh_b, ())
    print(f"  不给点对时: 顶点数 {len(separate.mesh.vertices)}（16，不合并），"
          f"连通分量 {len(separate.report.components)}")
    assert len(separate.mesh.vertices) == 16

    # 指定一对位置不同的点：拒绝（不按距离近吸附），原网格保持不变。
    print("  错误点对（A[1] 与 B[0] 坐标不同）：")
    try:
        merge_meshes(mesh_a, mesh_b, ((1, 0),))
    except MergePointError as exc:
        print(f"    -> {type(exc).__name__} 于 {exc.location}：已拒绝，不留下半成品")
    assert len(mesh_b.tets) == 5

    elapsed = time.perf_counter() - t0
    print()
    print(f"演示完成，实际耗时 {elapsed:.3f} 秒。")


if __name__ == "__main__":
    main()
