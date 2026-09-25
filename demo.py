#!/usr/bin/env python3
"""固定输入演示：四面体定向审计（A/B）与有理平面剖切（C）。

真实计算，无 sleep、无预录结果、无网络。三部分：

A. 正常结果：单位立方体 5 四面体剖分中故意翻转两个单元 ->
   audit 定位负向单元 -> repair_orientation 修复（只交换最后两顶点）->
   再审计：精确总体积 1、12 个外法向边界面、单一连通分量。
   总体积另用*独立*的面通量公式 V = (1/6)·Σ det(p,q,r) 核对
   （该公式不经过任何单元体积，也不调用本库核心）。

B. 实际触发的边界：三单元共享一个面 -> 报告非流形位置；
   加入共面退化单元 -> repair_orientation 抛 ZeroVolumeError 拒绝修复；
   重复单元 -> TetMesh.build 抛 DuplicateTetrahedraError。

C. 新增：平面 section 四个场景——两个共享面四面体、穿顶点、
   平面沿共享内面、穿四边形截面。每个交点精确验证平面方程，
   片内点用重心坐标验证位于来源单元，面积用解析棱锥
   V = 底面积·高/3 独立核对（不调用被测核心）。

运行：./.venv/bin/python demo.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import (
    DuplicateTetrahedraError,
    MeshInputError,
    Plane,
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

    # ============ C. 有理平面剖切 section（新增能力） ============
    print()
    print("=" * 68)
    print("C. 平面剖切：有理 n·x = d，交点身份来自原顶点对")
    print("=" * 68)

    def on_plane(coord, n, d) -> bool:
        return sum(n[k] * coord[k] for k in range(3)) == d

    def barycentric(coord, tet_pts):
        """4x4 体积比重心坐标（独立于被测核心）。"""

        def vol6(a, b, c, d4):
            m = [[Fraction(1), *a], [Fraction(1), *b],
                 [Fraction(1), *c], [Fraction(1), *d4]]
            prev, sign = Fraction(1), 1
            for k in range(3):
                if m[k][k] == 0:
                    for i in range(k + 1, 4):
                        if m[i][k]:
                            m[k], m[i] = m[i], m[k]
                            sign = -sign
                            break
                piv = m[k][k]
                bm = [row[:] for row in m]
                for i in range(k + 1, 4):
                    for j in range(k + 1, 4):
                        bm[i][j] = (m[i][j] * piv - m[i][k] * m[k][j]) / prev
                m, prev = bm, piv
            return sign * m[3][3]

        p = tuple(Fraction(x) for x in coord)
        q = [tuple(Fraction(x) for x in t) for t in tet_pts]
        tot = vol6(*q)
        return (
            vol6(p, q[1], q[2], q[3]) / tot,
            vol6(q[0], p, q[2], q[3]) / tot,
            vol6(q[0], q[1], p, q[3]) / tot,
            vol6(q[0], q[1], q[2], p) / tot,
        )

    def pyramid_alpha(tri_coords, apex, n, d) -> Fraction:
        """解析棱锥：面积/|n| = |det6|/(2|n·apex−d|)，独立核对。"""

        def vol6(a, b, c, d4):
            m = [[Fraction(1), *a], [Fraction(1), *b],
                 [Fraction(1), *c], [Fraction(1), *d4]]
            prev, sign = Fraction(1), 1
            for k in range(3):
                if m[k][k] == 0:
                    for i in range(k + 1, 4):
                        if m[i][k]:
                            m[k], m[i] = m[i], m[k]
                            sign = -sign
                            break
                piv = m[k][k]
                bm = [row[:] for row in m]
                for i in range(k + 1, 4):
                    for j in range(k + 1, 4):
                        bm[i][j] = (m[i][j] * piv - m[i][k] * m[k][j]) / prev
                m, prev = bm, piv
            return sign * m[3][3]

        det6 = abs(vol6(tuple(apex), *[tuple(c) for c in tri_coords]))
        height_num = abs(
            sum(n[k] * Fraction(apex[k]) for k in range(3)) - d
        )
        return det6 / (2 * height_num)

    def verify_report(name, mesh, plane, expected_tets):
        rep = section(mesh, plane)
        cmap = {p.id: p.coord for p in rep.points}
        print(f"[{name}] 平面 {tuple(plane.normal)}·x = {plane.offset}")
        for p in rep.points:
            assert on_plane(p.coord, plane.normal, plane.offset), p
        for poly in rep.polygons:
            shape = "四边形" if poly.is_quad else "三角形"
            print(
                f"  tets[{poly.tet_index}] {shape}片 顶点={list(poly.vertices)}"
                f" α={poly.area_coefficient}（面积=α·|n|）"
            )
            assert poly.tet_index in expected_tets
            for tri in poly.triangles:
                qs = [cmap[v] for v in tri]
                center = tuple(sum(p[k] for p in qs) / 3 for k in range(3))
                lams = barycentric(
                    center, [mesh.vertices[i] for i in mesh.tets[poly.tet_index]]
                )
                assert all(l >= 0 for l in lams) and sum(lams) == 1
                sign3 = sum(
                    plane.normal[k]
                    * (
                        (qs[1][(k + 1) % 3] - qs[0][(k + 1) % 3])
                        * (qs[2][(k + 2) % 3] - qs[0][(k + 2) % 3])
                        - (qs[1][(k + 2) % 3] - qs[0][(k + 2) % 3])
                        * (qs[2][(k + 1) % 3] - qs[0][(k + 1) % 3])
                    )
                    for k in range(3)
                )
                assert sign3 > 0  # 朝向统一为给定法向
        for f in rep.coplanar_faces:
            print(
                f"  共面原面 key={f.key} 来源单元={list(f.tet_indices)}"
                f"（全局只一份）α={f.area_coefficient}"
            )
        for c in rep.contact_points:
            print(f"  零面积点接触 {c.id} 来源单元={list(c.tet_indices)}")
        for c in rep.contact_edges:
            print(f"  零面积棱接触 {tuple(c.endpoints)} 来源单元={list(c.tet_indices)}")
        print(f"  全部 {len(rep.points)} 个交点精确满足平面方程；"
              f"片内重心在来源单元内；总 α={rep.total_area_coefficient}")
        return rep

    # C1. 两个共享面四面体，平面同时穿过它们（共享棱交点同一身份）。
    bi_verts = (
        (0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4), (0, 0, -2),
    )
    bi_tets = [(0, 1, 2, 3), (1, 0, 2, 4)]
    bi_mesh = TetMesh.build(bi_verts, bi_tets)
    plane1 = Plane.build((0, 1, 0), 1)
    r1 = verify_report("C1 共享面双四面体 y=1", bi_mesh, plane1, {0, 1})
    shared = ("edge", (1, 2))
    assert len([p for p in r1.points if p.id == shared]) == 1
    assert all(shared in poly.vertices for poly in r1.polygons)
    cmap1 = {p.id: p.coord for p in r1.points}
    # 解析棱锥独立核对两片面积（顶点 3 在 +y 侧、顶点 4 在 −y 侧）。
    assert pyramid_alpha(
        [cmap1[v] for v in r1.polygons[0].triangles[0]],
        bi_verts[3], plane1.normal, plane1.offset,
    ) == Fraction(9, 2) == r1.polygons[0].area_coefficient
    assert pyramid_alpha(
        [cmap1[v] for v in r1.polygons[1].triangles[0]],
        bi_verts[4], plane1.normal, plane1.offset,
    ) == Fraction(9, 4) == r1.polygons[1].area_coefficient
    print("  解析棱锥核对：两片 α = 9/2、9/4，与核心一致。")

    # C2. 平面穿过原顶点（z=0 含立方体底面四点）：点/棱接触单列。
    cube = TetMesh.build(VERTICES, TETS)
    r2 = verify_report("C2 穿顶点 z=0", cube, Plane.build((0, 0, 1), 0), set())
    assert {c.id for c in r2.contact_points} == {("vertex", 0), ("vertex", 2)}
    assert len(r2.contact_edges) == 1 and not r2.polygons
    assert sorted(f.key for f in r2.coplanar_faces) == [(0, 1, 2), (0, 2, 3)]

    # C3. 平面沿共享内面 (0,2,5)：共面原面只保留一份、记录 t0 与 t4。
    plane3 = Plane.build((1, -1, -1), 0)
    r3 = verify_report("C3 沿共享内面", cube, plane3, set())
    face = next(f for f in r3.coplanar_faces if f.key == (0, 2, 5))
    assert face.tet_indices == (0, 4)
    assert len([f for f in r3.coplanar_faces if f.key == (0, 2, 5)]) == 1
    print("  共面内面 (0,2,5) 仅一份，来源单元 (0,4)；其余单元只棱接触。")

    # C4. 穿四边形截面 z=1/2：4 个三角片 + 中心单元的四边形片。
    plane4 = Plane.build((0, 0, 1), Fraction(1, 2))
    r4 = verify_report("C4 四边形截面 z=1/2", cube, plane4, {0, 1, 2, 3, 4})
    quads = [p for p in r4.polygons if p.is_quad]
    assert len(quads) == 1 and quads[0].tet_index == 4
    assert len(quads[0].triangles) == 2
    # 总截面是单位正方形，α 精确为 1。
    assert r4.total_area_coefficient == Fraction(1)
    # 多单元共享棱 (2,5) 的交点全局唯一。
    assert len([p for p in r4.points if p.id == ("edge", (2, 5))]) == 1
    print("  四边形片属 tets[4]，按全局交点 ID 定对角线；总截面 α = 1（单位正方形）。")

    # C5. 负向输入单元：先按原规则（只交换最后两顶点）处理，原网格不变。
    neg_tets = [list(t) for t in TETS]
    neg_tets[1] = [neg_tets[1][0], neg_tets[1][1], neg_tets[1][3], neg_tets[1][2]]
    neg_mesh = TetMesh.build(VERTICES, neg_tets)
    r5a = section(neg_mesh, plane4)
    r5b = section(cube, plane4)
    assert [p.vertices for p in r5a.polygons] == [p.vertices for p in r5b.polygons]
    assert neg_mesh.tets[1] == tuple(neg_tets[1])  # 原网格未被修改/分裂
    print("[C5] 负向单元先修正再剖切，截面与正向网格一致；原网格顶点/棱未改动。")

    elapsed = time.perf_counter() - t0
    print()
    print(f"演示完成，实际耗时 {elapsed:.3f} 秒。")


if __name__ == "__main__":
    main()
