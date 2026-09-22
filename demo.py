"""固定输入演示：四面体网格定向与边界审计。

运行：python demo.py（约 1 秒内完成，纯标准库、无网络）。

内容全部为真实计算：
  1. 一致定向的单位立方体 6 四面体剖分：解析体积 + 独立面通量核对；
  2. 翻转一个单元：审计定位、按规则修复（仅交换最后两顶点）、核对恢复；
  3. 非流形网格：修复单元符号后仍报告三面关联缺陷；
  4. 拒绝边界（实际触发的异常）：零体积单元拒绝修复，NaN 坐标拒绝构造。
"""

from fractions import Fraction

from tetmesh_audit import (
    DegenerateTetrahedronError,
    InvalidMeshError,
    TetMesh,
    audit,
    repair_orientation,
    samples,
)


def independent_flux(vertices, oriented_triangles):
    """独立通量公式 V = 1/6 Σ det(p,q,r)（逐项展开，不经过被测库）。"""
    def det3(a, b, c):
        return (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        )

    total = 0
    for tri in oriented_triangles:
        total += det3(vertices[tri[0]], vertices[tri[1]], vertices[tri[2]])
    return Fraction(total, 6)


def section(title):
    print()
    print("=" * 68)
    print(title)
    print("=" * 68)


def main():
    # ---- 1. 正常结果：一致立方体 ----
    section("1. 正常网格：单位立方体 6 四面体扇形剖分")
    mesh = samples.cube_tetmesh()
    res = audit(mesh)
    print(res.summary())
    print(f"  逐单元 6 倍有向体积: {list(res.six_volumes)}")
    print(f"  总有向体积（精确）  : {res.total_signed_volume}  解析值 = 1")
    print(f"  边界定向三角面数    : {len(res.boundary_triangles)}")
    print(f"  内面数（恰两单元）  : {len(res.interior_faces)}")
    print(f"  连通分量数          : {len(res.components)}")
    tris = [bt.vertices for bt in res.boundary_triangles]
    flux = independent_flux(mesh.vertices, tris)
    print(f"  独立面通量核对 V    : {flux}  (与审计总体积 {'相等' if flux == res.total_signed_volume else '不一致'})")
    assert flux == res.total_signed_volume == 1
    assert res.is_consistent and res.defects == ()

    # ---- 2. 修复一个翻转单元 ----
    section("2. 修复演示：第 2 号单元被奇置换翻转")
    bad = samples.cube_with_flipped_tet(2)
    before = audit(bad)
    print(f"  修复前: 负向单元 {before.negative_tets}，"
          f"总体积 {before.total_signed_volume}")
    bad_faces = [d.face for d in before.defects
                 if d.kind == "inconsistent_interior_face"]
    print(f"  定向冲突内面位置: {bad_faces}")
    repaired = repair_orientation(bad)
    old = bad.tetrahedra[2]
    new = repaired.repaired_mesh.tetrahedra[2]
    print(f"  修复动作: tet[2] {old} -> {new}（仅交换最后两顶点，顶点不动）")
    after = repaired.audit_after
    print(f"  修复后: {after.summary()}")
    print(f"  修复后顶点坐标与修复前逐点相等: "
          f"{repaired.repaired_mesh.vertices == bad.vertices}")
    flux2 = independent_flux(repaired.repaired_mesh.vertices,
                             [bt.vertices for bt in after.boundary_triangles])
    print(f"  修复后独立面通量核对: {flux2}")
    assert after.is_consistent and flux2 == 1

    # ---- 3. 非流形：修复符号不能消除拓扑缺陷 ----
    section("3. 非流形网格：面 (0,1,2) 关联 3 个单元")
    nm_mesh = samples.nonmanifold_mesh()
    nm_rep = repair_orientation(nm_mesh)
    nma = nm_rep.audit_after
    print(f"  符号修复翻转的单元: {nm_rep.flipped_tets}")
    for d in nma.defects:
        print(f"  [缺陷] {d.kind}: {d.detail}")
    print(f"  修复后仍一致? {nma.is_consistent}（非流形不被符号修复消除）")
    assert any(d.kind == "non_manifold_face" for d in nma.defects)

    # ---- 4a. 实际触发的拒绝异常：零体积单元 ----
    section("4. 拒绝/异常边界（实际触发并捕获）")
    deg = samples.coplanar_degenerate_mesh()
    deg_res = audit(deg)
    print("4a. 共面退化单元：")
    print(f"  审计定位零体积单元: {deg_res.zero_volume_tets}")
    try:
        repair_orientation(deg)
    except DegenerateTetrahedronError as exc:
        print(f"  repair_orientation 拒绝修复 -> {type(exc).__name__}: {exc}")
        print("  （原子失败：未产生部分修复，原网格不变）")

    # ---- 4b. 实际触发的输入异常：NaN 坐标 ----
    print("4b. 非法数值输入：")
    bad_verts = [list(v) for v in samples.CUBE_VERTICES]
    bad_verts[5][2] = float("nan")
    try:
        TetMesh(bad_verts, samples.CUBE_TETS)
    except InvalidMeshError as exc:
        print(f"  TetMesh 拒绝构造 -> {type(exc).__name__}: {exc}")

    section("演示完成：以上数值均为当场精确计算，无预录结果。")


if __name__ == "__main__":
    main()
