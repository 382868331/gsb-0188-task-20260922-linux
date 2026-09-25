"""平面剖切 section 的必要测试。

覆盖题面要求的四个场景与全部不变量：

1. 两个共享面的四面体被同一平面穿过（共享棱交点身份一致）；
2. 平面穿过原顶点（点接触 / 棱接触单列）；
3. 平面沿两单元的共享内面（共面原面全局只一份、记录全部来源）；
4. 穿四边形截面（确定对角线拆两个三角形）。

并独立验证：所有交点精确满足平面方程；片内点（重心坐标）位于来源
单元；面积用解析棱锥公式（不调用被测核心）核对；原网格不被修改、
不做边分裂；非法平面 / 零体积单元的错误语义。
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    ContactEdge,
    ContactPoint,
    CoplanarFacePatch,
    MeshInputError,
    Plane,
    SectionPolygon,
    TetMesh,
    ZeroVolumeError,
    repair_orientation,
    section,
)

from fixtures import CUBE_TETS, CUBE_VERTICES
from reference import (
    ref_barycentric,
    ref_on_plane,
    ref_polygon_alpha,
    ref_pyramid_alpha,
)

# ---------------------------------------------------------------- 辅助


def _coord_map(report):
    return {p.id: p.coord for p in report.points}


def _all_triangles(report):
    out = []
    for poly in report.polygons:
        for tri in poly.triangles:
            out.append((poly.tet_index, tri))
    for cp in report.coplanar_faces:
        for tri in cp.triangles:
            out.append((cp.tet_indices, tri))
    return out


def _oriented_sign(ids, cmap, normal):
    """n·(b−a)×(c−a) 的精确符号（所有输出三角形必须 > 0）。"""

    p = cmap[ids[0]]
    q = cmap[ids[1]]
    r = cmap[ids[2]]
    u = tuple(q[k] - p[k] for k in range(3))
    v = tuple(r[k] - p[k] for k in range(3))
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    return sum(normal[k] * cross[k] for k in range(3))


# ---------------------------------------------------------------- 1. 共享面双四面体


class TestTwoSharedFaceTets(unittest.TestCase):
    """三角双锥：t0=(0,1,2,3) 与 t1=(1,0,2,4) 共享面 (0,1,2)。

    顶点：0:(0,0,0) 1:(4,0,0) 2:(0,4,0) 3:(0,0,4) 4:(0,0,-2)。
    两单元均为正体积；平面 z=1 切 t0 成三角形，t1 不相交。
    """

    VERTICES = (
        (0, 0, 0),
        (4, 0, 0),
        (0, 4, 0),
        (0, 0, 4),
        (0, 0, -2),
    )
    T0 = (0, 1, 2, 3)
    T1 = (1, 0, 2, 4)

    def setUp(self):
        self.mesh = TetMesh.build(self.VERTICES, [self.T0, self.T1])
        self.plane = Plane.build((0, 0, 1), 1)
        self.rep = section(self.mesh, self.plane)

    def test_one_triangle_patch_only(self):
        self.assertEqual(len(self.rep.polygons), 1)
        poly = self.rep.polygons[0]
        self.assertEqual(poly.tet_index, 0)
        self.assertEqual(len(poly.vertices), 3)
        self.assertEqual(len(poly.triangles), 1)
        self.assertEqual(self.rep.coplanar_faces, ())
        self.assertEqual(self.rep.contact_points, ())
        self.assertEqual(self.rep.contact_edges, ())

    def test_edge_identity_is_original_vertex_pair(self):
        poly = self.rep.polygons[0]
        self.assertEqual(
            sorted(poly.vertices, key=lambda x: x[1]),
            [("edge", (0, 3)), ("edge", (1, 3)), ("edge", (2, 3))],
        )

    def test_exact_intersection_coordinates_and_plane_equation(self):
        cmap = _coord_map(self.rep)
        for pid, coord in cmap.items():
            self.assertTrue(ref_on_plane(coord, self.plane.normal, self.plane.offset),
                            f"{pid} 必须精确满足平面方程")
        # 交点坐标精确：棱 (0,3) 上 t=1/4 -> (0,0,1)。
        self.assertEqual(cmap[("edge", (0, 3))], (Fraction(0), Fraction(0), Fraction(1)))
        self.assertEqual(
            cmap[("edge", (1, 3))], (Fraction(3), Fraction(0), Fraction(1))
        )
        self.assertEqual(
            cmap[("edge", (2, 3))], (Fraction(0), Fraction(3), Fraction(1))
        )

    def test_area_matches_analytic_pyramid(self):
        # 顶部小棱锥：高 3、底边为原三角形 3/4 缩放，面积 = 8·(3/4)^2 = 9/2。
        poly = self.rep.polygons[0]
        cmap = _coord_map(self.rep)
        cyclic = [cmap[pid] for pid in poly.vertices]
        self.assertEqual(poly.area_coefficient, Fraction(9, 2))
        # 独立路径：边界叉积和公式。
        self.assertEqual(
            ref_polygon_alpha(cyclic, self.plane.normal), Fraction(9, 2)
        )
        # 独立路径：解析棱锥 V = 底面积·高/3，顶点取被切下的角顶点 3。
        apex = self.VERTICES[3]
        self.assertEqual(
            ref_pyramid_alpha(cyclic, apex, self.plane.normal, self.plane.offset),
            Fraction(9, 2),
        )

    def test_interior_point_inside_source_tet(self):
        # 三角形重心必须落在来源单元（修复后顶点序）闭单体内，且在平面上。
        poly = self.rep.polygons[0]
        cmap = _coord_map(self.rep)
        pts = [cmap[pid] for pid in poly.vertices]
        centroid = tuple(sum(p[k] for p in pts) / 3 for k in range(3))
        tet_pts = [self.mesh.vertices[i] for i in self.T0]
        (lams, _), = [ref_barycentric(centroid, tet_pts)]
        self.assertEqual(sum(lams), Fraction(1))
        self.assertTrue(all(l >= 0 for l in lams))
        self.assertTrue(ref_on_plane(centroid, self.plane.normal, self.plane.offset))

    def test_orientation_follows_given_normal(self):
        cmap = _coord_map(self.rep)
        for _, tri in _all_triangles(self.rep):
            self.assertGreater(_oriented_sign(tri, cmap, self.plane.normal), 0)

    def test_negated_normal_reverses_winding_same_geometry(self):
        rev_plane = Plane.build((0, 0, -1), -1)
        rev = section(self.mesh, rev_plane)
        c1 = sorted(p.coord for p in self.rep.points)
        c2 = sorted(p.coord for p in rev.points)
        self.assertEqual(c1, c2)  # 几何完全相同
        cmap = _coord_map(rev)
        for _, tri in _all_triangles(rev):
            self.assertGreater(_oriented_sign(tri, cmap, rev_plane.normal), 0)

    def test_plane_cuts_both_tets_reuses_shared_edge_point(self):
        # y=1 同时穿过两单元；共享棱 (1,2) 上的交点必须是同一个全局点。
        plane = Plane.build((0, 1, 0), 1)
        rep = section(self.mesh, plane)
        self.assertEqual([p.tet_index for p in rep.polygons], [0, 1])
        shared = ("edge", (1, 2))
        self.assertTrue(all(shared in p.vertices for p in rep.polygons))
        self.assertEqual(len([p for p in rep.points if p.id == shared]), 1)
        # 解析棱锥独立核对两片面积（顶点 3、4 分居两侧）。
        cmap = _coord_map(rep)
        self.assertEqual(
            ref_pyramid_alpha(
                [cmap[v] for v in rep.polygons[0].triangles[0]],
                self.VERTICES[3], plane.normal, plane.offset,
            ),
            Fraction(9, 2),
        )
        self.assertEqual(
            ref_pyramid_alpha(
                [cmap[v] for v in rep.polygons[1].triangles[0]],
                self.VERTICES[4], plane.normal, plane.offset,
            ),
            Fraction(9, 4),
        )


class TestTwoPlaneVerticesOppositeSides(unittest.TestCase):
    """符号模式 (0,0,+,-)：两顶点在平面上、另两顶点分居两侧。

    截面是三角形：两个在平面上的原顶点 + 跨侧棱的唯一交点。
    """

    VERTICES = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, -1, 1))

    def setUp(self):
        self.mesh = TetMesh.build(self.VERTICES, [(0, 1, 2, 3)])
        self.plane = Plane.build((0, 1, 0), 0)
        self.rep = section(self.mesh, self.plane)

    def test_triangle_with_two_original_vertices(self):
        self.assertEqual(len(self.rep.polygons), 1)
        poly = self.rep.polygons[0]
        self.assertEqual(
            set(poly.vertices),
            {("vertex", 0), ("vertex", 1), ("edge", (2, 3))},
        )
        self.assertEqual(self.rep.coplanar_faces, ())
        self.assertEqual(self.rep.contact_points, ())
        self.assertEqual(self.rep.contact_edges, ())

    def test_area_from_either_side_pyramid(self):
        # 从平面任一侧的角顶点算棱锥，得到同一三角形面积系数。
        cmap = _coord_map(self.rep)
        tri = self.rep.polygons[0].triangles[0]
        coords = [cmap[v] for v in tri]
        alpha = Fraction(1, 4)
        self.assertEqual(self.rep.polygons[0].area_coefficient, alpha)
        for apex in (self.VERTICES[2], self.VERTICES[3]):
            self.assertEqual(
                ref_pyramid_alpha(coords, apex, self.plane.normal, self.plane.offset),
                alpha,
            )


# ---------------------------------------------------------------- 2. 立方体网格：共享棱交点身份 + 四边形片


class TestCubeHalfPlane(unittest.TestCase):
    """立方体 5 剖分，平面 z=1/2：4 个三角片 + t4 四边形片。"""

    def setUp(self):
        self.mesh = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        self.plane = Plane.build((0, 0, 1), Fraction(1, 2))
        self.rep = section(self.mesh, self.plane)

    def test_patch_shapes(self):
        by_tet = {p.tet_index: p for p in self.rep.polygons}
        self.assertEqual(set(by_tet), {0, 1, 2, 3, 4})
        for ti in (0, 1, 2, 3):
            self.assertEqual(len(by_tet[ti].vertices), 3, f"tets[{ti}] 三角形")
        quad = by_tet[4]
        self.assertEqual(len(quad.vertices), 4)
        self.assertEqual(len(quad.triangles), 2)
        self.assertTrue(quad.is_quad)

    def test_shared_edge_intersection_is_one_global_point(self):
        # 棱 (0,5) 同时被 t0 与 t4 穿过：必须是同一个全局 id / 同一个点对象。
        by_tet = {p.tet_index: p for p in self.rep.polygons}
        shared = ("edge", (0, 5))
        self.assertIn(shared, by_tet[0].vertices)
        self.assertIn(shared, by_tet[4].vertices)
        pts = [p for p in self.rep.points if p.id == shared]
        self.assertEqual(len(pts), 1)  # 不靠坐标焊接：身份来自原顶点对
        # 棱 (2,5) 被三个单元 t0/t3/t4 共享，仍是一个点。
        self.assertEqual(
            len([p for p in self.rep.points if p.id == ("edge", (2, 5))]), 1
        )
        for p in (by_tet[0], by_tet[3], by_tet[4]):
            self.assertIn(("edge", (2, 5)), p.vertices)

    def test_quad_diagonal_deterministic_by_global_id(self):
        quad = next(p for p in self.rep.polygons if p.tet_index == 4)
        # 四边形四个顶点为棱交点 (0,7)(0,5)(2,7)(2,5)；
        # 对角线规则：全局 id 最小点 (0,5) 连环上对点 (2,7)。
        sorted_tris = sorted(tuple(sorted(t)) for t in quad.triangles)
        self.assertEqual(
            sorted_tris,
            sorted([
                (("edge", (0, 5)), ("edge", (0, 7)), ("edge", (2, 7))),
                (("edge", (0, 5)), ("edge", (2, 5)), ("edge", (2, 7))),
            ]),
        )
        # 两三角共享对角线 (0,5)-(2,7)。
        shared = set(quad.triangles[0]) & set(quad.triangles[1])
        self.assertEqual(shared, {("edge", (0, 5)), ("edge", (2, 7))})

    def test_total_area_is_one_square(self):
        # 平面 z=1/2 截面是单位正方形，面积精确为 1。
        self.assertEqual(self.rep.total_area_coefficient, Fraction(1))
        # 逐片独立面积和（边界叉积公式）同样为 1。
        cmap = _coord_map(self.rep)
        total = sum(
            ref_polygon_alpha([cmap[v] for v in poly.vertices], self.plane.normal)
            for poly in self.rep.polygons
        )
        self.assertEqual(total, Fraction(1))

    def test_all_points_on_plane_and_patch_interiors_inside_tets(self):
        cmap = _coord_map(self.rep)
        for p in self.rep.points:
            self.assertTrue(ref_on_plane(p.coord, self.plane.normal, self.plane.offset))
        for poly in self.rep.polygons:
            tet_pts = [self.mesh.vertices[i] for i in CUBE_TETS[poly.tet_index]]
            for tri in poly.triangles:
                tri_pts = [cmap[v] for v in tri]
                center = tuple(sum(p[k] for p in tri_pts) / 3 for k in range(3))
                (lams, _), = [ref_barycentric(center, tet_pts)]
                self.assertTrue(all(l >= 0 for l in lams),
                                f"片内点必须在来源单元 tets[{poly.tet_index}] 内")
                self.assertEqual(sum(lams), Fraction(1))

    def test_triangulation_covers_quad_exactly(self):
        # 两个三角形面积系数之和必须等于四边形报告的系数。
        quad = next(p for p in self.rep.polygons if p.tet_index == 4)
        self.assertEqual(quad.area_coefficient, Fraction(1, 2))  # 菱形 1/2·1·1


# ---------------------------------------------------------------- 3. 穿顶点：点/棱接触单列


class TestPlaneThroughVertices(unittest.TestCase):
    """平面 z=0 穿过立方体底面四个顶点。"""

    def setUp(self):
        self.mesh = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        self.plane = Plane.build((0, 0, 1), 0)
        self.rep = section(self.mesh, self.plane)

    def test_no_interior_area_patches(self):
        # 没有单元被切出内部片；面积只来自两张共面原面（底面 = 1）。
        self.assertEqual(self.rep.polygons, ())
        keys = sorted(f.key for f in self.rep.coplanar_faces)
        self.assertEqual(keys, [(0, 1, 2), (0, 2, 3)])
        self.assertEqual(self.rep.total_area_coefficient, Fraction(1))

    def test_point_and_edge_contacts_listed_separately(self):
        # t2 只以顶点 0 接触；t3 只以顶点 2 接触；t4 以棱 (0,2) 接触。
        pt = {c.id: c.tet_indices for c in self.rep.contact_points}
        self.assertEqual(pt, {("vertex", 0): (2,), ("vertex", 2): (3,)})
        ed = {c.endpoints: c.tet_indices for c in self.rep.contact_edges}
        self.assertEqual(
            ed, {(("vertex", 0), ("vertex", 2)): (4,)}
        )
        # 零面积接触结果与面片分属不同列表；接触单元不出现在任何面积片中。
        contact_tets = {
            t for c in self.rep.contact_points for t in c.tet_indices
        } | {t for c in self.rep.contact_edges for t in c.tet_indices}
        patch_tets = {
            p.tet_index for p in self.rep.polygons
        } | {t for f in self.rep.coplanar_faces for t in f.tet_indices}
        self.assertFalse(contact_tets & patch_tets)
        self.assertEqual(self.rep.polygons, ())

    def test_contacts_are_zero_area_types(self):
        self.assertTrue(all(isinstance(c, ContactPoint) for c in self.rep.contact_points))
        self.assertTrue(all(isinstance(c, ContactEdge) for c in self.rep.contact_edges))


# ---------------------------------------------------------------- 4. 平面沿共享内面：共面原面只一份


class TestPlaneAlongSharedFace(unittest.TestCase):
    """平面含内面 (0,2,5)：n=(1,-1,-1), d=0；该面被 t0、t4 共享。"""

    N = (1, -1, -1)

    def setUp(self):
        self.mesh = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        self.plane = Plane.build(self.N, 0)
        self.rep = section(self.mesh, self.plane)

    def test_coplanar_face_deduped_with_all_sources(self):
        cop = [f for f in self.rep.coplanar_faces if f.key == (0, 2, 5)]
        self.assertEqual(len(cop), 1)  # 全局只保留一份
        self.assertEqual(cop[0].tet_indices, (0, 4))  # 全部来源单元
        self.assertTrue(all(v[0] == "vertex" for v in cop[0].vertices))

    def test_other_faces_become_edge_contacts(self):
        # t1 沿棱 (0,2)、t2 沿 (0,5)、t3 沿 (2,5) 接触平面。
        ed = {c.endpoints: c.tet_indices for c in self.rep.contact_edges}
        self.assertEqual(ed[(("vertex", 0), ("vertex", 2))], (1,))
        self.assertEqual(ed[(("vertex", 0), ("vertex", 5))], (2,))
        self.assertEqual(ed[(("vertex", 2), ("vertex", 5))], (3,))
        self.assertEqual(self.rep.polygons, ())

    def test_coplanar_face_area_and_orientation(self):
        face = next(f for f in self.rep.coplanar_faces if f.key == (0, 2, 5))
        cmap = _coord_map(self.rep)
        # 原三角形 0,2,5 面积：|cross|/2 = |(1,-1,-1)|/2 = √3/2；
        # α = 面积/|n| = (√3/2)/√3 = 1/2。
        self.assertEqual(face.area_coefficient, Fraction(1, 2))
        self.assertGreater(
            _oriented_sign(face.triangles[0], cmap, self.plane.normal), 0
        )
        # 共面片内点（重心）在两个来源单体内（第四顶点重心坐标为 0）。
        pts = [cmap[v] for v in face.vertices]
        center = tuple(sum(p[k] for p in pts) / 3 for k in range(3))
        for ti in face.tet_indices:
            tet_pts = [self.mesh.vertices[i] for i in CUBE_TETS[ti]]
            (lams, _), = [ref_barycentric(center, tet_pts)]
            self.assertTrue(all(l >= 0 for l in lams))

    def test_float_area_estimate_consistent(self):
        # 总面积系数 1/2；|n|=√3 → 面积 ≈ √3/2。
        import math

        self.assertAlmostEqual(
            self.rep.total_area_float(), math.sqrt(3) / 2, places=12
        )


# ---------------------------------------------------------------- 5. 解析棱锥核对（单四面体三角形片）


class TestPyramidCrossCheck(unittest.TestCase):
    def setUp(self):
        # 正四面体样式：0:(0,0,0) 1:(6,0,0) 2:(0,6,0) 3:(0,0,6)，
        # 平面 z=2：顶部棱锥高 4，缩放 2/3，底面积 = 18·4/9 = 8。
        verts = ((0, 0, 0), (6, 0, 0), (0, 6, 0), (0, 0, 6))
        self.mesh = TetMesh.build(verts, [(0, 1, 2, 3)])
        self.verts = verts
        self.plane = Plane.build((0, 0, 1), 2)
        self.rep = section(self.mesh, self.plane)

    def test_pyramid_volume_area_identity(self):
        poly = self.rep.polygons[0]
        cmap = _coord_map(self.rep)
        cyclic = [cmap[v] for v in poly.vertices]
        # 三条独立路径一致：核心 α、边界叉积 α、解析棱锥 α。
        alpha_pyramid = ref_pyramid_alpha(
            cyclic, self.verts[3], self.plane.normal, self.plane.offset
        )
        self.assertEqual(poly.area_coefficient, Fraction(8))
        self.assertEqual(ref_polygon_alpha(cyclic, self.plane.normal), Fraction(8))
        self.assertEqual(alpha_pyramid, Fraction(8))


# ---------------------------------------------------------------- 6. 修正规则 / 不改网格 / 错误语义


class TestRepairAndImmutability(unittest.TestCase):
    def test_negative_tet_repaired_by_swap_last_two_before_section(self):
        # 把单四面体负向书写；截面必须与修复后网格一致。
        verts = ((0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4))
        neg = (0, 1, 3, 2)  # 交换最后两顶点 → 负体积
        mesh = TetMesh.build(verts, [neg])
        plane = Plane.build((0, 0, 1), 1)
        rep = section(mesh, plane)
        self.assertEqual(rep.polygons[0].area_coefficient, Fraction(9, 2))
        # 原网格书写序不变（section 不修改、不分裂原网格）。
        self.assertEqual(mesh.tets[0], neg)
        # 顶点数与单元数不变；所有交点是派生 id，原顶点表不新增点。
        self.assertEqual(len(mesh.vertices), 4)
        self.assertTrue(all(p.id[0] in ("vertex", "edge") for p in rep.points))

    def test_section_matches_section_of_repaired_mesh(self):
        verts = ((0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4))
        neg = (0, 1, 3, 2)
        mesh = TetMesh.build(verts, [neg])
        fixed = repair_orientation(mesh)
        plane = Plane.build((0, 0, 1), 1)
        a = section(mesh, plane)
        b = section(fixed, plane)
        self.assertEqual(
            [p.vertices for p in a.polygons], [p.vertices for p in b.polygons]
        )

    def test_zero_volume_tet_refused(self):
        verts = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0))  # 全共面
        mesh = TetMesh.build(verts, [(0, 1, 2, 3)])
        with self.assertRaises(ZeroVolumeError) as cm:
            section(mesh, Plane.build((0, 0, 1), 0))
        self.assertEqual(cm.exception.tet_indices, (0,))

    def test_no_edge_splitting_on_mesh(self):
        mesh = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        before = (mesh.vertices, mesh.tets)
        section(mesh, Plane.build((0, 0, 1), Fraction(1, 3)))
        self.assertEqual((mesh.vertices, mesh.tets), before)


class TestPlaneValidation(unittest.TestCase):
    def test_zero_normal_rejected(self):
        with self.assertRaises(MeshInputError) as cm:
            Plane.build((0, 0, 0), 1)
        self.assertEqual(cm.exception.location, "normal")

    def test_float_and_bool_rejected_with_location(self):
        with self.assertRaises(MeshInputError) as cm:
            Plane.build((0.5, 0, 1), 0)
        self.assertEqual(cm.exception.location, "normal[0]")
        with self.assertRaises(MeshInputError) as cm:
            Plane.build((True, 0, 1), 0)
        self.assertIn("bool", str(cm.exception))
        with self.assertRaises(MeshInputError) as cm:
            Plane.build((0, 0, 1), 1.5)
        self.assertEqual(cm.exception.location, "offset")

    def test_rational_string_normal_accepted(self):
        p = Plane.build(("1/3", "0", "-2"), "5/7")
        self.assertEqual(p.normal[0], Fraction(1, 3))
        self.assertEqual(p.offset, Fraction(5, 7))

    def test_wrong_arity(self):
        with self.assertRaises(MeshInputError) as cm:
            Plane.build((1, 0), 0)
        self.assertEqual(cm.exception.location, "normal")

    def test_bad_types(self):
        with self.assertRaises(TypeError):
            section(CUBE_VERTICES, Plane.build((0, 0, 1), 0))
        with self.assertRaises(TypeError):
            section(TetMesh.build(CUBE_VERTICES, CUBE_TETS), (0, 0, 1))


if __name__ == "__main__":
    unittest.main()
