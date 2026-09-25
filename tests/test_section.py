"""平面截面必要测试：共享棱身份、穿顶点、沿共享面、四边形截面、
解析棱锥面积核对、平面方程与来源单元包含验证、负向单元修正、
零面积接触、退化单元、非法平面参数与输入不可变。"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    MeshInputError,
    TetMesh,
    section,
)

from reference import ref_barycentric

# 两个共享面 (0,1,2) 的四面体：顶点 3 在 z>0 侧，顶点 4 在 z<0 侧。
TWO_VERTS = (
    (0, 0, 0),  # 0
    (1, 0, 0),  # 1
    (0, 1, 0),  # 2
    (0, 0, 1),  # 3
    (0, 0, -1),  # 4
)
TWO_TETS = (
    (0, 1, 2, 3),  # 正体积
    (0, 2, 1, 4),  # 正体积（与上者共享面 (0,1,2)）
)

# 方底棱锥：底面 z=0 的 [-1,1]^2，顶点 (0,0,1)，两个正体积四面体。
PYRAMID_VERTS = (
    (-1, -1, 0),  # 0
    (1, -1, 0),  # 1
    (1, 1, 0),  # 2
    (-1, 1, 0),  # 3
    (0, 0, 1),  # 4
)
PYRAMID_TETS = (
    (0, 1, 2, 4),
    (0, 2, 3, 4),
)


def assert_all_on_plane(tc, report):
    """每个截面点精确满足平面方程 normal·p = offset。"""
    for p in report.points:
        val = sum(report.normal[k] * p.coords[k] for k in range(3))
        tc.assertEqual(val, report.offset, f"点 {p.origin} 不在平面上")
    tc.assertEqual(report.plane_residuals(), (Fraction(0),) * len(report.points))


def assert_patch_centroids_in_sources(tc, mesh, report):
    """每片面片的多边形质心（有理）落在每个来源单元内（含边界）。"""
    for patch in report.patches:
        pts = [report.points[pid].coords for pid in patch.polygon]
        centroid = tuple(
            sum(p[k] for p in pts) / len(pts) for k in range(3)
        )
        for ti in patch.source_tets:
            bary = ref_barycentric(mesh.vertices, mesh.tets[ti], centroid)
            tc.assertIsNotNone(bary)
            for lam in bary:
                tc.assertGreaterEqual(lam, Fraction(0))
                tc.assertLessEqual(lam, Fraction(1))


class TestTwoTetsSharedEdgeIdentity(unittest.TestCase):
    """两个共享面的四面体被同一平面截：共享棱交点身份一致。"""

    def setUp(self):
        self.mesh = TetMesh.build(TWO_VERTS, TWO_TETS)
        self.report = section(self.mesh, (1, 0, 0), Fraction(1, 3))

    def test_two_triangle_patches(self):
        kinds = [p.kind for p in self.report.patches]
        self.assertEqual(kinds, ["triangle", "triangle"])
        self.assertEqual(
            [p.source_tets for p in self.report.patches], [(0,), (1,)]
        )

    def test_shared_edge_points_have_same_global_ids(self):
        # 两片共享原网格棱 (0,1) 与 (1,2) 上的交点：全局 ID 必须相同，
        # 身份就是原顶点对，不靠坐标近似焊接。
        poly0, poly1 = (set(p.polygon) for p in self.report.patches)
        shared = poly0 & poly1
        self.assertEqual(len(shared), 2)
        origins = sorted(self.report.points[pid].origin for pid in shared)
        self.assertEqual(origins, [("edge", 0, 1), ("edge", 1, 2)])
        # 每片恰有 3 个不同的点，合并后共 4 个全局点。
        self.assertEqual(len(self.report.points), 4)
        self.assertTrue(all(o[0] == "edge" for o in
                            (p.origin for p in self.report.points)))

    def test_exact_coordinates_and_area(self):
        coords = {p.origin: p.coords for p in self.report.points}
        self.assertEqual(
            coords[("edge", 0, 1)],
            (Fraction(1, 3), Fraction(0), Fraction(0)),
        )
        self.assertEqual(
            coords[("edge", 1, 2)],
            (Fraction(1, 3), Fraction(2, 3), Fraction(0)),
        )
        # 直角边 2/3 的三角形：精确面积平方 = (2·(2/9))²/4 = 4/81。
        for patch in self.report.patches:
            self.assertEqual(patch.area2, Fraction(4, 81))

    def test_orientation_follows_plane_normal(self):
        # 所有输出三角形的叉积法向与 (1,0,0) 点积为正。
        for patch in self.report.patches:
            for tri in patch.triangles:
                a, b, c = (self.report.points[p].coords for p in tri)
                u = tuple(b[k] - a[k] for k in range(3))
                v = tuple(c[k] - a[k] for k in range(3))
                nx = u[1] * v[2] - u[2] * v[1]
                self.assertGreater(nx, 0)

    def test_points_on_plane_and_centroids_inside(self):
        assert_all_on_plane(self, self.report)
        assert_patch_centroids_in_sources(self, self.mesh, self.report)


class TestPlaneThroughVertex(unittest.TestCase):
    """平面穿过原顶点：统一分类为“在平面上”，身份是原顶点本身。"""

    def setUp(self):
        self.mesh = TetMesh.build(TWO_VERTS[:4], [(0, 1, 2, 3)])

    def test_section_through_two_vertices(self):
        # 平面 2x - y = 0 过顶点 0 和 3，并切棱 (1,2)。
        report = section(self.mesh, (2, -1, 0), 0)
        self.assertEqual(len(report.patches), 1)
        patch = report.patches[0]
        origins = sorted(
            report.points[pid].origin for pid in patch.polygon
        )
        self.assertEqual(
            origins, [("edge", 1, 2), ("vertex", 0), ("vertex", 3)]
        )
        # 顶点身份点的坐标就是原顶点坐标对象（不新建、不近似）。
        for p in report.points:
            if p.origin[0] == "vertex":
                self.assertIs(
                    p.coords, self.mesh.vertices[p.origin[1]]
                )
        assert_all_on_plane(self, report)
        assert_patch_centroids_in_sources(self, self.mesh, report)

    def test_touch_only_at_vertex_is_zero_area_contact(self):
        # 平面 x+y+z = 0 只在顶点 0 处碰到单元。
        report = section(self.mesh, (1, 1, 1), 0)
        self.assertEqual(report.patches, ())
        self.assertEqual(report.edge_contacts, ())
        self.assertEqual(len(report.point_contacts), 1)
        contact = report.point_contacts[0]
        self.assertEqual(contact.vertex, 0)
        self.assertEqual(contact.source_tets, (0,))
        self.assertEqual(
            report.points[contact.point_id].origin, ("vertex", 0)
        )

    def test_touch_only_along_edge_is_zero_area_contact(self):
        # 平面 x + y = 0 含棱 (0,3)，另两个顶点严格同侧。
        report = section(self.mesh, (1, 1, 0), 0)
        self.assertEqual(report.patches, ())
        self.assertEqual(report.point_contacts, ())
        self.assertEqual(len(report.edge_contacts), 1)
        contact = report.edge_contacts[0]
        self.assertEqual(contact.edge, (0, 3))
        self.assertEqual(contact.source_tets, (0,))
        self.assertEqual(len(report.points), 2)
        assert_all_on_plane(self, report)


class TestCoplanarSharedFace(unittest.TestCase):
    """平面与共享面重合：完全共面的原面只保留一份，记录所有来源。"""

    def setUp(self):
        self.mesh = TetMesh.build(TWO_VERTS, TWO_TETS)
        self.report = section(self.mesh, (0, 0, 1), 0)

    def test_single_coplanar_patch_with_all_sources(self):
        self.assertEqual(len(self.report.patches), 1)
        patch = self.report.patches[0]
        self.assertEqual(patch.kind, "coplanar_face")
        self.assertEqual(patch.source_tets, (0, 1))
        self.assertEqual(len(patch.polygon), 3)
        self.assertEqual(len(patch.triangles), 1)
        # 边界点全部是原顶点身份。
        origins = sorted(
            self.report.points[pid].origin for pid in patch.polygon
        )
        self.assertEqual(
            origins, [("vertex", 0), ("vertex", 1), ("vertex", 2)]
        )
        # 三角形 (0,0,0),(1,0,0),(0,1,0)：面积 1/2，精确面积平方 1/4。
        self.assertEqual(patch.area2, Fraction(1, 4))
        self.assertEqual(self.report.point_contacts, ())
        self.assertEqual(self.report.edge_contacts, ())

    def test_coplanar_patch_on_plane_and_inside_both_tets(self):
        assert_all_on_plane(self, self.report)
        assert_patch_centroids_in_sources(self, self.mesh, self.report)


class TestQuadSection(unittest.TestCase):
    """2-2 分布给四边形：按全局交点 ID 选确定对角线。"""

    def setUp(self):
        self.mesh = TetMesh.build(TWO_VERTS[:4], [(0, 1, 2, 3)])
        # 平面 x + z = 1/2：顶点 1、3 在正侧，0、2 在负侧。
        self.report = section(self.mesh, (1, 0, 1), Fraction(1, 2))

    def test_quad_split_into_two_triangles(self):
        self.assertEqual(len(self.report.patches), 1)
        patch = self.report.patches[0]
        self.assertEqual(patch.kind, "quad")
        self.assertEqual(len(patch.polygon), 4)
        self.assertEqual(len(patch.triangles), 2)
        self.assertEqual(len(self.report.points), 4)
        self.assertTrue(
            all(p.origin[0] == "edge" for p in self.report.points)
        )

    def test_diagonal_contains_min_global_id(self):
        patch = self.report.patches[0]
        t1, t2 = (set(t) for t in patch.triangles)
        diagonal = t1 & t2
        self.assertEqual(len(diagonal), 2)
        # 契约：对角线含四边形四个全局 ID 中的最小者。
        self.assertIn(min(patch.polygon), diagonal)
        # 两个三角形合起来恰好覆盖多边形边界点。
        self.assertEqual(t1 | t2, set(patch.polygon))

    def test_exact_area_and_geometry(self):
        # 截面是矩形 (1/√2)·(1/2)，精确面积平方 = 1/8。
        patch = self.report.patches[0]
        self.assertEqual(patch.area2, Fraction(1, 8))
        assert_all_on_plane(self, self.report)
        assert_patch_centroids_in_sources(self, self.mesh, self.report)


class TestPyramidAnalyticArea(unittest.TestCase):
    """解析棱锥核对：z = 1/2 截方底棱锥，截面是边长 1 的正方形。"""

    def setUp(self):
        self.mesh = TetMesh.build(PYRAMID_VERTS, PYRAMID_TETS)
        self.report = section(self.mesh, (0, 0, 1), Fraction(1, 2))

    def test_total_area_matches_analytic_square(self):
        # 高 h 处截面边长 2(1-h)：h=1/2 时边长 1，面积 1。
        self.assertEqual(len(self.report.patches), 2)
        for patch in self.report.patches:
            self.assertEqual(patch.kind, "triangle")
            # 每个四面体贡献直角边 1 的三角形：精确面积平方 1/4。
            self.assertEqual(patch.area2, Fraction(1, 4))
        self.assertEqual(self.report.total_area(), 1.0)

    def test_points_on_plane_and_inside_source_tets(self):
        assert_all_on_plane(self, self.report)
        assert_patch_centroids_in_sources(self, self.mesh, self.report)


class TestNegativeAndDegenerateTets(unittest.TestCase):
    @staticmethod
    def _canon(cycle):
        # 环序规范化：同一朝向的循环旋转视为同一多边形/三角形。
        m = cycle.index(min(cycle))
        return tuple(cycle[m:] + cycle[:m])

    def test_negative_tet_repaired_by_original_rule_before_section(self):
        # 单元 1 写成负体积排列；截面结果必须与正向写法一致。
        neg_mesh = TetMesh.build(TWO_VERTS, [(0, 1, 2, 3), (0, 1, 2, 4)])
        pos_mesh = TetMesh.build(TWO_VERTS, TWO_TETS)
        neg_report = section(neg_mesh, (1, 0, 0), Fraction(1, 3))
        pos_report = section(pos_mesh, (1, 0, 0), Fraction(1, 3))
        self.assertEqual(neg_report.repaired_tets, (1,))
        self.assertEqual(pos_report.repaired_tets, ())
        self.assertEqual(neg_report.points, pos_report.points)
        # 修正规则只交换最后两顶点，环序可能差一个循环旋转，
        # 规范化后多边形、三角剖分、来源与精确面积必须完全一致。
        def canon_patches(report):
            return sorted(
                (
                    p.kind,
                    self._canon(list(p.polygon)),
                    tuple(sorted(self._canon(list(t)) for t in p.triangles)),
                    p.source_tets,
                    p.area2,
                )
                for p in report.patches
            )

        self.assertEqual(canon_patches(neg_report), canon_patches(pos_report))
        # 输入网格不被修改：负向排列原样保留。
        self.assertEqual(neg_mesh.tets[1], (0, 1, 2, 4))

    def test_degenerate_tet_skipped_and_listed(self):
        verts = (*TWO_VERTS, (1, 1, 0))  # 5 与 0,1,2 共面
        mesh = TetMesh.build(verts, [*TWO_TETS, (0, 1, 2, 5)])
        report = section(mesh, (1, 0, 0), Fraction(1, 3))
        self.assertEqual(report.degenerate_tets, (2,))
        self.assertEqual(len(report.patches), 2)
        self.assertTrue(
            all(2 not in p.source_tets for p in report.patches)
        )


class TestSectionInputValidation(unittest.TestCase):
    def setUp(self):
        self.mesh = TetMesh.build(TWO_VERTS[:4], [(0, 1, 2, 3)])

    def test_rejects_zero_normal(self):
        with self.assertRaises(MeshInputError) as cm:
            section(self.mesh, (0, 0, 0), 1)
        self.assertEqual(cm.exception.location, "plane.normal")

    def test_rejects_float_and_bad_rational(self):
        with self.assertRaises(MeshInputError) as cm:
            section(self.mesh, (0, 0, 0.5), 1)
        self.assertEqual(cm.exception.location, "plane.normal[2]")
        with self.assertRaises(MeshInputError) as cm:
            section(self.mesh, (0, 0, 1), "1/0")
        self.assertEqual(cm.exception.location, "plane.offset")
        with self.assertRaises(MeshInputError) as cm:
            section(self.mesh, (0, 0), 1)
        self.assertEqual(cm.exception.location, "plane.normal")

    def test_rejects_non_mesh(self):
        with self.assertRaises(TypeError):
            section("not a mesh", (0, 0, 1), 0)

    def test_rational_string_and_int_params_accepted(self):
        report = section(self.mesh, ("0", "0", "1"), "1/2")
        self.assertEqual(report.normal, (Fraction(0), Fraction(0), Fraction(1)))
        self.assertEqual(report.offset, Fraction(1, 2))
        self.assertEqual(len(report.patches), 1)


class TestInputMeshUnchanged(unittest.TestCase):
    def test_section_does_not_mutate_mesh_or_split_edges(self):
        mesh = TetMesh.build(PYRAMID_VERTS, PYRAMID_TETS)
        tets_before = mesh.tets
        verts_before = mesh.vertices
        report = section(mesh, (0, 0, 1), Fraction(1, 2))
        # 单元与顶点完全不动（对象恒等），没有边分裂、没有新顶点入网格。
        self.assertIs(mesh.tets, tets_before)
        self.assertIs(mesh.vertices, verts_before)
        self.assertEqual(len(mesh.vertices), 5)
        # 截面点是独立的新对象，不回写网格。
        self.assertTrue(
            all(p.coords not in mesh.vertices or p.origin[0] == "vertex"
                for p in report.points)
        )


if __name__ == "__main__":
    unittest.main()
