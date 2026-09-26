"""守恒重映射测试：手算核对 + 守恒性 + 缺口保留 + 非法输入。

所有期望值均为手算精确 Fraction，不调用被测核心得到。
"""

from __future__ import annotations

import unittest
from fractions import Fraction as F

from orient_audit import (
    MeshInputError,
    RemapResult,
    TetMesh,
    intersection_volume,
    remap,
    tet_volume,
)

from fixtures import CUBE_TETS, CUBE_VERTICES

# 单位立方体的 6 四面体剖分（绕主对角线 (0,0,0)-(1,1,1)），每单元体积 1/6。
CUBE6_TETS = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)


def box_vertices(x0, x1, y0=0, y1=1, z0=0, z1=1):
    """长方体 8 顶点，编号与 CUBE_VERTICES 相同。"""

    return (
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    )


class IntersectionVolumeTest(unittest.TestCase):
    def test_self_intersection_equals_tet_volume(self):
        tet = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        mesh = TetMesh.build([tet[0], tet[1], tet[2], tet[3]], [(0, 1, 2, 3)])
        verts = tuple(mesh.vertices[k] for k in mesh.tets[0])
        self.assertEqual(tet_volume(verts), F(1, 6))
        self.assertEqual(intersection_volume(verts, verts), F(1, 6))

    def test_inner_tet_fully_contained(self):
        outer = ((0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2))
        inner = ((0, 0, 0), (F(1, 2), 0, 0), (0, F(1, 2), 0), (0, 0, F(1, 2)))
        # 内四面体体积 (1/2)^3 / 6 = 1/48
        self.assertEqual(intersection_volume(outer, inner), F(1, 48))

    def test_disjoint_is_zero(self):
        a = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        b = ((2, 0, 0), (3, 0, 0), (2, 1, 0), (2, 0, 1))
        self.assertEqual(intersection_volume(a, b), F(0))

    def test_partial_overlap_hand_computed(self):
        # A: 角四面体 x,y,z>=0, x+y+z<=1，体积 1/6。
        a = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        # C: 边长 1/2 的平行四面体，三个顶点越出平面 x+y+z=1；
        # 交集是边长 1/4 的小四面体，体积 (1/4)^3/6 = 1/384。
        c = (
            (F(1, 4), F(1, 4), F(1, 4)),
            (F(3, 4), F(1, 4), F(1, 4)),
            (F(1, 4), F(3, 4), F(1, 4)),
            (F(1, 4), F(1, 4), F(3, 4)),
        )
        self.assertEqual(intersection_volume(a, c), F(1, 384))

    def test_orientation_irrelevant(self):
        a = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        rev = (a[1], a[0], a[3], a[2])
        self.assertEqual(intersection_volume(a, rev), F(1, 6))


class RemapConstantFieldTest(unittest.TestCase):
    """同一立方体的两种四面体剖分（5 剖分 -> 6 剖分），常量场。"""

    def setUp(self):
        self.source = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        self.target = TetMesh.build(CUBE_VERTICES, CUBE6_TETS)
        self.result = remap(self.source, self.target, [F(7, 2)] * 5)

    def test_constant_field_preserved(self):
        for avg in self.result.target_average_density:
            self.assertEqual(avg, F(7, 2))

    def test_full_coverage_no_gap(self):
        for gap in self.result.target_uncovered_volume:
            self.assertEqual(gap, F(0))
        for cov in self.result.target_covered_volume:
            self.assertEqual(cov, F(1, 6))

    def test_mass_conserved_exactly(self):
        # 源总质量 = 密度 * 立方体体积 = 7/2
        self.assertEqual(self.result.total_target_mass, F(7, 2))
        self.assertEqual(self.result.covered_source_mass, F(7, 2))
        self.assertEqual(self.result.mass_error, F(0))

    def test_weights_sparse_and_cover_target(self):
        seen = {(w.source_index, w.target_index) for w in self.result.weights}
        self.assertEqual(len(seen), len(self.result.weights))
        for w in self.result.weights:
            self.assertGreater(w.volume, 0)
        # 每个目标单元的权重和 = 1/6
        for ti in range(6):
            s = sum(
                (w.volume for w in self.result.weights if w.target_index == ti),
                F(0),
            )
            self.assertEqual(s, F(1, 6))

    def test_inputs_not_modified(self):
        verts_s = self.source.vertices
        tets_s = self.source.tets
        remap(self.source, self.target, [1] * 5)
        self.assertIs(self.source.vertices, verts_s)
        self.assertIs(self.source.tets, tets_s)


class RemapPiecewiseTest(unittest.TestCase):
    """跨界分段常量：源为 [0,1/2] 与 [1/2,1] 两个半立方体（各 5 剖分），
    密度 1 与 3；目标为中间条 [1/4,3/4]x[0,1]x[0,1] 的 5 剖分。
    总质量手算 = 1/4*1 + 1/4*3 = 1；单个跨界四面体另有精确手算。"""

    def setUp(self):
        left = TetMesh.build(box_vertices(0, F(1, 2)), CUBE_TETS)
        right = TetMesh.build(box_vertices(F(1, 2), 1), CUBE_TETS)
        self.source = TetMesh.build(
            list(left.vertices) + list(right.vertices),
            list(left.tets) + [tuple(k + 8 for k in t) for t in right.tets],
        )
        self.target = TetMesh.build(box_vertices(F(1, 4), F(3, 4)), CUBE_TETS)
        self.result = remap(self.source, self.target, [1] * 5 + [3] * 5)

    def test_piecewise_average_density(self):
        # 5 剖分关于 x=1/2 不对称，逐单元手算不可行；核验不变量：
        # 每单元被满覆盖、平均密度介于两侧源密度之间。
        for ti, avg in enumerate(self.result.target_average_density):
            self.assertEqual(self.result.target_uncovered_volume[ti], F(0))
            self.assertGreaterEqual(avg, F(1))
            self.assertLessEqual(avg, F(3))

    def test_piecewise_mass(self):
        # 目标总体积 1/2，左半密度 1 右半密度 3：质量 = 1/4*1 + 1/4*3 = 1
        self.assertEqual(self.result.total_target_mass, F(1))
        self.assertEqual(self.result.mass_error, F(0))
        for gap in self.result.target_uncovered_volume:
            self.assertEqual(gap, F(0))

    def test_single_crossing_tet_hand_check(self):
        # 单个跨界目标四面体，手算：左部体积 1/24（密度1），右部 1/24（密度3）。
        target = TetMesh.build(
            [(F(1, 4), 0, 0), (F(3, 4), 0, 0), (F(1, 4), 1, 0), (F(1, 4), 0, 1)],
            [(0, 1, 2, 3)],
        )
        result = remap(self.source, target, [1] * 5 + [3] * 5)
        # 手算：平面 x=1/2 把该四面体切成左侧五面体（体积 7/96，密度 1）
        # 与右侧小四面体（体积 1/96，密度 3），总质量 7/96 + 3/96 = 5/48。
        self.assertEqual(result.target_mass[0], F(7, 96) * 1 + F(1, 96) * 3)
        # 单元总体积 1/12，平均密度 = (5/48) / (1/12) = 5/4。
        self.assertEqual(result.target_average_density[0], F(5, 4))
        self.assertEqual(result.target_uncovered_volume[0], F(0))


class RemapGapTest(unittest.TestCase):
    """目标网格大于源区域：缺口保留，不归一到满覆盖。"""

    def setUp(self):
        self.source = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        # 目标：[0,1]x[0,1]x[0,2]，上半完全无覆盖。
        self.target = TetMesh.build(box_vertices(0, 1, 0, 1, 0, 2), CUBE_TETS)
        self.result = remap(self.source, self.target, [2] * 5)

    def test_uncovered_volume_kept(self):
        total_gap = sum(self.result.target_uncovered_volume, F(0))
        self.assertEqual(total_gap, F(1))  # 上半单位立方体完全未覆盖
        self.assertGreater(max(self.result.target_uncovered_volume), F(0))

    def test_no_renormalization(self):
        # 拉伸后角单元体积 1/3、中心单元 2/3；平均密度 = 覆盖质量 /
        # 整个单元体积，覆盖部分密度 2，未覆盖部分按 0，不归一。
        cell_volumes = (F(1, 3), F(1, 3), F(1, 3), F(1, 3), F(2, 3))
        for ti, avg in enumerate(self.result.target_average_density):
            cov = self.result.target_covered_volume[ti]
            self.assertEqual(avg, 2 * cov / cell_volumes[ti])
        # 总质量仍是源总质量 2，没有被放大到满覆盖。
        self.assertEqual(self.result.total_target_mass, F(2))
        self.assertEqual(self.result.mass_error, F(0))


class RemapInputErrorTest(unittest.TestCase):
    def setUp(self):
        self.mesh = TetMesh.build(CUBE_VERTICES, CUBE_TETS)

    def test_density_length_mismatch(self):
        with self.assertRaises(MeshInputError):
            remap(self.mesh, self.mesh, [1, 2])

    def test_float_density_rejected(self):
        with self.assertRaises(MeshInputError):
            remap(self.mesh, self.mesh, [0.5] * 5)

    def test_bool_density_rejected(self):
        with self.assertRaises(MeshInputError):
            remap(self.mesh, self.mesh, [True] * 5)

    def test_zero_volume_tet_rejected(self):
        degenerate = TetMesh.build(
            [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)], [(0, 1, 2, 3)]
        )
        with self.assertRaises(MeshInputError):
            remap(degenerate, self.mesh, [1])
        with self.assertRaises(MeshInputError):
            remap(self.mesh, degenerate, [1] * 5)

    def test_string_density_accepted(self):
        result = remap(self.mesh, self.mesh, ["5/2"] * 5)
        self.assertEqual(result.total_target_mass, F(5, 2))


if __name__ == "__main__":
    unittest.main()
