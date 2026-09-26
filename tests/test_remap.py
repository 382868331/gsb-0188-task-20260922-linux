"""守恒重映射的手算核验测试。

参考值全部由解析几何手算给出，不调用被测核心得到期望值。
"""

from __future__ import annotations

import unittest
from fractions import Fraction as F

from orient_audit import MeshInputError, TetMesh
from orient_audit.remap import RemapReport, remap, tet_tet_intersection_volume


def cube_vertices(scale=F(1), shift=(F(0), F(0), F(0))):
    return [
        (F(x) * scale + shift[0], F(y) * scale + shift[1], F(z) * scale + shift[2])
        for x in (0, 1) for y in (0, 1) for z in (0, 1)
    ]


def idx(x, y, z):
    return x * 4 + y * 2 + z


# 单位立方体的 6 四面体剖分（主对角线 000-111）
T6 = [
    (idx(0, 0, 0), idx(1, 0, 0), idx(1, 1, 0), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(1, 0, 0), idx(1, 0, 1), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 1, 0), idx(1, 1, 0), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 1, 0), idx(0, 1, 1), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 0, 1), idx(1, 0, 1), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 0, 1), idx(0, 1, 1), idx(1, 1, 1)),
]

# 同一立方体的 5 四面体剖分：4 个角四面体 + 1 个中心四面体
T5 = [
    (idx(0, 0, 0), idx(1, 0, 0), idx(0, 1, 0), idx(0, 0, 1)),
    (idx(1, 1, 0), idx(1, 0, 0), idx(0, 1, 0), idx(1, 1, 1)),
    (idx(1, 0, 1), idx(1, 0, 0), idx(0, 0, 1), idx(1, 1, 1)),
    (idx(0, 1, 1), idx(0, 1, 0), idx(0, 0, 1), idx(1, 1, 1)),
    (idx(1, 0, 0), idx(0, 1, 0), idx(0, 0, 1), idx(1, 1, 1)),
]


def box_tets(x0, x1):
    """[x0,x1]x[0,1]x[0,1] 长方体的 6 四面体剖分（顶点坐标直接给出）。"""
    verts = [
        (F(x), F(y), F(z)) for x in (x0, x1) for y in (0, 1) for z in (0, 1)
    ]
    return verts, T6


class IntersectionVolumeTest(unittest.TestCase):
    def setUp(self):
        self.V = cube_vertices()

    def test_disjoint(self):
        V2 = cube_vertices(shift=(F(2), F(0), F(0)))
        for a in T6:
            for b in T6:
                self.assertEqual(
                    tet_tet_intersection_volume(self.V, a, V2, b), F(0)
                )

    def test_self_overlap_is_partition(self):
        # 同一剖分两两交叠：总和恰为立方体体积 1
        total = sum(
            tet_tet_intersection_volume(self.V, a, self.V, b) for a in T6 for b in T6
        )
        self.assertEqual(total, F(1))

    def test_shifted_cube_by_half(self):
        # 平移 (1/2,1/2,1/2) 的立方体：交叠为 [1/2,1]^3，体积 1/8
        V2 = cube_vertices(shift=(F(1, 2), F(1, 2), F(1, 2)))
        total = sum(
            tet_tet_intersection_volume(self.V, a, V2, b) for a in T6 for b in T6
        )
        self.assertEqual(total, F(1, 8))

    def test_contained_cube(self):
        # [0,1/3]^3 完全在内部：交叠 = 1/27
        V3 = cube_vertices(scale=F(1, 3))
        total = sum(
            tet_tet_intersection_volume(self.V, a, V3, b) for a in T6 for b in T6
        )
        self.assertEqual(total, F(1, 27))

    def test_symmetry(self):
        V2 = cube_vertices(shift=(F(1, 3), F(1, 2), F(0)))
        for a in T6:
            for b in T6:
                self.assertEqual(
                    tet_tet_intersection_volume(self.V, a, V2, b),
                    tet_tet_intersection_volume(V2, b, self.V, a),
                )


class RemapConstantFieldTest(unittest.TestCase):
    """同一立方体的两种不同四面体剖分之间转常量场。"""

    def setUp(self):
        V = cube_vertices()
        self.src = TetMesh.build(V, T6)
        self.tgt = TetMesh.build(V, T5)

    def test_constant_field_exact(self):
        rep = remap(self.src, [F(7, 3)] * 6, self.tgt)
        self.assertEqual(rep.mass_error, 0)
        self.assertEqual(sum(rep.uncovered_volumes), 0)
        # 每个目标单元平均密度恰为常量
        for avg in rep.target_averages:
            self.assertEqual(avg, F(7, 3))
        # 总质量 = 密度 * 体积
        self.assertEqual(rep.target_mass, F(7, 3))
        self.assertEqual(rep.source_mass, F(7, 3))

    def test_conservation_across_tessellations(self):
        rho = [F(0), F(1), F(2), F(3), F(5), F(8)]
        rep = remap(self.src, rho, self.tgt)
        self.assertEqual(rep.mass_error, 0)
        self.assertEqual(rep.target_mass, rep.covered_source_mass)
        self.assertEqual(rep.covered_source_mass, rep.source_mass)
        # 权重稀疏且非负，按目标汇总等于目标体积
        self.assertTrue(all(w > 0 for _, _, w in rep.weights))
        agg = [F(0)] * 5
        for ti, si, w in rep.weights:
            agg[ti] += w
        self.assertEqual(tuple(agg), rep.target_volumes)
        # 源网格对象只读：内容未变
        self.assertEqual(self.src.tets, tuple(tuple(t) for t in T6))

    def test_reverse_direction_consistent(self):
        rho = [F(1), F(2), F(4), F(8), F(16)]
        fwd = remap(self.tgt, rho, self.src)
        self.assertEqual(fwd.mass_error, 0)
        self.assertEqual(fwd.target_mass, fwd.source_mass)


class RemapPiecewiseConstantTest(unittest.TestCase):
    """跨 x=1/2 的分段常量场：源网格贴合分界面，目标不贴合。"""

    def setUp(self):
        lv, lt = box_tets(0, F(1, 2))
        rv, rt = box_tets(F(1, 2), 1)
        verts = lv + rv
        tets = [tuple(t) for t in lt] + [tuple(i + 8 for i in t) for t in rt]
        self.src = TetMesh.build(verts, tets)
        self.tgt = TetMesh.build(cube_vertices(), T6)
        # 左半密度 1，右半密度 2
        self.rho = [F(1)] * 6 + [F(2)] * 6

    def test_total_mass(self):
        rep = remap(self.src, self.rho, self.tgt)
        self.assertEqual(rep.mass_error, 0)
        # 总质量 = 1 * 1/2 + 2 * 1/2 = 3/2
        self.assertEqual(rep.source_mass, F(3, 2))
        self.assertEqual(rep.target_mass, F(3, 2))

    def test_crossing_tet_average_by_hand(self):
        # 目标单元 0 = (000,100,110,111)：被 x=1/2 截断。
        # 手算：x<=1/2 部分体积 1/48，x>=1/2 部分体积 7/48，
        # 质量 = 1/48*1 + 7/48*2 = 15/48 = 5/16，平均 = (5/16)/(1/6) = 15/8。
        rep = remap(self.src, self.rho, self.tgt)
        self.assertEqual(rep.target_volumes[0], F(1, 6))
        self.assertEqual(rep.target_masses[0], F(5, 16))
        self.assertEqual(rep.target_averages[0], F(15, 8))

    def test_fully_one_side_tets(self):
        # 目标单元 5 = (000,001,011,111) 与单元 0 关于中心对称：
        # x<=1/2 部分体积 7/48，x>=1/2 部分 1/48，
        # 质量 = 7/48*1 + 1/48*2 = 9/48 = 3/16，平均 = (3/16)/(1/6) = 9/8。
        rep = remap(self.src, self.rho, self.tgt)
        self.assertEqual(rep.target_masses[5], F(3, 16))
        self.assertEqual(rep.target_averages[5], F(9, 8))


class RemapPartialCoverageTest(unittest.TestCase):
    """目标区域大于源区域：保留缺口，不归一。"""

    def setUp(self):
        self.src = TetMesh.build(cube_vertices(), T6)
        # 目标：[0,2]x[0,1]x[0,1]，源只覆盖左半
        verts = [(F(x), F(y), F(z)) for x in (0, 2) for y in (0, 1) for z in (0, 1)]
        self.tgt = TetMesh.build(verts, T6)

    def test_gap_preserved(self):
        rep = remap(self.src, [F(3)] * 6, self.tgt)
        # 每个目标单元体积 1/3，只有一半被覆盖
        self.assertEqual(sum(rep.uncovered_volumes), F(1))
        self.assertEqual(sum(rep.covered_volumes), F(1))
        # 目标总质量 = 源被覆盖部分质量 = 3 * 1
        self.assertEqual(rep.target_mass, F(3))
        self.assertEqual(rep.covered_source_mass, F(3))
        self.assertEqual(rep.mass_error, 0)
        # 覆盖部分平均密度仍为 3（不归一、不稀释）
        for avg, cov in zip(rep.target_averages, rep.covered_volumes):
            if cov > 0:
                self.assertEqual(avg, F(3))

    def test_no_overlap_at_all(self):
        far = TetMesh.build(cube_vertices(shift=(F(5), F(0), F(0))), T6)
        rep = remap(self.src, [F(3)] * 6, far)
        self.assertEqual(rep.weights, ())
        self.assertEqual(rep.target_mass, 0)
        self.assertEqual(sum(rep.uncovered_volumes), 1)
        self.assertTrue(all(a is None for a in rep.target_averages))
        self.assertEqual(rep.mass_error, 0)


class RemapInputValidationTest(unittest.TestCase):
    def setUp(self):
        V = cube_vertices()
        self.mesh = TetMesh.build(V, T6)

    def test_float_density_rejected(self):
        with self.assertRaises(MeshInputError):
            remap(self.mesh, [0.5] * 6, self.mesh)

    def test_bool_density_rejected(self):
        with self.assertRaises(MeshInputError):
            remap(self.mesh, [True] * 6, self.mesh)

    def test_density_length_mismatch(self):
        with self.assertRaises(MeshInputError):
            remap(self.mesh, [F(1)] * 5, self.mesh)

    def test_zero_volume_tet_rejected(self):
        V = cube_vertices()
        bad = TetMesh.build(V, T6[:-1] + [(0, 1, 2, 3)])  # 末单元共面
        with self.assertRaises(MeshInputError):
            remap(bad, [F(1)] * 6, self.mesh)
        with self.assertRaises(MeshInputError):
            remap(self.mesh, [F(1)] * 6, bad)

    def test_non_mesh_rejected(self):
        with self.assertRaises(MeshInputError):
            remap("not a mesh", [F(1)] * 6, self.mesh)

    def test_report_is_frozen(self):
        rep = remap(self.mesh, [F(1)] * 6, self.mesh)
        self.assertIsInstance(rep, RemapReport)
        with self.assertRaises(Exception):
            rep.mass_error = F(1)  # type: ignore[misc]


if __name__ == "__main__":
    unittest.main()
