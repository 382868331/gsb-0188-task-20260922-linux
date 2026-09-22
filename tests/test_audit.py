"""必要测试：解析校验、定向修复、面关联/非流形、解析体积与面通量核对。"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    MAX_TETS,
    DuplicateTetrahedraError,
    MeshInputError,
    TetMesh,
    ZeroVolumeError,
    audit,
    repair_orientation,
)
from orient_audit.mesh import FACE_LOCAL

from fixtures import CUBE_INTERIOR_KEYS, CUBE_TETS, CUBE_VERTICES
from reference import ref_flux_volume, ref_outward, ref_signed_volume6


def _flip_last_two(tet):
    return (tet[0], tet[1], tet[3], tet[2])


class TestParsingAndValidation(unittest.TestCase):
    def test_rejects_float_coordinate(self):
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build([(0, 0, 0.1)], [])
        self.assertEqual(cm.exception.location, "verts[0][2]")

    def test_rejects_bool_coordinate_and_index(self):
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build([(0, 0, True)], [])
        self.assertIn("bool", str(cm.exception))
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build([(0, 0, 0), (1, 0, 0)], [(0, 1, 1, False)])
        self.assertEqual(cm.exception.location, "tets[0][3]")

    def test_rejects_nan_infinity_string(self):
        # 即便以字符串混入，也必须在解析处被拒绝并可定位。
        for bad in ("nan", "inf", "-Infinity", "1/0"):
            with self.subTest(bad=bad):
                with self.assertRaises(MeshInputError) as cm:
                    TetMesh.build([(0, 0, bad)], [])
                self.assertEqual(cm.exception.location, "verts[0][2]")

    def test_rejects_wrong_dimension(self):
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build([(0, 0)], [])
        self.assertEqual(cm.exception.location, "verts[0]")

    def test_rejects_bad_tet_size_and_duplicate_within_tet(self):
        v = CUBE_VERTICES
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build(v, [(0, 1, 2)])
        self.assertEqual(cm.exception.location, "tets[0]")
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build(v, [(0, 0, 1, 2)])
        self.assertEqual(cm.exception.location, "tets[0]")

    def test_rejects_index_out_of_range(self):
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build(CUBE_VERTICES, [(0, 1, 2, 99)])
        self.assertEqual(cm.exception.location, "tets[0][3]")

    def test_rejects_too_many_tets(self):
        # 构造 501 个互不相同的下标集合，确保不是被重复检查提前拦下。
        verts = [(i, 0, 0) for i in range(2004)]
        tets = [(4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3) for i in range(501)]
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build(verts, tets)
        self.assertEqual(cm.exception.location, "tets")

    def test_boundary_exactly_500_accepted(self):
        verts = [(i, 0, 0) for i in range(2000)]
        tets = [(4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3) for i in range(500)]
        mesh = TetMesh.build(verts, tets)
        self.assertEqual(mesh.tets[-1], (1996, 1997, 1998, 1999))

    def test_fraction_and_rational_string_coordinates(self):
        mesh = TetMesh.build(
            [("1/3", Fraction(2, 7), 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
            [(0, 1, 2, 3)],
        )
        self.assertEqual(mesh.vertices[0][0], Fraction(1, 3))
        self.assertEqual(mesh.vertices[0][1], Fraction(2, 7))


class TestDuplicateTets(unittest.TestCase):
    def test_duplicate_exact_copy_reported_separately(self):
        with self.assertRaises(DuplicateTetrahedraError) as cm:
            TetMesh.build(CUBE_VERTICES, (*CUBE_TETS, CUBE_TETS[0]))
        groups = cm.exception.duplicates
        self.assertEqual(len(groups), 1)
        key, positions = groups[0]
        self.assertEqual(key, tuple(sorted(CUBE_TETS[0])))
        self.assertEqual(positions, (0, 5))

    def test_duplicate_permutation_reported(self):
        perm = (CUBE_TETS[0][2], CUBE_TETS[0][0], CUBE_TETS[0][3], CUBE_TETS[0][1])
        with self.assertRaises(DuplicateTetrahedraError) as cm:
            TetMesh.build(CUBE_VERTICES, (CUBE_TETS[0], perm))
        self.assertEqual(cm.exception.duplicates[0][1], (0, 1))


class TestCubeAnalyticVolume(unittest.TestCase):
    """单立方体分割：解析体积核验 + 面通量独立核对。"""

    def setUp(self):
        self.mesh = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        self.report = audit(self.mesh)

    def test_each_tet_volume_matches_independent_reference(self):
        # 被测核心是 3x3 余子式；参考是 4x4 Bareiss，路径独立。
        expected = (
            Fraction(1, 6),
            Fraction(1, 6),
            Fraction(1, 6),
            Fraction(1, 6),
            Fraction(1, 3),
        )
        for ti, want in enumerate(expected):
            got = ref_signed_volume6(CUBE_VERTICES, CUBE_TETS[ti]) / 6
            self.assertEqual(got, want, f"参考体积 tets[{ti}]")
        self.assertEqual(
            self.report.total_signed_volume, Fraction(1), "总精确体积必须为 1"
        )

    def test_no_negative_or_zero(self):
        self.assertEqual(self.report.negative_tets, ())
        self.assertEqual(self.report.zero_tets, ())
        self.assertTrue(self.report.consistently_oriented)
        self.assertTrue(self.report.is_manifold)

    def test_interior_faces_exactly_four_and_opposite(self):
        self.assertEqual(len(self.report.interior_faces), 4)
        keys = sorted(f.key for f in self.report.interior_faces)
        self.assertEqual(keys, sorted(CUBE_INTERIOR_KEYS))
        self.assertTrue(all(f.opposite for f in self.report.interior_faces))

    def test_boundary_is_twelve_triangles_and_flux_volume_one(self):
        self.assertEqual(len(self.report.boundary_faces), 12)
        triples = [f.ordered for f in self.report.boundary_faces]
        # 面通量公式独立给出总体积 1（不经过单元体积）。
        self.assertEqual(ref_flux_volume(CUBE_VERTICES, triples), Fraction(1))

    def test_boundary_faces_are_outward_independent_check(self):
        for bf in self.report.boundary_faces:
            self.assertTrue(
                ref_outward(CUBE_VERTICES, CUBE_TETS[bf.tet_index], bf.ordered),
                f"边界面 {bf.ordered} 应朝外",
            )

    def test_single_connected_component(self):
        self.assertEqual(self.report.components, (tuple(range(5)),))


class TestRepairOrientation(unittest.TestCase):
    def test_flip_negative_tet_only_swap_last_two(self):
        flipped = [list(t) for t in CUBE_TETS]
        flipped[1] = list(_flip_last_two(CUBE_TETS[1]))  # 单元 1 反向
        flipped[3] = list(_flip_last_two(CUBE_TETS[3]))  # 单元 3 反向
        mesh = TetMesh.build(CUBE_VERTICES, flipped)
        report = audit(mesh)
        self.assertEqual(report.negative_tets, (1, 3))

        fixed = repair_orientation(mesh)
        # 原网格不被修改。
        self.assertEqual(mesh.tets[1], tuple(flipped[1]))
        # 修复恰为交换位 2、3；其余单元完全不动。
        self.assertEqual(fixed.tets[1], CUBE_TETS[1])
        self.assertEqual(fixed.tets[3], CUBE_TETS[3])
        for ti in (0, 2, 4):
            self.assertEqual(fixed.tets[ti], CUBE_TETS[ti])
        # 顶点对象恒等：修复不移动顶点。
        self.assertTrue(
            all(a is b for a, b in zip(fixed.vertices, mesh.vertices, strict=True))
        )

        fixed_report = audit(fixed)
        self.assertEqual(fixed_report.negative_tets, ())
        self.assertEqual(fixed_report.total_signed_volume, Fraction(1))
        # 修复后所有内面方向恢复一致。
        self.assertTrue(
            all(f.opposite for f in fixed_report.interior_faces)
        )
        self.assertEqual(ref_flux_volume(
            CUBE_VERTICES, [f.ordered for f in fixed_report.boundary_faces]
        ), Fraction(1))

    def test_repair_is_idempotent(self):
        mesh = TetMesh.build(
            CUBE_VERTICES, [_flip_last_two(t) for t in CUBE_TETS]
        )
        once = repair_orientation(mesh)
        twice = repair_orientation(once)
        self.assertEqual(once.tets, twice.tets)

    def test_repair_after_even_permutation_keeps_positive(self):
        # 偶置换（3-循环）不改变符号，修复不应触碰该单元。
        a, b, c, d = CUBE_TETS[0]
        even = (b, c, a, d)  # 前三位循环，偶置换
        mesh = TetMesh.build(CUBE_VERTICES, [even, *CUBE_TETS[1:]])
        report = audit(mesh)
        self.assertNotIn(0, report.negative_tets)
        fixed = repair_orientation(mesh)
        self.assertEqual(fixed.tets[0], even)


class TestZeroVolume(unittest.TestCase):
    def test_coplanar_tet_detected_and_repair_refused(self):
        # 4 个共面顶点（z 全为 0）。
        verts = (*CUBE_VERTICES, (1, 1, 0))  # 下标 8 与 2 共面于 z=0
        tets = [*CUBE_TETS, (0, 1, 2, 8)]
        mesh = TetMesh.build(verts, tets)
        report = audit(mesh)
        self.assertEqual(report.zero_tets, (5,))
        # 零体积不影响其他单元的定向判定。
        self.assertEqual(report.negative_tets, ())

        with self.assertRaises(ZeroVolumeError) as cm:
            repair_orientation(mesh)
        self.assertEqual(cm.exception.tet_indices, (5,))
        # 拒绝修复时输入网格保持原状（无部分变更）。
        self.assertEqual(mesh.tets[-1], (0, 1, 2, 8))

    def test_audit_accepts_zero_volume_and_reports_faces(self):
        # audit 不抛异常，仍完成面关联与体积统计。
        verts = (*CUBE_VERTICES, (1, 1, 0))
        mesh = TetMesh.build(verts, [*CUBE_TETS, (0, 1, 2, 8)])
        report = audit(mesh)
        self.assertTrue(report.has_zero_volume)
        # 退化总体积仍精确为 1。
        self.assertEqual(report.total_signed_volume, Fraction(1))


class TestNonManifold(unittest.TestCase):
    def test_three_tets_sharing_one_face_is_non_manifold(self):
        # 共享三点面 (0,2,5)：立方体中心 tets[4] 被再加一份不同的第四点
        # 版本，需要避免重复单元。构造三个单元共用面 (0,2,5)，第四点
        # 分别取 1、7、3（均与面不共面且四点集合不同）。
        verts = (
            (0, 0, 0),  # 0
            (1, 1, 0),  # 2 -> 1
            (1, 0, 1),  # 5 -> 2  面三点
            (1, 0, 0),  # 1 -> 3
            (0, 1, 1),  # 7 -> 4
            (0, 1, 0),  # 3 -> 5
        )
        tets = [
            (0, 1, 2, 3),  # 第四点 3
            (0, 1, 2, 4),  # 第四点 4
            (0, 1, 2, 5),  # 第四点 5
        ]
        mesh = TetMesh.build(verts, tets)
        report = audit(mesh)
        self.assertEqual(len(report.non_manifold_faces), 1)
        nm = report.non_manifold_faces[0]
        self.assertEqual(nm.key, (0, 1, 2))
        self.assertEqual(nm.tet_indices, (0, 1, 2))
        self.assertFalse(report.is_manifold)
        # 非流形不阻止定向修复：逐单元翻转照常工作。
        fixed = repair_orientation(mesh)
        self.assertEqual(audit(fixed).negative_tets, ())

    def test_two_tets_same_induced_orientation_is_inconsistent(self):
        # 两单元共享面 (0,1,2)，第四点 3、4 分居面两侧。
        # 两个单元都按“同一种顶点书写序”给出时，处在负侧的单元
        # 体积为负，且对共享面诱导出与正侧单元*相同*的朝向。
        verts = (
            (0, 0, 0),
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (0, 0, -1),
        )
        t1 = (0, 1, 2, 3)  # 正侧，正体积，面诱导序 (2,1,0)
        t2 = (0, 1, 2, 4)  # 负侧，负体积，面诱导序同样为 (2,1,0)
        mesh = TetMesh.build(verts, [t1, t2])
        report = audit(mesh)
        self.assertEqual(report.negative_tets, (1,))
        self.assertEqual(len(report.interior_faces), 1)
        self.assertFalse(report.interior_faces[0].opposite)
        self.assertEqual(report.inconsistent_faces[0].key, (0, 1, 2))
        # 修复翻转后内面恢复相反朝向。
        fixed = repair_orientation(mesh)
        self.assertTrue(all(f.opposite for f in audit(fixed).interior_faces))


class TestComponents(unittest.TestCase):
    def test_disjoint_tet_sets_form_two_components(self):
        # 两个互不接触的四面体。
        verts = [
            (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (5, 5, 5), (6, 5, 5), (5, 6, 5), (5, 5, 6),
        ]
        mesh = TetMesh.build(verts, [(0, 1, 2, 3), (4, 5, 6, 7)])
        report = audit(mesh)
        self.assertEqual(report.components, ((0,), (1,)))


class TestFaceWindingTable(unittest.TestCase):
    def test_face_local_are_all_odd_permutations_of_complement(self):
        # 契约固定：四个诱导面均为外法向朝向（相对 (i,j,k,l) 奇置换）。
        for local in FACE_LOCAL:
            self.assertEqual(len(set(local)), 3)
        self.assertEqual(tuple(sorted(sum(FACE_LOCAL, ()))), (0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3))


if __name__ == "__main__":
    unittest.main()
