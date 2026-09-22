"""核心验收：解析体积、独立通量核对、定向/边界/非流形审计与修复。"""

import unittest
from fractions import Fraction

from tetmesh_audit import audit, repair_orientation, samples
from tetmesh_audit import DegenerateTetrahedronError
from tetmesh_audit.geometry import FACE_REMOVAL_ORDER, outward_triangle

try:  # 兼容 `unittest discover -s tests`（顶层导入）与包内导入
    from tests import reference
except ImportError:  # pragma: no cover
    import reference


def _verts_as_frac(verts):
    return [tuple(Fraction(x) for x in v) for v in verts]


class CubeAnalyticTests(unittest.TestCase):
    """单立方体 6 四面体分割：解析体积 + 独立面通量核对。"""

    def setUp(self):
        self.mesh = samples.cube_tetmesh()
        self.result = audit(self.mesh)
        self.frac_verts = _verts_as_frac(samples.CUBE_VERTICES)

    def test_each_tet_six_volume_independently(self):
        # 参考实现（Leibniz 展开）逐单元核对，不调用被测行列式。
        for i, t in enumerate(samples.CUBE_TETS):
            pts = [self.frac_verts[j] for j in t]
            expected = reference.tet_signed_six_volume(*pts)
            self.assertEqual(self.result.six_volumes[i], expected)
            self.assertEqual(expected, 1, f"tet {i} six-volume should be +1")

    def test_total_volume_analytic(self):
        self.assertEqual(self.result.total_signed_six_volume, 6)
        self.assertEqual(self.result.total_signed_volume, Fraction(1, 1))
        self.assertEqual(
            self.result.total_signed_volume, reference.box_volume_analytic(1)
        )

    def test_boundary_triangle_count_and_no_duplicate_faces(self):
        # 立方体 6 个正方形面，每个由 2 个三角形覆盖 → 12 个边界三角面。
        self.assertEqual(len(self.result.boundary_triangles), 12)
        keys = [tuple(sorted(bt.vertices)) for bt in self.result.boundary_triangles]
        self.assertEqual(len(keys), len(set(keys)))
        # 扇形剖分：6 个径向内面（每对相邻四面体一个），无重复键；
        # 总面关联数 6*4=24 = 12 边界 + 6*2 内面。
        self.assertEqual(len(self.result.interior_faces), 6)

    def test_flux_formula_independent_check(self):
        # 独立参考：1/6 Σ det(p,q,r)，只读取审计输出的定向三角面，
        # 行列式与求和均在 reference 中独立完成。
        tris = [bt.vertices for bt in self.result.boundary_triangles]
        got = reference.flux_volume(self.frac_verts, tris)
        self.assertEqual(got, Fraction(1, 1))
        self.assertEqual(got, self.result.total_signed_volume)
        self.assertEqual(self.result.boundary_det_sum, 6)

    def test_boundary_outward_normals(self):
        # 每个边界三角面的外法向必须指向立方体外（按所在面抽查坐标分量）。
        for bt in self.result.boundary_triangles:
            p = self.frac_verts[bt.vertices[0]]
            q = self.frac_verts[bt.vertices[1]]
            r = self.frac_verts[bt.vertices[2]]
            u = (q[0]-p[0], q[1]-p[1], q[2]-p[2])
            v = (r[0]-p[0], r[1]-p[1], r[2]-p[2])
            n = (
                u[1]*v[2]-u[2]*v[1],
                u[2]*v[0]-u[0]*v[2],
                u[0]*v[1]-u[1]*v[0],
            )
            cx = (p[0]+q[0]+r[0]) / 3
            cy = (p[1]+q[1]+r[1]) / 3
            cz = (p[2]+q[2]+r[2]) / 3
            # 形心落在哪个坐标面，法向对应分量必须朝外。
            if cx == 0:
                self.assertLess(n[0], 0)
            elif cx == 1:
                self.assertGreater(n[0], 0)
            if cy == 0:
                self.assertLess(n[1], 0)
            elif cy == 1:
                self.assertGreater(n[1], 0)
            if cz == 0:
                self.assertLess(n[2], 0)
            elif cz == 1:
                self.assertGreater(n[2], 0)

    def test_interior_faces_induced_opposite_and_reference_incidence(self):
        ref_inc = reference.face_incidence(samples.CUBE_TETS)
        manifold = {k: v for k, v in ref_inc.items() if len(v) == 2}
        self.assertEqual(len(manifold), 6)
        for f in self.result.interior_faces:
            self.assertTrue(
                f.consistent,
                msg=f"interior face {f.vertices} not consistently oriented",
            )
            # 独立判据：相对排序面键的排列奇偶性必须不同。
            key = tuple(sorted(f.vertices))

            def parity(tri):
                ranks = [key.index(v) for v in tri]
                inv = sum(
                    1 for a in range(3) for b in range(a + 1, 3)
                    if ranks[a] > ranks[b]
                )
                return inv % 2

            tr0, tr1 = f.induced
            self.assertNotEqual(parity(tr0), parity(tr1))
            # 两面诱导方向相反（外法向互指），三元组字面值自然不同。
            self.assertNotEqual(tr0, tr1)
        self.assertTrue(self.result.is_consistent)
        self.assertEqual(self.result.defects, ())
        self.assertEqual(len(self.result.components), 1)

    def test_negative_volume_equivalence_to_odd_permutation(self):
        # 独立枚举：对固定四点的全部 24 种排列，参考行列式符号
        # 必须与审计为每个排列给出的符号一致。
        from itertools import permutations
        from tetmesh_audit import TetMesh

        base = samples.CUBE_TETS[0]
        pts = [self.frac_verts[j] for j in base]
        for perm in permutations(base):
            mesh = TetMesh(samples.CUBE_VERTICES, [perm])
            res = audit(mesh)
            ref_sign = reference.tet_signed_six_volume(
                *[self.frac_verts[j] for j in perm]
            )
            self.assertEqual(
                (res.six_volumes[0] > 0) - (res.six_volumes[0] < 0),
                (ref_sign > 0) - (ref_sign < 0),
            )


class OrientationRepairTests(unittest.TestCase):
    def test_flipped_tet_detected_and_repaired_by_last_two_swap(self):
        mesh = samples.cube_with_flipped_tet(2)
        res = audit(mesh)
        self.assertEqual(res.negative_tets, (2,))
        # 修复前：受影响的两个内面方向不再相反。
        bad = [d for d in res.defects if d.kind == "inconsistent_interior_face"]
        self.assertTrue(bad)

        repaired = repair_orientation(mesh)
        self.assertEqual(repaired.flipped_tets, (2,))
        # 精确验证修复动作：仅交换单元 2 的最后两个顶点，其余不变。
        old = mesh.tetrahedra[2]
        new = repaired.repaired_mesh.tetrahedra[2]
        self.assertEqual(new, (old[0], old[1], old[3], old[2]))
        for i in range(6):
            if i != 2:
                self.assertEqual(repaired.repaired_mesh.tetrahedra[i], mesh.tetrahedra[i])
        # 顶点未移动。
        self.assertEqual(repaired.repaired_mesh.vertices, mesh.vertices)

        after = repaired.audit_after
        self.assertEqual(after.negative_tets, ())
        self.assertEqual(after.zero_volume_tets, ())
        self.assertTrue(after.is_consistent)
        self.assertEqual(after.total_signed_volume, Fraction(1, 1))
        self.assertEqual(after.defects, ())
        # 修复后面通量恒等式恢复。
        tris = [bt.vertices for bt in after.boundary_triangles]
        self.assertEqual(
            reference.flux_volume(_verts_as_frac(samples.CUBE_VERTICES), tris),
            Fraction(1, 1),
        )

    def test_repair_is_functional_original_unchanged(self):
        mesh = samples.cube_with_flipped_tet(0)
        before_tets = mesh.tetrahedra
        repair_orientation(mesh)
        self.assertEqual(mesh.tetrahedra, before_tets)

    def test_all_tets_flipped_repair_restores_volume(self):
        tets = [(t[0], t[1], t[3], t[2]) for t in samples.CUBE_TETS]
        from tetmesh_audit import TetMesh

        mesh = TetMesh(samples.CUBE_VERTICES, tets)
        res = audit(mesh)
        self.assertEqual(len(res.negative_tets), 6)
        self.assertEqual(res.total_signed_volume, Fraction(-1, 1))
        repaired = repair_orientation(mesh)
        self.assertEqual(len(repaired.flipped_tets), 6)
        self.assertTrue(repaired.audit_after.is_consistent)
        self.assertEqual(repaired.audit_after.total_signed_volume, Fraction(1, 1))


class NonManifoldTests(unittest.TestCase):
    def test_three_incident_tets_flagged_with_location(self):
        mesh = samples.nonmanifold_mesh()
        res = audit(mesh)
        nm = [d for d in res.defects if d.kind == "non_manifold_face"]
        self.assertEqual(len(nm), 1)
        self.assertEqual(nm[0].face, (0, 1, 2))
        self.assertEqual(len(nm[0].tet_indices), 3)
        # 参考独立统计同样得到该面关联 3 个单元。
        inc = reference.face_incidence(mesh.tetrahedra)
        self.assertEqual(len(inc[(0, 1, 2)]), 3)
        self.assertFalse(res.is_consistent)

    def test_repair_fixes_signs_but_reports_nonmanifold(self):
        mesh = samples.nonmanifold_mesh()
        repaired = repair_orientation(mesh)  # 非流形不阻止符号修复
        after = repaired.audit_after
        self.assertEqual(after.negative_tets, ())
        nm = [d for d in after.defects if d.kind == "non_manifold_face"]
        self.assertEqual(len(nm), 1)
        self.assertEqual(nm[0].face, (0, 1, 2))
        self.assertFalse(after.is_consistent)

    def test_same_side_overlap_is_inconsistent_interface_not_nonmanifold(self):
        mesh = samples.overlapping_same_side_mesh()
        res = audit(mesh)
        kinds = [d.kind for d in res.defects]
        # 共享面恰 2 个单元但两第四顶点在面的同一侧：诱导方向相同。
        self.assertEqual(kinds, ["inconsistent_interior_face"])
        # 两个单元各自都是正体积，故没有负向单元可翻转来消除此缺陷。
        self.assertEqual(res.negative_tets, ())
        self.assertFalse(
            any(d.kind == "non_manifold_face" for d in res.defects)
        )
        # 修复符号不能消除该结构缺陷。
        repaired = repair_orientation(mesh)
        self.assertEqual(repaired.flipped_tets, ())
        self.assertTrue(
            any(d.kind == "inconsistent_interior_face"
                for d in repaired.audit_after.defects)
        )


class DegenerateTests(unittest.TestCase):
    def test_coplanar_tet_zero_volume_exact(self):
        mesh = samples.coplanar_degenerate_mesh()
        res = audit(mesh)
        self.assertIn(6, res.zero_volume_tets)
        # 参考独立确认行列式为零。
        t = mesh.tetrahedra[6]
        pts = [_verts_as_frac(mesh.vertices)[j] for j in t]
        self.assertEqual(reference.tet_signed_six_volume(*pts), 0)
        # 立方体其余部分体积仍精确为 1。
        self.assertEqual(res.total_signed_volume, Fraction(1, 1))
        defects = [d for d in res.defects if d.kind == "zero_volume_tet"]
        self.assertEqual(defects[0].tet_indices, (6,))
        self.assertFalse(res.is_consistent)

    def test_repair_refused_atomically(self):
        mesh = samples.coplanar_degenerate_mesh()
        original = mesh.tetrahedra
        with self.assertRaises(DegenerateTetrahedronError):
            repair_orientation(mesh)
        # 原子性：失败后原网格状态不变，仍可审计。
        res = audit(mesh)
        self.assertIn(6, res.zero_volume_tets)
        self.assertEqual(mesh.tetrahedra, original)

    def test_nearly_coplanar_but_exact_rational_is_not_zero(self):
        # 看似“接近共面”但有理坐标精确非零：不得被容差误判为零。
        from tetmesh_audit import TetMesh

        verts = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, Fraction(1, 10**9))]
        mesh = TetMesh(verts, [(0, 1, 2, 3)])
        res = audit(mesh)
        self.assertEqual(res.six_volumes[0], Fraction(1, 10**9))
        self.assertEqual(res.zero_volume_tets, ())


class ComponentTests(unittest.TestCase):
    def test_two_disjoint_cubes_two_components(self):
        mesh = samples.two_disjoint_cubes()
        res = audit(mesh)
        self.assertEqual(len(res.components), 2)
        self.assertEqual([len(c.tet_indices) for c in res.components], [6, 6])
        self.assertEqual(res.total_signed_volume, Fraction(2, 1))
        # 各分量体积独立为 1。
        self.assertEqual(res.components[0].signed_six_volume, 6)
        self.assertEqual(res.components[1].signed_six_volume, 6)
        # 24 个边界三角面（两个立方体互不为内面）。
        self.assertEqual(len(res.boundary_triangles), 24)

    def test_scaled_fractional_cube_exact(self):
        mesh = samples.scaled_cube_tetmesh(Fraction(2, 3))
        res = audit(mesh)
        self.assertEqual(res.total_signed_volume, Fraction(8, 27))
        self.assertEqual(
            res.total_signed_volume, reference.box_volume_analytic(Fraction(2, 3))
        )
        tris = [bt.vertices for bt in res.boundary_triangles]
        self.assertEqual(
            reference.flux_volume(_verts_as_frac(mesh.vertices), tris),
            Fraction(8, 27),
        )


class SingleTetBoundaryTests(unittest.TestCase):
    def test_single_tet_four_boundary_faces_flux(self):
        from tetmesh_audit import TetMesh

        verts = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        mesh = TetMesh(verts, [(0, 1, 2, 3)])
        res = audit(mesh)
        self.assertEqual(len(res.boundary_triangles), 4)
        self.assertEqual(res.total_signed_volume, Fraction(1, 6))
        tris = [bt.vertices for bt in res.boundary_triangles]
        self.assertEqual(
            reference.flux_volume(_verts_as_frac(verts), tris), Fraction(1, 6)
        )

    def test_outward_triangle_rule_matches_boundary_operator(self):
        # 对正体积单元，外法向三角面 = (-1)^i * 删除序；逐面独立核对。
        t = (0, 1, 2, 3)
        vol = 1
        for lf, rem in enumerate(FACE_REMOVAL_ORDER):
            got = outward_triangle(t, lf, vol)
            base = (t[rem[0]], t[rem[1]], t[rem[2]])
            expect = base if lf % 2 == 0 else (base[2], base[1], base[0])
            self.assertEqual(got, expect)


class FluxIdentityTests(unittest.TestCase):
    """边界通量恒等式仅在闭流形且定向一致时成立（精确区分边界错误）。"""

    def test_identity_breaks_with_flipped_tet(self):
        mesh = samples.cube_with_flipped_tet(2)
        res = audit(mesh)
        # 几何外法向通量给出各单元绝对体积之和，与有向总体积不再相等。
        self.assertEqual(res.boundary_det_sum, 6)
        self.assertEqual(res.total_signed_six_volume, 4)
        self.assertNotEqual(res.boundary_det_sum, res.total_signed_six_volume)

    def test_identity_breaks_with_nonmanifold_face(self):
        mesh = samples.nonmanifold_mesh()
        res = audit(mesh)
        # 三个单元共用一面，通量无法成对抵消。
        self.assertNotEqual(res.boundary_det_sum, res.total_signed_six_volume)

    def test_identity_holds_for_two_components(self):
        mesh = samples.two_disjoint_cubes()
        res = audit(mesh)
        self.assertEqual(res.boundary_det_sum, res.total_signed_six_volume)
        self.assertEqual(res.boundary_det_sum, 12)


if __name__ == "__main__":
    unittest.main()
