"""内部边分裂 split_interior_edge 的必要测试。

期望值来源（不调用被测实现生成）：

- 插点坐标 P = A + t(B-A)、子体积 (1-t)·V 与 t·V 由测试按解析关系手算，
  体积再用 tests/reference.py 的 4x4 Bareiss 独立实现复核
  （被测核心走 3x3 余子式，是另一条算法路径）；
- 外四面体/重心剖分的总体积 32/3 为解析值，并用独立面通量公式核对
  分裂前后边界面；
- 非流形边、夹伤边、零体积星、500 上限等夹具均为小规模手构造枚举，
  原因码在 fixtures.py 注释中给出拓扑依据。
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    MAX_TETS,
    EdgeSplitError,
    MeshInputError,
    ProtectedFaceError,
    TetMesh,
    audit,
    repair_orientation,
    split_interior_edge,
)

from fixtures import (
    BARY_INTERIOR_EDGE,
    BARY_LABELS,
    BARY_STAR_TETS,
    BARY_TETS,
    BARY_VERTICES,
    DEGEN_STAR_TETS,
    DEGEN_STAR_VERTICES,
    PINCH_TETS,
    PINCH_VERTICES,
    PINWHEEL_TETS,
    PINWHEEL_VERTICES,
    ring_star_vertices_tets,
)
from reference import ref_flux_volume, ref_signed_volume6


# 外四面体（0..3）6 倍体积：底面积 8、高 4，V = 32/3，V6 = 64。
OUTER_V6 = Fraction(64)
# 重心剖分每个子单元 6 倍体积 = 64/4 = 16。
BARY_CHILD_V6 = Fraction(16)
T_PARAM = Fraction(1, 3)


def _bary_mesh(labels=None, protected_faces=None):
    return TetMesh.build(
        BARY_VERTICES,
        BARY_TETS,
        tet_labels=labels,
        protected_faces=protected_faces,
    )


def _ref_total_v6(vertices, tets):
    """用独立 Bareiss 参考对一组单元求 6 倍有向体积和。"""

    return sum((ref_signed_volume6(vertices, tet) for tet in tets), Fraction(0))


class TestSplitHappyPath(unittest.TestCase):
    def setUp(self):
        self.mesh = _bary_mesh()
        self.before = audit(self.mesh)

    def test_baseline_barycentric_mesh(self):
        # 4 个等体积子单元，外边界恰 4 面，6 个内面全部反向，单连通分量。
        self.assertEqual(self.before.total_signed_volume6, OUTER_V6)
        self.assertEqual(len(self.before.boundary_faces), 4)
        self.assertEqual(len(self.before.interior_faces), 6)
        self.assertTrue(self.before.consistently_oriented)
        self.assertTrue(self.before.is_manifold)

    def test_split_geometry_and_counts(self):
        result = split_interior_edge(self.mesh, *BARY_INTERIOR_EDGE, T_PARAM)
        new_mesh = result.mesh
        after = audit(new_mesh)

        # 元数据。
        self.assertEqual(result.edge, BARY_INTERIOR_EDGE)
        self.assertEqual(result.new_vertex_index, len(BARY_VERTICES))
        self.assertEqual(result.parameter, T_PARAM)
        self.assertEqual(result.star_tet_indices, BARY_STAR_TETS)
        self.assertEqual(len(new_mesh.vertices), len(BARY_VERTICES) + 1)
        self.assertEqual(after.n_tets, len(BARY_TETS) + len(BARY_STAR_TETS))

        # 插点坐标：A + t(B-A)，a=4=(1,1,1)，b=0=(0,0,0)，t=1/3。
        a, b = BARY_INTERIOR_EDGE
        want_p = tuple(
            Fraction(BARY_VERTICES[a][k])
            + T_PARAM * (Fraction(BARY_VERTICES[b][k]) - Fraction(BARY_VERTICES[a][k]))
            for k in range(3)
        )
        self.assertEqual(new_mesh.vertices[result.new_vertex_index], want_p)
        self.assertEqual(want_p, (Fraction(2, 3),) * 3)

        # 审计：无负向/零体积，仍为流形且定向一致。
        self.assertEqual(after.negative_tets, ())
        self.assertEqual(after.zero_tets, ())
        self.assertTrue(after.is_manifold)
        self.assertTrue(after.consistently_oriented)

        # 内面数：6 - 3（3 个含分裂边的内面消失）+ 6（3 对新侧面）
        #       + 3（切面）= 12。
        self.assertEqual(len(after.interior_faces), 12)
        # 边界面仍恰为外四面体 4 个面。
        self.assertEqual(len(after.boundary_faces), 4)
        outer_faces = {(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)}
        self.assertEqual(
            {tuple(sorted(f.ordered)) for f in after.boundary_faces},
            outer_faces,
        )
        # 单一连通分量（7 个单元）。
        self.assertEqual(after.components, (tuple(range(7)),))

    def test_volume_conservation_and_split_ratios_independent(self):
        result = split_interior_edge(self.mesh, *BARY_INTERIOR_EDGE, T_PARAM)
        nm = result.mesh

        # 被测核心总 6 倍体积精确守恒。
        self.assertEqual(audit(nm).total_signed_volume6, OUTER_V6)

        # 独立 Bareiss 路径：对“旧顶点 + 新点”的完整坐标表重算全部新单元。
        all_vertices = (*BARY_VERTICES, nm.vertices[result.new_vertex_index])
        self.assertEqual(
            _ref_total_v6(all_vertices, nm.tets), OUTER_V6
        )

        # 逐对父子单元：T_a = (1-t)·V6，T_b = t·V6（独立参考）。
        p_index = result.new_vertex_index
        for parent, (ia, ib) in result.child_tet_pairs:
            parent_v6 = ref_signed_volume6(BARY_VERTICES, BARY_TETS[parent])
            self.assertEqual(parent_v6, BARY_CHILD_V6)
            va = ref_signed_volume6(all_vertices, nm.tets[ia])
            vb = ref_signed_volume6(all_vertices, nm.tets[ib])
            self.assertEqual(va, (1 - T_PARAM) * parent_v6)
            self.assertEqual(vb, T_PARAM * parent_v6)
            # 子四点集合 = 父集合把 a（或 b）换成新点。
            a, b = BARY_INTERIOR_EDGE
            self.assertEqual(set(nm.tets[ia]), set(BARY_TETS[parent]) - {a} | {p_index})
            self.assertEqual(set(nm.tets[ib]), set(BARY_TETS[parent]) - {b} | {p_index})

        # 未参与分裂的单元 0 原样保留。
        self.assertEqual(result.untouched_tet_indices, (0,))
        self.assertEqual(nm.tets[0], BARY_TETS[0])

    def test_boundary_flux_unchanged_independent(self):
        # 边界面几何逐点不变 → 独立面通量公式在分裂前后给出同一解析体积。
        before_triples = [f.ordered for f in self.before.boundary_faces]
        self.assertEqual(
            ref_flux_volume(BARY_VERTICES, before_triples), Fraction(32, 3)
        )
        result = split_interior_edge(self.mesh, *BARY_INTERIOR_EDGE, T_PARAM)
        after = audit(result.mesh)
        after_triples = [f.ordered for f in after.boundary_faces]
        self.assertEqual(
            ref_flux_volume(result.mesh.vertices, after_triples), Fraction(32, 3)
        )
        # 边界面的三点键集合逐点相同。
        self.assertEqual(
            result.mesh.boundary_face_keys(), self.mesh.boundary_face_keys()
        )

    def test_material_labels_inherited(self):
        mesh = _bary_mesh(labels=BARY_LABELS)
        result = split_interior_edge(mesh, *BARY_INTERIOR_EDGE, T_PARAM)
        labels = result.mesh.tet_labels
        # 未分裂单元 0 标签不变。
        self.assertEqual(labels[0], "copper")
        # 单元 1、2、3 各产生两个同标签子单元。
        for parent, (ia, ib) in result.child_tet_pairs:
            self.assertEqual(labels[ia], BARY_LABELS[parent])
            self.assertEqual(labels[ib], BARY_LABELS[parent])
        self.assertEqual(
            labels,
            ("copper", "rubber", "rubber", "rubber", "rubber", "glass", "glass"),
        )

    def test_default_parameter_is_midpoint(self):
        result = split_interior_edge(self.mesh, *BARY_INTERIOR_EDGE)
        self.assertEqual(result.parameter, Fraction(1, 2))
        self.assertEqual(
            result.mesh.vertices[result.new_vertex_index], (Fraction(1, 2),) * 3
        )
        # 中点：两个子单元体积相等 = 父的一半（V6 = 8）。
        all_vertices = (
            *BARY_VERTICES,
            result.mesh.vertices[result.new_vertex_index],
        )
        for _parent, (ia, ib) in result.child_tet_pairs:
            self.assertEqual(
                ref_signed_volume6(all_vertices, result.mesh.tets[ia]), Fraction(8)
            )
            self.assertEqual(
                ref_signed_volume6(all_vertices, result.mesh.tets[ib]), Fraction(8)
            )

    def test_rational_string_parameter_accepted(self):
        result = split_interior_edge(
            self.mesh, *BARY_INTERIOR_EDGE, "1/4"
        )
        self.assertEqual(result.parameter, Fraction(1, 4))

    def test_reversed_edge_with_complementary_parameter_same_point(self):
        a, b = BARY_INTERIOR_EDGE  # (4, 0)
        r1 = split_interior_edge(self.mesh, a, b, Fraction(1, 3))
        r2 = split_interior_edge(self.mesh, b, a, Fraction(2, 3))
        self.assertEqual(
            r1.mesh.vertices[r1.new_vertex_index],
            r2.mesh.vertices[r2.new_vertex_index],
        )

    def test_star_size_four_on_ring(self):
        verts, tets = ring_star_vertices_tets(4)
        mesh = TetMesh.build(verts, tets)
        result = split_interior_edge(mesh, 0, 1, Fraction(2, 5))
        after = audit(result.mesh)
        self.assertEqual(after.n_tets, 8)
        self.assertTrue(after.is_manifold)
        self.assertEqual(after.total_signed_volume6, audit(mesh).total_signed_volume6)
        self.assertEqual(len(result.star_tet_indices), 4)

    def test_existing_vertices_identity_preserved(self):
        result = split_interior_edge(self.mesh, *BARY_INTERIOR_EDGE)
        self.assertTrue(
            all(x is y for x, y in zip(result.mesh.vertices, self.mesh.vertices))
        )

    def test_split_preserves_orientation_sign_for_negative_tet(self):
        # 故意把一个星内单元反向：分裂后两个子单元也都为负，
        # 且*有向*体积（含符号）依然精确守恒。
        flipped = [list(t) for t in BARY_TETS]
        flipped[1] = [flipped[1][0], flipped[1][1], flipped[1][3], flipped[1][2]]
        mesh = TetMesh.build(BARY_VERTICES, flipped)
        signed_before = audit(mesh).total_signed_volume6
        result = split_interior_edge(mesh, *BARY_INTERIOR_EDGE, Fraction(1, 3))
        after = audit(result.mesh)
        self.assertEqual(after.total_signed_volume6, signed_before)
        # 单元 1 的两个子单元均为负；其余全正。
        neg = set(after.negative_tets)
        ia, ib = dict(result.child_tet_pairs)[1]
        self.assertIn(ia, neg)
        self.assertIn(ib, neg)
        self.assertEqual(len(neg), 2)
        # 与 repair_orientation 协作：修复后体积恢复 32/3。
        fixed = repair_orientation(result.mesh)
        self.assertEqual(audit(fixed).total_signed_volume, Fraction(32, 3))


class TestSplitBoundaryProtection(unittest.TestCase):
    def test_boundary_edge_rejected_on_cube(self):
        from fixtures import CUBE_TETS, CUBE_VERTICES

        mesh = TetMesh.build(CUBE_VERTICES, CUBE_TETS)
        original_tets = mesh.tets
        # 立方体底边 (0,1) 同时属于两个边界三角面。
        with self.assertRaises(ProtectedFaceError) as cm:
            split_interior_edge(mesh, 0, 1, Fraction(1, 2))
        err = cm.exception
        self.assertEqual(err.reason, "protected_face")
        self.assertEqual(err.edge, (0, 1))
        self.assertEqual(len(err.boundary_faces), 2)
        self.assertEqual(err.explicit_faces, ())
        # 原子：原网格不被修改。
        self.assertIs(mesh.tets, original_tets)
        self.assertEqual(mesh.tets, tuple(CUBE_TETS))
        self.assertEqual(len(mesh.vertices), 8)

    def test_boundary_edge_on_ring_rejected(self):
        # ring3 中链环顶点之间的边 (2,3) 是边界面上的边。
        verts, tets = ring_star_vertices_tets(3)
        mesh = TetMesh.build(verts, tets)
        with self.assertRaises(ProtectedFaceError) as cm:
            split_interior_edge(mesh, 2, 3)
        self.assertEqual(len(cm.exception.boundary_faces), 2)

    def test_explicit_protected_interior_face_rejected(self):
        # 内面 (0,3,4) 含内部边 (4,0)；声明保护后必须拒绝，
        # 而同一操作在无保护时成功。
        protected = [(0, 3, 4)]
        mesh_prot = _bary_mesh(protected_faces=protected)
        with self.assertRaises(ProtectedFaceError) as cm:
            split_interior_edge(mesh_prot, 4, 0)
        self.assertEqual(cm.exception.explicit_faces, ((0, 3, 4),))
        self.assertEqual(cm.exception.boundary_faces, ())

        mesh_plain = _bary_mesh()
        split_interior_edge(mesh_plain, 4, 0)  # 不抛异常


class TestSplitRejections(unittest.TestCase):
    def test_edge_not_in_mesh(self):
        # ring4 的两个不相邻链环顶点从不同属一个单元。
        verts, tets = ring_star_vertices_tets(4)
        mesh = TetMesh.build(verts, tets)
        with self.assertRaises(EdgeSplitError) as cm:
            split_interior_edge(mesh, 2, 4)
        self.assertEqual(cm.exception.reason, "edge_not_found")
        self.assertEqual(cm.exception.edge, (2, 4))

    def test_non_manifold_edge_rejected(self):
        mesh = TetMesh.build(PINWHEEL_VERTICES, PINWHEEL_TETS)
        # 夹具本身无零体积单元：失败原因精确为非流形边。
        self.assertEqual(audit(mesh).zero_tets, ())
        with self.assertRaises(EdgeSplitError) as cm:
            split_interior_edge(mesh, 0, 1)
        self.assertEqual(cm.exception.reason, "non_manifold_edge")
        self.assertEqual(cm.exception.tet_indices, tuple(range(5)))
        self.assertIsNotNone(cm.exception.edge)

    def test_pinched_disconnected_link_rejected(self):
        mesh = TetMesh.build(PINCH_VERTICES, PINCH_TETS)
        self.assertEqual(audit(mesh).zero_tets, ())
        with self.assertRaises(EdgeSplitError) as cm:
            split_interior_edge(mesh, 0, 1)
        self.assertEqual(cm.exception.reason, "disconnected_link")
        self.assertEqual(cm.exception.tet_indices, tuple(range(6)))

    def test_zero_volume_star_rejected(self):
        mesh = TetMesh.build(DEGEN_STAR_VERTICES, DEGEN_STAR_TETS)
        with self.assertRaises(EdgeSplitError) as cm:
            split_interior_edge(mesh, 0, 1)
        self.assertEqual(cm.exception.reason, "zero_volume_star")
        self.assertEqual(cm.exception.tet_indices, (0,))

    def test_tet_limit_boundary(self):
        # 250 环：500 个单元恰好达到上限，接受；251 环：502，拒绝。
        v250, t250 = ring_star_vertices_tets(250)
        result = split_interior_edge(TetMesh.build(v250, t250), 0, 1)
        self.assertEqual(len(result.mesh.tets), MAX_TETS)

        v251, t251 = ring_star_vertices_tets(251)
        with self.assertRaises(EdgeSplitError) as cm:
            split_interior_edge(TetMesh.build(v251, t251), 0, 1)
        self.assertEqual(cm.exception.reason, "too_many_tets")

    def test_parameter_out_of_open_range(self):
        mesh = _bary_mesh()
        for bad in (0, 1, -1, Fraction(5, 4), Fraction(-1, 7)):
            with self.subTest(bad=bad):
                with self.assertRaises(EdgeSplitError) as cm:
                    split_interior_edge(mesh, 4, 0, bad)
                self.assertEqual(cm.exception.reason, "t_out_of_range")

    def test_parameter_float_rejected_with_location(self):
        mesh = _bary_mesh()
        with self.assertRaises(MeshInputError) as cm:
            split_interior_edge(mesh, 4, 0, 0.5)  # type: ignore[arg-type]
        self.assertEqual(cm.exception.location, "t")

    def test_parameter_bad_string_rejected(self):
        mesh = _bary_mesh()
        for bad in ("1/0", "nan", "inf", "abc"):
            with self.subTest(bad=bad):
                with self.assertRaises(MeshInputError) as cm:
                    split_interior_edge(mesh, 4, 0, bad)
                self.assertEqual(cm.exception.location, "t")

    def test_invalid_endpoints(self):
        mesh = _bary_mesh()
        n = len(BARY_VERTICES)
        cases = [
            (0, 0, "b", lambda: split_interior_edge(mesh, 0, 0)),
            (0, n, "b", lambda: split_interior_edge(mesh, 0, n)),
            (-1, 0, "a", lambda: split_interior_edge(mesh, -1, 0)),
            (True, 0, "a", lambda: split_interior_edge(mesh, True, 0)),  # type: ignore[arg-type]
            (0, 1.0, "b", lambda: split_interior_edge(mesh, 0, 1.0)),  # type: ignore[arg-type]
            ("4", 0, "a", lambda: split_interior_edge(mesh, "4", 0)),  # type: ignore[arg-type]
        ]
        for _a, _b, loc, call in cases:
            with self.subTest(loc=loc):
                with self.assertRaises(MeshInputError) as cm:
                    call()
                self.assertEqual(cm.exception.location, loc)

    def test_requires_tetmesh(self):
        with self.assertRaises(TypeError):
            split_interior_edge(BARY_VERTICES, 4, 0)  # type: ignore[arg-type]

    def test_failure_is_atomic(self):
        # 每一种失败后，原网格在审计意义下与失败前完全一致。
        mesh = _bary_mesh(labels=BARY_LABELS)
        before = audit(mesh)
        failing_calls = [
            lambda: split_interior_edge(mesh, 0, 1),       # 边界边
            lambda: split_interior_edge(mesh, 4, 0, 0),   # t=0
            lambda: split_interior_edge(mesh, 1, 3),       # 边界面上的边
        ]
        for call in failing_calls:
            with self.assertRaises(Exception):
                call()
            after = audit(mesh)
            self.assertEqual(after.total_signed_volume6, before.total_signed_volume6)
            self.assertEqual(len(after.boundary_faces), len(before.boundary_faces))
            self.assertEqual(mesh.tets, tuple(BARY_TETS))
            self.assertEqual(mesh.tet_labels, tuple(BARY_LABELS))
            self.assertEqual(len(mesh.vertices), len(BARY_VERTICES))


    def test_split_coexists_with_unrelated_non_manifold_face(self):
        # ring3 的帽面 {0,2,4} 原本 1 重；追加两个分居两侧的单元共用它，
        # 使该面变成 3 重非流形面（不含分裂边 0-1）。分裂内部边必须成功，
        # 且既不引入也不消除这个既有的非流形缺陷。
        rv, rt = ring_star_vertices_tets(3)
        verts = (*rv, (0, 0, 1), (0, 0, -1))  # 新点 5、6
        tets = (*rt, (0, 4, 2, 5), (0, 4, 6, 2))
        mesh = TetMesh.build(verts, tets)
        old_nm = {f.key for f in audit(mesh).non_manifold_faces}
        self.assertEqual(old_nm, {(0, 2, 4)})

        result = split_interior_edge(mesh, 0, 1, Fraction(1, 3))
        new_rep = audit(result.mesh)
        self.assertEqual({f.key for f in new_rep.non_manifold_faces}, old_nm)
        # 体积与边界面集合仍精确守恒（即便输入含若干负向单元）。
        self.assertEqual(new_rep.total_signed_volume6, audit(mesh).total_signed_volume6)
        self.assertEqual(
            result.mesh.boundary_face_keys(), mesh.boundary_face_keys()
        )


class TestTetMeshLabelsAndProtectedFaces(unittest.TestCase):
    def test_default_labels_are_zero(self):
        mesh = TetMesh.build(BARY_VERTICES, BARY_TETS)
        self.assertEqual(mesh.tet_labels, (0, 0, 0, 0))
        self.assertEqual(mesh.protected_faces, ())

    def test_labels_length_mismatch(self):
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build(BARY_VERTICES, BARY_TETS, tet_labels=(1, 2))
        self.assertEqual(cm.exception.location, "tet_labels")

    def test_unhashable_label(self):
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build(
                BARY_VERTICES, BARY_TETS, tet_labels=([1], 0, 0, 0)
            )
        self.assertEqual(cm.exception.location, "tet_labels[0]")

    def test_protected_face_validation(self):
        bad_faces = [
            ([(0, 1)], "protected_faces[0]"),
            ([(0, 1, 1)], "protected_faces[0]"),
            ([(0, 1, 9)], "protected_faces[0][2]"),
            ([(0, 1, True)], "protected_faces[0][2]"),
        ]
        for faces, loc in bad_faces:
            with self.subTest(loc=loc):
                with self.assertRaises(MeshInputError) as cm:
                    TetMesh.build(BARY_VERTICES, BARY_TETS, protected_faces=faces)
                self.assertEqual(cm.exception.location, loc)

    def test_protected_faces_deduplicated_and_sorted(self):
        mesh = TetMesh.build(
            BARY_VERTICES,
            BARY_TETS,
            protected_faces=[(4, 0, 3), (0, 3, 4)],  # 同一面两种排列
        )
        self.assertEqual(mesh.protected_faces, ((0, 3, 4),))

    def test_repair_keeps_labels_and_protected_faces(self):
        # repair_orientation 不增删单元，附属数据随位置原样保留。
        mesh = _bary_mesh(
            labels=BARY_LABELS, protected_faces=[(0, 3, 4)]
        )
        fixed = repair_orientation(mesh)
        self.assertEqual(fixed.tet_labels, tuple(BARY_LABELS))
        self.assertEqual(fixed.protected_faces, ((0, 3, 4),))


if __name__ == "__main__":
    unittest.main()
