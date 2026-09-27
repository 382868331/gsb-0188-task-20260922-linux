"""必要测试：显式点对合并、拒绝自动吸附、合并后复检与原子失败。"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    DuplicateTetrahedraError,
    MergePointError,
    TetMesh,
    ZeroVolumeError,
    audit,
    merge_meshes,
)

# 两个单位立方体：A=[0,1]^3，B 为 A 关于 x=1 的镜像剖分（x∈[1,2]），
# 保证接合面 x=1 两侧三角剖分共形（同一条对角），且 B 已按正体积给出。
A_VERTS = (
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
)
A_TETS = (
    (0, 1, 2, 5),
    (0, 2, 3, 7),
    (0, 4, 5, 7),
    (2, 5, 6, 7),
    (0, 2, 7, 5),
)
B_VERTS = (
    (2, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 0),
    (2, 0, 1), (1, 0, 1), (1, 1, 1), (2, 1, 1),
)
B_TETS = (
    (0, 1, 5, 2),
    (0, 2, 7, 3),
    (0, 4, 7, 5),
    (2, 5, 7, 6),
    (0, 2, 5, 7),
)
# 接合面四点两侧编号相同：1,2,5,6。
FACE_PAIRS = ((1, 1), (2, 2), (5, 5), (6, 6))


def _build_pair():
    return TetMesh.build(A_VERTS, A_TETS), TetMesh.build(B_VERTS, B_TETS)


class TestMergeHappyPath(unittest.TestCase):
    def test_glued_box_geometry_and_maps(self):
        a, b = _build_pair()
        result = merge_meshes(a, b, FACE_PAIRS)

        # 新旧编号关系：A 恒等；B 共享点指向 A 编号，其余追加为 8..11。
        self.assertEqual(result.vertex_map_a, tuple(range(8)))
        self.assertEqual(result.vertex_map_b, (8, 1, 2, 9, 10, 5, 6, 11))
        self.assertEqual(result.merged_pairs, tuple(FACE_PAIRS))
        self.assertEqual(result.n_vertices_a, 8)
        self.assertEqual(len(result.mesh.vertices), 12)
        self.assertEqual(len(result.mesh.tets), 10)

        # 合并后重审：总体积 2、单一连通分量、无非流形。
        report = result.report
        self.assertEqual(report.total_signed_volume, Fraction(2))
        self.assertFalse(report.has_non_manifold)
        self.assertFalse(report.has_zero_volume)
        self.assertEqual(report.components, (tuple(range(10)),))
        # 接合面两个三角成为新内面；连同各自 4 个内部面共 10 个内面。
        self.assertEqual(len(report.boundary_faces), 20)

        # 与直接对结果网格重新审计一致（报告不是缓存的旧结果）。
        self.assertEqual(audit(result.mesh), report)

        # 反查新编号来源。
        self.assertEqual(result.origin_of(2), ("a", 2))
        self.assertEqual(result.origin_of(9), ("b", 3))

    def test_inputs_not_mutated_and_vertices_reused(self):
        a, b = _build_pair()
        snap = (a.vertices, a.tets, b.vertices, b.tets)
        result = merge_meshes(a, b, FACE_PAIRS)
        # 两份原网格顶点/单元完全不动（失败或成功都不回写）。
        self.assertEqual((a.vertices, a.tets, b.vertices, b.tets), snap)
        # A 顶点对象原样进入新网格。
        self.assertTrue(
            all(x is y for x, y in zip(result.mesh.vertices[:8], a.vertices, strict=True))
        )


class TestNoAutoSnap(unittest.TestCase):
    def test_coincident_but_unlisted_points_stay_separate(self):
        # 一对点都不列：8+8 个顶点，拼接处仍是边界（两份独立网格）。
        a, b = _build_pair()
        result = merge_meshes(a, b, ())
        self.assertEqual(len(result.mesh.vertices), 16)
        self.assertEqual(result.vertex_map_b, tuple(range(8, 16)))
        self.assertEqual(len(result.report.components), 2)
        # 未合并的重合点不会被吸在一起：不存在跨越两份网格（单元 0..4 与
        # 5..9）的内面；各自的 4 个剖分内面仍在，但属于两个独立分量。
        cross = [
            f for f in result.report.interior_faces
            if f.tet_indices[0] < 5 <= f.tet_indices[1]
        ]
        self.assertEqual(cross, [])
        self.assertEqual(result.report.total_signed_volume, Fraction(2))

    def test_partial_pairs_only_those_points_merged(self):
        # 只合并两个点（棱端点）：其余重合点仍独立。
        a, b = _build_pair()
        result = merge_meshes(a, b, ((1, 1), (2, 2)))
        self.assertEqual(result.vertex_map_b, (8, 1, 2, 9, 10, 11, 12, 13))
        self.assertEqual(len(result.mesh.vertices), 14)

    def test_position_mismatch_pair_rejected_no_snapping(self):
        a, b = _build_pair()
        before = (a.vertices, a.tets, b.vertices, b.tets)
        with self.assertRaises(MergePointError) as cm:
            merge_meshes(a, b, ((1, 0),))  # A[1]=(1,0,0) vs B[0]=(2,0,0)
        self.assertEqual(cm.exception.location, "merge_pairs[0]")
        self.assertIn("不按距离自动吸附", str(cm.exception))
        # 失败后输入原样。
        self.assertEqual((a.vertices, a.tets, b.vertices, b.tets), before)


class TestMergePairValidation(unittest.TestCase):
    def test_bad_types_and_locations(self):
        a, b = _build_pair()
        cases = [
            ("not pairs", "merge_pairs"),
            ([(1,)], "merge_pairs[0]"),
            ([(1, 0, 2)], "merge_pairs[0]"),
            ([(True, 0)], "merge_pairs[0][0]"),
            ([(1, False)], "merge_pairs[0][1]"),
            ([(99, 0)], "merge_pairs[0][0]"),
            ([(1, 99)], "merge_pairs[0][1]"),
        ]
        for pairs, loc in cases:
            with self.subTest(pairs=pairs):
                with self.assertRaises(MergePointError) as cm:
                    merge_meshes(a, b, pairs)
                self.assertEqual(cm.exception.location, loc)

    def test_repeated_index_rejected(self):
        a, b = _build_pair()
        with self.assertRaises(MergePointError) as cm:
            merge_meshes(a, b, ((1, 1), (1, 2)))  # A 点 1 重复
        self.assertEqual(cm.exception.location, "merge_pairs[1][0]")
        with self.assertRaises(MergePointError) as cm:
            merge_meshes(a, b, ((1, 1), (2, 1)))  # B 点 1 重复
        self.assertEqual(cm.exception.location, "merge_pairs[1][1]")

    def test_non_mesh_arguments(self):
        a, _ = _build_pair()
        with self.assertRaises(TypeError):
            merge_meshes(A_VERTS, a, FACE_PAIRS)


class TestPostMergeRechecks(unittest.TestCase):
    def test_cross_mesh_duplicate_tet_rejected(self):
        # A/B 同位置同剖分（B 直接取 A 的坐标），合并四个共享面点后，
        # 跨网格出现完全相同的单元 -> 重复单元复检必须拦下。
        a = TetMesh.build(A_VERTS, A_TETS)
        b = TetMesh.build(A_VERTS, A_TETS)
        pairs = tuple((i, i) for i in range(8))  # 全部点一一对应
        before = (a.vertices, a.tets, b.vertices, b.tets)
        with self.assertRaises(DuplicateTetrahedraError) as cm:
            merge_meshes(a, b, pairs)
        # 5 组重复，且原网格不变。
        self.assertEqual(len(cm.exception.duplicates), 5)
        self.assertEqual((a.vertices, a.tets, b.vertices, b.tets), before)

    def test_zero_volume_after_merge_rejected(self):
        a = TetMesh.build(A_VERTS, A_TETS)
        # B 合法 + 一个共面零体积单元；合并点不涉及退化点。
        # B 的 0,2,3 都在 z=0；追加 (2,2,0) 与之共面。
        b_verts = (*B_VERTS, (2, 2, 0))
        b_tets = (*B_TETS, (0, 2, 3, 8))
        b = TetMesh.build(b_verts, b_tets)
        before = (a.vertices, a.tets, b.vertices, b.tets)
        with self.assertRaises(ZeroVolumeError) as cm:
            merge_meshes(a, b, FACE_PAIRS)
        # 零体积单元是 A 的 5 个之后的第 6 个，新编号 10。
        self.assertEqual(cm.exception.tet_indices, (10,))
        self.assertEqual((a.vertices, a.tets, b.vertices, b.tets), before)

    def test_non_manifold_after_merge_is_reported_not_raised(self):
        # 两个单四面体共享同一个三点面：合法 2 个单元；再加第三个单元
        # （来自 B 的第二个四面体）也引用同一面 -> 合并后非流形被报告。
        verts_a = (
            (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
        )
        tets_a = ((0, 1, 2, 3),)
        verts_b = (
            (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, -1), (0, 0, -2),
        )
        tets_b = ((0, 1, 2, 3), (0, 1, 2, 4))
        a = TetMesh.build(verts_a, tets_a)
        b = TetMesh.build(verts_b, tets_b)
        result = merge_meshes(a, b, ((0, 0), (1, 1), (2, 2)))
        self.assertEqual(len(result.mesh.tets), 3)
        nm = result.report.non_manifold_faces
        self.assertEqual(len(nm), 1)
        self.assertEqual(nm[0].key, (0, 1, 2))
        self.assertEqual(nm[0].tet_indices, (0, 1, 2))
        self.assertFalse(result.report.is_manifold)


if __name__ == "__main__":
    unittest.main()
