"""必要测试：显式重合点合并。

覆盖：
- 只有 merge_pairs 明确给出且坐标精确相等的点才会合并；
- 未列出的重合点保持独立，很近但不相等的配对被拒绝（无自动吸附）；
- 新旧顶点 / 单元编号映射可双向核对；
- 合并后复查重复单元、零体积、非流形，失败时两份原网格保持不变；
- 合并后单元总数仍受 MAX_TETS 约束；负向单元不阻止合并。
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    MAX_TETS,
    DuplicateTetrahedraError,
    MeshInputError,
    MergeResult,
    NonManifoldError,
    TetMesh,
    ZeroVolumeError,
    merge_meshes,
)

# 第一份：正四面体 (0,0,0),(1,0,0),(0,1,0),(0,0,1)。
A_VERTS = (
    (0, 0, 0),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)
A_TETS = ((0, 1, 2, 3),)

# 第二份：前三个点与第一份的 1,2,3 坐标完全相同，第四点在对角内侧。
B_VERTS = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 1),
)
B_TETS = ((0, 1, 2, 3),)

# (第一份下标, 第二份下标)：共享面三点。
PAIRS = ((1, 0), (2, 1), (3, 2))


def _build_pair(verts_a=A_VERTS, tets_a=A_TETS, verts_b=B_VERTS, tets_b=B_TETS):
    return TetMesh.build(verts_a, tets_a), TetMesh.build(verts_b, tets_b)


class TestExplicitMerge(unittest.TestCase):
    def test_only_explicit_exact_coincident_points_merge(self):
        first, second = _build_pair()
        result = merge_meshes(first, second, PAIRS)
        self.assertIsInstance(result, MergeResult)

        mesh = result.mesh
        self.assertEqual(len(mesh.vertices), 5)
        self.assertEqual(mesh.tets, ((0, 1, 2, 3), (1, 2, 3, 4)))
        self.assertEqual(mesh.vertices[4], (Fraction(1), Fraction(1), Fraction(1)))

        mapping = result.mapping
        self.assertEqual(mapping.first_vertex_map, (0, 1, 2, 3))
        self.assertEqual(mapping.second_vertex_map, (1, 2, 3, 4))
        self.assertEqual(
            mapping.merged_vertices, ((1, 1, 0), (2, 2, 1), (3, 3, 2))
        )
        self.assertEqual(mapping.vertex_origins[1], ((0, 1), (1, 0)))
        self.assertEqual(mapping.vertex_origins[4], ((1, 3),))
        self.assertEqual(mapping.first_tet_map, (0,))
        self.assertEqual(mapping.second_tet_map, (1,))
        self.assertEqual(mapping.tet_origins, ((0, 0), (1, 0)))

        # 合并后是两单元共面的合法流形，内面诱导方向相反。
        report = result.report
        self.assertTrue(report.is_manifold)
        self.assertEqual(report.zero_tets, ())
        self.assertEqual(len(report.interior_faces), 1)
        self.assertTrue(report.interior_faces[0].opposite)
        self.assertEqual(report.components, ((0, 1),))

    def test_unlisted_coincident_points_stay_separate(self):
        first, second = _build_pair()
        result = merge_meshes(first, second, ())
        # 即便有三个点坐标完全相同，未显式配对也必须各自保留。
        self.assertEqual(len(result.mesh.vertices), 8)
        self.assertEqual(len(result.mesh.tets), 2)
        self.assertEqual(result.mapping.merged_vertices, ())
        # 没有任何点被识别为连接，两个单元互不相连。
        self.assertEqual(result.report.components, ((0,), (1,)))

    def test_near_but_not_equal_pair_is_refused_without_snapping(self):
        almost_verts = (
            (1, 0, Fraction(1, 10**9)),  # 与 A 的 (1,0,0) 极近但不相等
            (0, 1, 0),
            (0, 0, 1),
            (1, 1, 1),
        )
        first = TetMesh.build(A_VERTS, A_TETS)
        second = TetMesh.build(almost_verts, B_TETS)
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, PAIRS)
        self.assertEqual(cm.exception.location, "merge_pairs[0]")
        self.assertIn("不精确相等", str(cm.exception))

    def test_duplicate_pair_index_is_refused(self):
        first, second = _build_pair()
        # 第一份的 1 被配对两次（一对一约束）。
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, ((1, 0), (1, 1)))
        self.assertEqual(cm.exception.location, "merge_pairs[1]")
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, ((1, 0), (2, 0)))
        self.assertEqual(cm.exception.location, "merge_pairs[1]")

    def test_out_of_range_and_bad_type_pairs_are_located(self):
        first, second = _build_pair()
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, ((9, 0),))
        self.assertEqual(cm.exception.location, "merge_pairs[0][0]")
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, ((1, 9),))
        self.assertEqual(cm.exception.location, "merge_pairs[0][1]")
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, ((1, True),))
        self.assertEqual(cm.exception.location, "merge_pairs[0][1]")
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, [(1,)])
        self.assertEqual(cm.exception.location, "merge_pairs[0]")

    def test_negative_orientation_does_not_block_merge(self):
        # 第二份四面体整体反向：共享面诱导同向，但这是定向问题不是拒绝条件。
        flipped_b = ((0, 3, 2, 1),)
        first = TetMesh.build(A_VERTS, A_TETS)
        second = TetMesh.build(B_VERTS, flipped_b)
        result = merge_meshes(first, second, PAIRS)
        self.assertEqual(result.report.negative_tets, (1,))
        self.assertFalse(result.report.interior_faces[0].opposite)


class TestMergeRechecks(unittest.TestCase):
    def setUp(self):
        self.first, self.second = _build_pair()

    def _snapshot(self, mesh):
        return mesh.vertices, mesh.tets

    def test_duplicate_tetrahedra_after_merge_rejected(self):
        # 第二份唯一的四面体配对后与第一份四点集合完全相同。
        dup_b_tets = ((0, 2, 1, 3),)  # 同四点集合的排列，构造单网格时合法
        first = TetMesh.build(A_VERTS, ((0, 1, 2, 3),))
        second = TetMesh.build(
            ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
            dup_b_tets,
        )
        before_a, before_b = self._snapshot(first), self._snapshot(second)
        with self.assertRaises(DuplicateTetrahedraError):
            merge_meshes(first, second, ((0, 0), (1, 1), (2, 2), (3, 3)))
        self.assertEqual(self._snapshot(first), before_a)
        self.assertEqual(self._snapshot(second), before_b)

    def test_zero_volume_after_merge_rejected_atomically(self):
        # 第二份新增一个四点共面（零体积）四面体；共享面三点照配。
        b_verts = (
            (1, 0, 0),  # 0
            (0, 1, 0),  # 1
            (0, 0, 1),  # 2
            (1, 1, 1),  # 3
            (1, 0, 0),  # 4 与 0 坐标相同（不同编号，单网格允许）
            (2, 0, 0),  # 5
        )
        # 退化体使用 0,4,5 三点共线再加第 1 点 -> 四点共面 z=0。
        deg_tets = (B_TETS[0], (0, 4, 5, 1))
        first = TetMesh.build(A_VERTS, A_TETS)
        second = TetMesh.build(b_verts, deg_tets)
        before_a, before_b = self._snapshot(first), self._snapshot(second)
        with self.assertRaises(ZeroVolumeError) as cm:
            merge_meshes(first, second, PAIRS)
        # 第一份的 1 个单元占用新下标 0，因此第二份退化体新下标为 2。
        self.assertEqual(cm.exception.tet_indices, (2,))
        self.assertEqual(self._snapshot(first), before_a)
        self.assertEqual(self._snapshot(second), before_b)

    def test_non_manifold_after_merge_rejected(self):
        # 第二份含两个共享同一面的四面体；与第一份配对后该面被引用 3 次。
        b_verts = (
            (1, 0, 0),  # 0 共享
            (0, 1, 0),  # 1 共享
            (0, 0, 1),  # 2 共享
            (1, 1, 1),  # 3
            (-1, -1, -1),  # 4
        )
        b_tets = ((0, 1, 2, 3), (0, 1, 2, 4))
        first = TetMesh.build(A_VERTS, A_TETS)
        second = TetMesh.build(b_verts, b_tets)
        before_a, before_b = self._snapshot(first), self._snapshot(second)
        with self.assertRaises(NonManifoldError) as cm:
            merge_meshes(first, second, PAIRS)
        faces = cm.exception.faces
        self.assertEqual(len(faces), 1)
        key, tet_indices = faces[0]
        self.assertEqual(key, (1, 2, 3))
        self.assertEqual(tet_indices, (0, 1, 2))
        self.assertEqual(self._snapshot(first), before_a)
        self.assertEqual(self._snapshot(second), before_b)

    def test_total_tet_limit_enforced_for_merged_result(self):
        verts_a = [(i, 0, 0) for i in range(4 * 251)]
        tets_a = [(4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3) for i in range(251)]
        verts_b = [(i, 10, 0) for i in range(4 * 251)]
        tets_b = [(4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3) for i in range(251)]
        first = TetMesh.build(verts_a, tets_a)
        second = TetMesh.build(verts_b, tets_b)
        with self.assertRaises(MeshInputError) as cm:
            merge_meshes(first, second, ())
        self.assertEqual(cm.exception.location, "tets")


if __name__ == "__main__":
    unittest.main()
