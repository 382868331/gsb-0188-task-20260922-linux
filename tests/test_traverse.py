"""线段穿越采样与积分的必要测试。

独立对照纪律：

- 仿射场总积分用 ``reference.ref_affine_integral``（解析原函数）核对；
- 区间并集用 ``reference.ref_clip_segment``（叉积外法向 Cyrus–Beck
  逐四面体裁线，独立代码路径）核对；
- 两者都不调用被测核心得到期望值。
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    InvalidTraversalMeshError,
    MeshInputError,
    TetMesh,
    repair_orientation,
    traverse_segment,
)

from fixtures import CUBE_TETS, CUBE_VERTICES
from reference import ref_affine_integral, ref_affine_value, ref_clip_segment

#: 仿射场 f(x,y,z) = 2x - y + 3z + 1/2。
AFFINE = (Fraction(2), Fraction(-1), Fraction(3), Fraction(1, 2))


def affine_values(vertices):
    return [ref_affine_value(AFFINE, v) for v in vertices]


def cube_mesh():
    return TetMesh.build(CUBE_VERTICES, CUBE_TETS)


def merged(intervals):
    """把 (lo, hi) 区间列表合并成不重叠的有序列表。"""

    out = []
    for lo, hi in sorted(intervals):
        if out and lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [tuple(x) for x in out]


def piece_intervals(result):
    return [(p.t0, p.t1) for p in result.pieces]


def gap_intervals(result):
    return [(g.t0, g.t1) for g in result.gaps]


def ref_union(mesh, p0, p1):
    """独立逐四面体裁线的区间并集。"""

    ivs = []
    for tet in mesh.tets:
        r = ref_clip_segment(mesh.vertices, tet, p0, p1)
        if r is not None and r[0] < r[1]:
            ivs.append(r)
    return merged(ivs)


class TestTraverseInterior(unittest.TestCase):
    def setUp(self):
        self.mesh = cube_mesh()
        self.p0 = (Fraction(1, 8), Fraction(1, 8), Fraction(1, 8))
        self.p1 = (Fraction(7, 8), Fraction(7, 8), Fraction(7, 8))
        self.result = traverse_segment(
            self.mesh, affine_values(CUBE_VERTICES), self.p0, self.p1
        )

    def test_pieces_cover_whole_segment_no_gaps(self):
        self.assertEqual(self.result.gaps, ())
        self.assertEqual(merged(piece_intervals(self.result)), [(0, 1)])
        # 分段有序、首尾相接、互不重叠。
        for a, b in zip(self.result.pieces, self.result.pieces[1:]):
            self.assertEqual(a.t1, b.t0)
        self.assertEqual(self.result.pieces[0].t0, 0)
        self.assertEqual(self.result.pieces[-1].t1, 1)

    def test_total_integral_matches_affine_analytic(self):
        expect = ref_affine_integral(AFFINE, self.p0, self.p1, 0, 1)
        self.assertEqual(self.result.total_integral, expect)
        self.assertEqual(
            sum(p.integral for p in self.result.pieces), self.result.total_integral
        )

    def test_endpoint_values_match_affine_field(self):
        def point(t):
            return tuple(
                self.p0[k] + t * (self.p1[k] - self.p0[k]) for k in range(3)
            )

        for piece in self.result.pieces:
            self.assertEqual(piece.value0, ref_affine_value(AFFINE, point(piece.t0)))
            self.assertEqual(piece.value1, ref_affine_value(AFFINE, point(piece.t1)))
            self.assertEqual(
                piece.integral,
                ref_affine_integral(AFFINE, self.p0, self.p1, piece.t0, piece.t1),
            )

    def test_interval_union_matches_independent_clipping(self):
        self.assertEqual(
            merged(piece_intervals(self.result)), ref_union(self.mesh, self.p0, self.p1)
        )

    def test_multiple_tets_crossed(self):
        self.assertGreaterEqual(len(set(self.result.tet_sequence)), 2)


class TestTraverseGapsAndDegenerate(unittest.TestCase):
    def setUp(self):
        self.mesh = cube_mesh()

    def test_diagonal_through_shared_vertices(self):
        # 线段恰过顶点 (0,0,0) 与 (1,1,1)：t=1/3 与 t=2/3。
        p0, p1 = (-1, -1, -1), (2, 2, 2)
        result = traverse_segment(self.mesh, affine_values(CUBE_VERTICES), p0, p1)
        self.assertEqual(gap_intervals(result), [(Fraction(0), Fraction(1, 3)),
                                                 (Fraction(2, 3), Fraction(1))])
        # 分段与间隙互补覆盖 [0,1]，同一长度不计两次。
        self.assertEqual(
            merged(piece_intervals(result) + gap_intervals(result)), [(0, 1)]
        )
        # 与独立裁线并集一致。
        self.assertEqual(
            merged(piece_intervals(result)), ref_union(self.mesh, p0, p1)
        )
        # 仿射场总积分 = 解析积分（网内部分 [1/3, 2/3]）。
        self.assertEqual(
            result.total_integral,
            ref_affine_integral(AFFINE, p0, p1, Fraction(1, 3), Fraction(2, 3)),
        )
        # 没有零长度分段（相切/过点不产生重复计数）。
        for p in result.pieces:
            self.assertLess(p.t0, p.t1)

    def test_segment_along_shared_edge_has_deterministic_owner(self):
        # 线段与棱 (0,0,0)-(1,1,0) 重合，该棱被单元 0、1、4 共享：
        # 确定规则 → 整段归属最小下标单元 0，只计一次。
        p0, p1 = (0, 0, 0), (1, 1, 0)
        result = traverse_segment(self.mesh, affine_values(CUBE_VERTICES), p0, p1)
        self.assertEqual(gap_intervals(result), [])
        self.assertEqual(piece_intervals(result), [(Fraction(0), Fraction(1))])
        self.assertEqual(result.tet_sequence, (0,))
        self.assertEqual(
            result.total_integral, ref_affine_integral(AFFINE, p0, p1, 0, 1)
        )

    def test_segment_fully_outside(self):
        result = traverse_segment(
            self.mesh, affine_values(CUBE_VERTICES), (2, 2, 2), (3, 3, 3)
        )
        self.assertEqual(result.pieces, ())
        self.assertEqual(gap_intervals(result), [(Fraction(0), Fraction(1))])
        self.assertEqual(result.total_integral, 0)

    def test_tangent_touch_produces_no_piece(self):
        # 终点恰好停在顶点 (0,0,0)：只相触，不产生分段。
        result = traverse_segment(
            self.mesh, affine_values(CUBE_VERTICES), (-1, 0, 0), (0, 0, 0)
        )
        self.assertEqual(result.pieces, ())
        self.assertEqual(gap_intervals(result), [(Fraction(0), Fraction(1))])

    def test_exit_and_reenter_nonconvex(self):
        # 两个互不相接的四面体：线段先穿一个、经间隙、再穿另一个。
        verts = (
            (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (3, 0, 0), (4, 0, 0), (3, 1, 0), (3, 0, 1),
        )
        tets = [(0, 1, 2, 3), (4, 5, 6, 7)]
        mesh = TetMesh.build(verts, tets)
        values = [1] * 8
        result = traverse_segment(mesh, values, (-1, 0, 0), (5, 0, 0))
        # 常值场 1：总积分 = 网内参数长度之和 = 1/6 + 1/6。
        self.assertEqual(result.total_integral, Fraction(1, 3))
        self.assertEqual(len(result.pieces), 2)
        self.assertEqual(len(result.gaps), 3)
        self.assertEqual(
            merged(piece_intervals(result) + gap_intervals(result)), [(0, 1)]
        )
        self.assertEqual(
            merged(piece_intervals(result)), ref_union(mesh, (-1, 0, 0), (5, 0, 0))
        )


class TestInvalidMeshRejected(unittest.TestCase):
    def test_negative_tet_rejected_then_repaired(self):
        broken = [list(t) for t in CUBE_TETS]
        broken[2] = [CUBE_TETS[2][0], CUBE_TETS[2][1], CUBE_TETS[2][3], CUBE_TETS[2][2]]
        mesh = TetMesh.build(CUBE_VERTICES, broken)
        with self.assertRaises(InvalidTraversalMeshError) as cm:
            traverse_segment(mesh, affine_values(CUBE_VERTICES), (0, 0, 0), (1, 1, 1))
        self.assertEqual(cm.exception.negative_tets, (2,))
        # 按原审计规则修复后即可穿越。
        fixed = repair_orientation(mesh)
        result = traverse_segment(
            fixed, affine_values(CUBE_VERTICES), (0, 0, 0), (1, 1, 1)
        )
        self.assertEqual(merged(piece_intervals(result)), [(Fraction(0), Fraction(1))])

    def test_zero_volume_tet_rejected(self):
        verts = (*CUBE_VERTICES, (1, 1, 0))
        mesh = TetMesh.build(verts, [*CUBE_TETS, (0, 1, 2, 8)])
        with self.assertRaises(InvalidTraversalMeshError) as cm:
            traverse_segment(mesh, [0] * 9, (0, 0, 0), (1, 1, 1))
        self.assertEqual(cm.exception.zero_tets, (5,))

    def test_non_manifold_rejected(self):
        verts = ((0, 0, 0), (1, 1, 0), (1, 0, 1), (1, 0, 0), (0, 1, 1), (0, 1, 0))
        mesh = TetMesh.build(verts, [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5)])
        with self.assertRaises(InvalidTraversalMeshError) as cm:
            traverse_segment(mesh, [0] * 6, (0, 0, 0), (1, 1, 1))
        self.assertEqual(cm.exception.non_manifold_keys, ((0, 1, 2),))


class TestTraverseInputValidation(unittest.TestCase):
    def setUp(self):
        self.mesh = cube_mesh()

    def test_values_length_mismatch(self):
        with self.assertRaises(MeshInputError) as cm:
            traverse_segment(self.mesh, [1, 2, 3], (0, 0, 0), (1, 1, 1))
        self.assertEqual(cm.exception.location, "values")

    def test_float_value_rejected(self):
        values = [0] * 8
        values[3] = 0.5
        with self.assertRaises(MeshInputError) as cm:
            traverse_segment(self.mesh, values, (0, 0, 0), (1, 1, 1))
        self.assertEqual(cm.exception.location, "values[3]")

    def test_float_endpoint_rejected(self):
        with self.assertRaises(MeshInputError) as cm:
            traverse_segment(self.mesh, [0] * 8, (0, 0, 0.1), (1, 1, 1))
        self.assertEqual(cm.exception.location, "start[2]")

    def test_bool_value_rejected(self):
        values = [0] * 8
        values[0] = True
        with self.assertRaises(MeshInputError):
            traverse_segment(self.mesh, values, (0, 0, 0), (1, 1, 1))

    def test_non_mesh_rejected(self):
        with self.assertRaises(TypeError):
            traverse_segment("not a mesh", [0] * 8, (0, 0, 0), (1, 1, 1))


if __name__ == "__main__":
    unittest.main()
