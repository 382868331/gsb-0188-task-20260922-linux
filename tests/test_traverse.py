"""线段穿越采样与积分的必要测试（第二轮新增）。

独立核对纪律：
- 区间并集用 tests/reference.py 的 ref_clip_intervals（逐四面体独立
  解半空间不等式，不经面邻接行走）核对；
- 总积分用 ref_affine_integral（仿射场解析积分）核对；
- 两个参考都不 import orient_audit。
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from orient_audit import (
    MeshInputError,
    TetMesh,
    TraversalError,
    repair_orientation,
    sample_segment,
)

from fixtures import CUBE_TETS, CUBE_VERTICES
from reference import (
    ref_affine_integral,
    ref_clip_intervals,
    ref_merge_intervals,
)

F = Fraction

# 仿射场 f(x,y,z) = 2x + 3y - z + 5 在立方体顶点上的精确取值。
ALPHA = (F(2), F(3), F(-1))
BETA = F(5)
CUBE_VALUES = tuple(ALPHA[0] * x + ALPHA[1] * y + ALPHA[2] * z + BETA
                    for x, y, z in CUBE_VERTICES)


def _cube_mesh():
    return TetMesh.build(CUBE_VERTICES, CUBE_TETS)


def _assert_consistent_with_reference(case, mesh, values, start, end):
    """公共断言：区间并集、间隙互补、总积分都与独立参考一致。"""

    tr = sample_segment(mesh, values, start, end)

    # 1) 命中段升序、互不重叠（同一长度不会计两次）。
    for h1, h2 in zip(tr.hits, tr.hits[1:]):
        case.assertLess(h1.t0, h1.t1)
        case.assertLessEqual(h1.t1, h2.t0)

    # 2) 命中区间必须精确等于独立逐四面体裁线结果的并集。
    ref = ref_clip_intervals(mesh.vertices, mesh.tets, start, end)
    ref_union = ref_merge_intervals([(lo, hi) for _, lo, hi in ref])
    hit_union = ref_merge_intervals([(h.t0, h.t1) for h in tr.hits])
    case.assertEqual(hit_union, ref_union)

    # 3) 每个命中段确实属于它声明的单元（参考中该单元覆盖该区间）。
    for h in tr.hits:
        case.assertTrue(
            any(ti == h.tet_index and lo <= h.t0 and h.t1 <= hi
                for ti, lo, hi in ref),
            f"hit {h} 不在单元 {h.tet_index} 的参考区间内",
        )

    # 4) 间隙与命中段互补地覆盖 [0, 1]（间隙不补 0、不混入积分）。
    events = [(h.t0, h.t1) for h in tr.hits] + [(g.t0, g.t1) for g in tr.gaps]
    merged = ref_merge_intervals(events)
    case.assertEqual(merged, [(F(0), F(1))])

    # 5) 总积分 == 仿射场在命中区间并集上的解析积分。
    expect = ref_affine_integral(ALPHA, BETA, start, end, hit_union)
    case.assertEqual(tr.total_integral, expect)
    return tr


class TestTraversal(unittest.TestCase):
    def test_crossing_sequence_and_values(self):
        # 斜穿立方体：起点在网格外，依次进入中心四面体再到角四面体。
        tr = _assert_consistent_with_reference(
            self, _cube_mesh(), CUBE_VALUES,
            (F(-1), F(-1), F(-1)), (F(2), F(2), F(2)),
        )
        self.assertEqual([h.tet_index for h in tr.hits], [4, 3])
        self.assertEqual([(h.t0, h.t1) for h in tr.hits],
                         [(F(1, 3), F(5, 9)), (F(5, 9), F(2, 3))])
        self.assertEqual([(g.t0, g.t1) for g in tr.gaps],
                         [(F(0), F(1, 3)), (F(2, 3), F(1))])
        # 端点值 == 仿射场在 P(t) 处的精确值。
        for h in tr.hits:
            for t, v in ((h.t0, h.value0), (h.t1, h.value1)):
                x = F(-1) + 3 * t
                self.assertEqual(v, ALPHA[0] * x + ALPHA[1] * x + ALPHA[2] * x + BETA)

    def test_axis_segment_inside_single_tet(self):
        # y=z=1/2 的线段全程落在中心四面体内（两端在共享棱上）。
        _assert_consistent_with_reference(
            self, _cube_mesh(), CUBE_VALUES,
            (F(-1), F(1, 2), F(1, 2)), (F(2), F(1, 2), F(1, 2)),
        )

    def test_through_shared_vertex_exactly(self):
        # 线恰过顶点 0：进入点同时属于 4 个单元，确定规则唯一归属。
        tr = _assert_consistent_with_reference(
            self, _cube_mesh(), CUBE_VALUES,
            (F(-1), F(0), F(0)), (F(1), F(1), F(1)),
        )
        self.assertEqual(tr.hits[0].t0, F(1, 2))

    def test_along_boundary_edge_single_owner(self):
        # 沿立方体棱 0->1（边界面与棱上行走）：归属唯一、长度不计两次。
        tr = _assert_consistent_with_reference(
            self, _cube_mesh(), CUBE_VALUES,
            (F(0), F(0), F(0)), (F(1), F(0), F(0)),
        )
        self.assertEqual(len(tr.hits), 1)
        self.assertEqual(tr.hits[0].tet_index, 0)
        self.assertEqual(sum(h.t1 - h.t0 for h in tr.hits), F(1))

    def test_along_boundary_face_diagonal(self):
        # 沿底面 z=0 的对角线（两个边界面的共享棱）：仍唯一归属。
        tr = _assert_consistent_with_reference(
            self, _cube_mesh(), CUBE_VALUES,
            (F(0), F(0), F(0)), (F(1), F(1), F(0)),
        )
        self.assertEqual(len(tr.hits), 1)
        self.assertEqual(sum(h.t1 - h.t0 for h in tr.hits), F(1))

    def test_fully_outside_is_single_gap(self):
        tr = sample_segment(_cube_mesh(), CUBE_VALUES,
                            (F(2), F(2), F(2)), (F(3), F(3), F(3)))
        self.assertEqual(tr.hits, ())
        self.assertEqual([(g.t0, g.t1) for g in tr.gaps], [(F(0), F(1))])
        self.assertEqual(tr.total_integral, F(0))

    def test_gap_not_zero_filled(self):
        # 网格外区间只出现在 gaps，积分只来自命中段。
        tr = _assert_consistent_with_reference(
            self, _cube_mesh(), CUBE_VALUES,
            (F(1, 2), F(1, 2), F(-1)), (F(1, 2), F(1, 2), F(2)),
        )
        self.assertTrue(tr.gaps)
        self.assertEqual(tr.total_integral,
                         sum((h.integral for h in tr.hits), F(0)))

    def test_integral_matches_affine_reference_many_directions(self):
        # 多方向扫描：每条都与解析积分 + 独立裁线核对。
        mesh = _cube_mesh()
        starts = [(F(-1), F(1, 4), F(1, 3)), (F(1, 2), F(-1), F(1, 5)),
                  (F(2), F(0), F(1)), (F(-1), F(-1), F(0))]
        ends = [(F(2), F(1, 4), F(1, 3)), (F(1, 2), F(2), F(1, 5)),
                (F(-1), F(1), F(0)), (F(2), F(2), F(1))]
        for a, b in zip(starts, ends, strict=True):
            with self.subTest(start=a, end=b):
                _assert_consistent_with_reference(self, mesh, CUBE_VALUES, a, b)


class TestTraversalRejectsInvalid(unittest.TestCase):
    def test_negative_tet_rejected_until_repaired(self):
        flipped = (CUBE_TETS[0], CUBE_TETS[1], CUBE_TETS[2], CUBE_TETS[3],
                   (0, 2, 5, 7))  # 中心单元交换最后两顶点 -> 负定向
        mesh = TetMesh.build(CUBE_VERTICES, flipped)
        with self.assertRaises(TraversalError) as cm:
            sample_segment(mesh, CUBE_VALUES, (F(0), F(0), F(0)), (F(1), F(1), F(1)))
        self.assertIn("tets[4]", str(cm.exception))
        # 按原审计规则修复后即可穿越。
        repaired = repair_orientation(mesh)
        tr = sample_segment(repaired, CUBE_VALUES, (F(0), F(0), F(0)), (F(1), F(1), F(1)))
        self.assertTrue(tr.hits)

    def test_zero_volume_tet_rejected(self):
        tets = CUBE_TETS + ((0, 1, 2, 3),)  # 共面退化
        mesh = TetMesh.build(CUBE_VERTICES, tets)
        with self.assertRaises(TraversalError) as cm:
            sample_segment(mesh, CUBE_VALUES, (F(0), F(0), F(0)), (F(1), F(1), F(1)))
        self.assertIn("tets[5]", str(cm.exception))

    def test_non_manifold_rejected(self):
        # 三个单元共享同一三点面 (0,2,5)。
        tets = CUBE_TETS + ((0, 2, 5, 6),)
        mesh = TetMesh.build(CUBE_VERTICES, tets)
        with self.assertRaises(TraversalError):
            sample_segment(mesh, CUBE_VALUES, (F(0), F(0), F(0)), (F(1), F(1), F(1)))


class TestTraversalInputValidation(unittest.TestCase):
    def test_values_length_mismatch(self):
        with self.assertRaises(MeshInputError) as cm:
            sample_segment(_cube_mesh(), (1, 2, 3), (F(0),) * 3, (F(1),) * 3)
        self.assertEqual(cm.exception.location, "values")

    def test_float_value_rejected(self):
        bad = list(CUBE_VALUES)
        bad[0] = 0.5
        with self.assertRaises(MeshInputError) as cm:
            sample_segment(_cube_mesh(), bad, (F(0),) * 3, (F(1),) * 3)
        self.assertEqual(cm.exception.location, "values[0]")

    def test_float_endpoint_rejected(self):
        with self.assertRaises(MeshInputError) as cm:
            sample_segment(_cube_mesh(), CUBE_VALUES, (0.0, 0, 0), (1, 1, 1))
        self.assertEqual(cm.exception.location, "start[0]")

    def test_coincident_endpoints_rejected(self):
        with self.assertRaises(MeshInputError) as cm:
            sample_segment(_cube_mesh(), CUBE_VALUES, (F(0),) * 3, (F(0),) * 3)
        self.assertEqual(cm.exception.location, "end")

    def test_non_mesh_rejected(self):
        with self.assertRaises(TypeError):
            sample_segment("not a mesh", CUBE_VALUES, (F(0),) * 3, (F(1),) * 3)


if __name__ == "__main__":
    unittest.main()

