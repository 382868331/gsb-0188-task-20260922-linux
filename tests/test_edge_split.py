"""内部边分裂测试：正常 / 边界 / 失败用例。

期望值来源纪律：体积用 tests/reference.py 的 4x4 Bareiss 独立行列式
（不调用被测核心），总体积另用面通量公式独立核对；夹具体积均先手算
（见各注释），不调用被测实现生成期望值。
"""

from __future__ import annotations

import itertools
import random
import unittest
from fractions import Fraction

from orient_audit import (
    EdgeSplitError,
    MeshInputError,
    TetMesh,
    audit,
    repair_orientation,
    split_internal_edge,
)
from orient_audit.mesh import build_face_incidence
from fixtures import CUBE_VERTICES, KUHN_LINK_KEYS, KUHN_TETS
from reference import ref_flux_volume, ref_signed_volume6


def _vol6(vertices, tet):
    """独立参考（4x4 Bareiss）的 6 倍有向体积。"""

    return ref_signed_volume6(vertices, tet)


def _ring_mesh(n_radii):
    """绕边 (0,1) 的 n_radii 四面体闭合环（确定性生成）。

    顶点 0=(0,0,0)、1=(0,0,1)，半径点 r_k=(k+1, (k+1)^2, 0)。
    单元 k（k=0..n-2）为 (0,1,r_k,r_{k+1})，闭合单元为 (0,1,r_{n-1},r_0)。
    手算：det(单元 k) = x_k·y_{k+1} - y_k·x_{k+1} = (k+1)(k+2) > 0；
    闭合单元 det = n - n^2 < 0。侧面 (0,1,r_k) 各被恰两个单元引用，
    构成单一闭合环，边 (0,1) 是内部边。
    """

    vertices = [(0, 0, 0), (0, 0, 1)] + [
        (k + 1, (k + 1) ** 2, 0) for k in range(n_radii)
    ]
    tets = [(0, 1, 2 + k, 3 + k) for k in range(n_radii - 1)]
    tets.append((0, 1, n_radii + 1, 2))
    return TetMesh.build(vertices, tets)


def _ring_mesh_volume6(n):
    """_ring_mesh 总 6 倍有向体积的手算闭式：sum j(j+1), j=1..n-1，加 n-n^2。"""

    return Fraction((n - 1) * n * (n + 1), 3) + (n - n * n)


#: 三环 + 一个互不相邻的孤立单元（用于验证非星单元原样保留）。
RING3_VERTICES = (
    (0, 0, 0),   # 0  边端点 a
    (0, 0, 1),   # 1  边端点 b
    (1, 1, 0),   # 2  半径 r0
    (2, 4, 0),   # 3  半径 r1
    (3, 9, 0),   # 4  半径 r2
    (10, 10, 10),  # 5..8 孤立单元
    (11, 10, 10),
    (10, 11, 10),
    (10, 10, 11),
)
# 手算 6 倍体积：T0=2, T1=6, T2=6（行列式 z 分量），T3=1；合计 15。
RING3_TETS = (
    (0, 1, 2, 3),
    (0, 1, 3, 4),
    (0, 1, 2, 4),
    (5, 6, 7, 8),
)


class TestKuhnFixtureSanity(unittest.TestCase):
    """夹具自检：Kuhn 6 剖分的解析性质（独立参考核对）。"""

    def setUp(self):
        self.mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        self.report = audit(self.mesh)

    def test_each_tet_volume_is_one_sixth_by_reference(self):
        for tet in KUHN_TETS:
            self.assertEqual(_vol6(CUBE_VERTICES, tet), Fraction(1))
        self.assertEqual(self.report.total_signed_volume, Fraction(1))

    def test_diagonal_is_the_unique_interior_edge_star(self):
        self.assertEqual(len(self.report.boundary_faces), 12)
        self.assertEqual(len(self.report.interior_faces), 6)
        keys = sorted(f.key for f in self.report.interior_faces)
        self.assertEqual(keys, sorted(KUHN_LINK_KEYS))
        self.assertTrue(all(f.opposite for f in self.report.interior_faces))
        self.assertTrue(self.report.consistently_oriented)
        self.assertTrue(self.report.is_manifold)
        self.assertEqual(self.report.components, (tuple(range(6)),))

    def test_boundary_flux_volume_is_one(self):
        triples = [f.ordered for f in self.report.boundary_faces]
        self.assertEqual(ref_flux_volume(CUBE_VERTICES, triples), Fraction(1))


class TestSplitNormal(unittest.TestCase):
    """正常分裂：体积份额、朝向继承、边界不动、标签继承、结果映射。"""

    def test_split_kuhn_diagonal_at_one_half(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        res = split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        m2 = res.mesh

        # 结果映射。
        self.assertEqual(res.edge, (0, 6))
        self.assertEqual(res.parameter, Fraction(1, 2))
        self.assertEqual(res.new_vertex_index, 8)
        self.assertEqual(res.star_tet_indices, (0, 1, 2, 3, 4, 5))
        self.assertEqual(
            res.child_tet_indices,
            ((0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11)),
        )
        # 新点坐标精确为对角线中点（手算）。
        self.assertEqual(
            m2.vertices[8],
            (Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)),
        )
        # 旧顶点对象恒等（不移动顶点）。
        self.assertTrue(
            all(x is y for x, y in zip(m2.vertices[:8], mesh.vertices))
        )

        # 每个父单元 6 倍体积 1（手算）→ 子单元各 1/2（独立参考核对）。
        for ia, ib in res.child_tet_indices:
            self.assertEqual(_vol6(m2.vertices, m2.tets[ia]), Fraction(1, 2))
            self.assertEqual(_vol6(m2.vertices, m2.tets[ib]), Fraction(1, 2))

        # 审计：12 单元、总量 1、边界 12 面不变、内面 6 -> 18、定向一致。
        rep = audit(m2)
        self.assertEqual(rep.n_tets, 12)
        self.assertEqual(rep.total_signed_volume, Fraction(1))
        self.assertEqual(rep.negative_tets, ())
        self.assertEqual(rep.zero_tets, ())
        self.assertEqual(len(rep.boundary_faces), 12)
        self.assertEqual(len(rep.interior_faces), 18)
        self.assertTrue(rep.consistently_oriented)
        self.assertTrue(rep.is_manifold)
        self.assertEqual(rep.components, (tuple(range(12)),))
        # 边界有序面多重集与分裂前完全一致。
        before = sorted(f.ordered for f in audit(mesh).boundary_faces)
        after = sorted(f.ordered for f in rep.boundary_faces)
        self.assertEqual(before, after)
        # 面通量公式独立核对总体积仍为 1。
        self.assertEqual(
            ref_flux_volume(m2.vertices, [f.ordered for f in rep.boundary_faces]),
            Fraction(1),
        )
        # 原网格未被修改。
        self.assertEqual(len(mesh.tets), 6)
        self.assertEqual(audit(mesh).total_signed_volume, Fraction(1))

    def test_split_at_one_third_gives_exact_rational_shares(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        res = split_internal_edge(mesh, (0, 6), Fraction(1, 3))
        m2 = res.mesh
        # 手算：父 6 倍体积 1 → 子单元 1/3 与 2/3（独立参考核对）。
        for ia, ib in res.child_tet_indices:
            self.assertEqual(_vol6(m2.vertices, m2.tets[ia]), Fraction(1, 3))
            self.assertEqual(_vol6(m2.vertices, m2.tets[ib]), Fraction(2, 3))
        self.assertEqual(audit(m2).total_signed_volume, Fraction(1))

    def test_parameter_accepts_rational_string(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        res = split_internal_edge(mesh, (0, 6), "2/5")
        self.assertEqual(res.parameter, Fraction(2, 5))
        m2 = res.mesh
        for ia, ib in res.child_tet_indices:
            self.assertEqual(_vol6(m2.vertices, m2.tets[ia]), Fraction(2, 5))
            self.assertEqual(_vol6(m2.vertices, m2.tets[ib]), Fraction(3, 5))

    def test_edge_endpoint_order_is_symmetric(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        forward = split_internal_edge(mesh, (0, 6), Fraction(3, 4)).mesh
        backward = split_internal_edge(mesh, (6, 0), Fraction(1, 4)).mesh
        # 同一几何点：t=3/4 从 0 出发 == t=1/4 从 6 出发。
        self.assertEqual(forward.vertices, backward.vertices)
        # 单元集合一致（子单元 A/B 角色互换，集合不变）。
        self.assertEqual(
            sorted(tuple(sorted(t)) for t in forward.tets),
            sorted(tuple(sorted(t)) for t in backward.tets),
        )

    def test_materials_inherited_by_children(self):
        mats = (10, 20, 30, 40, 50, 60)
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS, materials=mats)
        res = split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        expected = tuple(v for m in mats for v in (m, m))
        self.assertEqual(res.mesh.materials, expected)

    def test_default_materials_are_zero(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        self.assertEqual(mesh.materials, (0,) * 6)
        res = split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        self.assertEqual(res.mesh.materials, (0,) * 12)

    def test_protected_boundary_face_preserved(self):
        mesh = TetMesh.build(
            CUBE_VERTICES, KUHN_TETS, protected_faces=[(0, 1, 2)]
        )
        res = split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        m2 = res.mesh
        self.assertEqual(m2.protected_faces, ((0, 1, 2),))
        # 保护面仍是边界面且三点集合不含新顶点。
        boundary_keys = {
            tuple(sorted(f.ordered)) for f in audit(m2).boundary_faces
        }
        self.assertIn((0, 1, 2), boundary_keys)

    def test_non_star_tets_preserved_in_order(self):
        mats = (5, 6, 7, 8)
        mesh = TetMesh.build(RING3_VERTICES, RING3_TETS, materials=mats)
        res = split_internal_edge(mesh, (0, 1), Fraction(1, 4))
        m2 = res.mesh
        self.assertEqual(res.star_tet_indices, (0, 1, 2))
        self.assertEqual(res.child_tet_indices, ((0, 1), (2, 3), (4, 5)))
        # 孤立单元原样落在下标 6，标签不变。
        self.assertEqual(m2.tets[6], (5, 6, 7, 8))
        self.assertEqual(m2.materials, (5, 5, 6, 6, 7, 7, 8))
        # 手算份额：父 6 倍体积 2,6,6 → t=1/4 得 1/2,3/2 与 3/2,9/2 两对。
        shares = (
            (Fraction(1, 2), Fraction(3, 2)),
            (Fraction(3, 2), Fraction(9, 2)),
            (Fraction(3, 2), Fraction(9, 2)),
        )
        for (ia, ib), (wa, wb) in zip(res.child_tet_indices, shares, strict=True):
            self.assertEqual(_vol6(m2.vertices, m2.tets[ia]), wa)
            self.assertEqual(_vol6(m2.vertices, m2.tets[ib]), wb)
        self.assertEqual(_vol6(m2.vertices, m2.tets[6]), Fraction(1))
        # 总量手算 15/6 = 5/2，分裂前后守恒。
        self.assertEqual(audit(m2).total_signed_volume, Fraction(5, 2))
        # 连通分量：环与孤立单元仍然分离。
        self.assertEqual(audit(m2).components, ((0, 1, 2, 3, 4, 5), (6,)))

    def test_negative_parent_orientation_inherited(self):
        # 翻转 T1（交换最后两顶点）使其体积为负；分裂后两个子单元
        # 必须继承负号与精确份额，已有定向冲突被保留而非放大。
        tets = [list(t) for t in RING3_TETS]
        tets[1] = [0, 1, 4, 3]  # 负体积版本
        mesh = TetMesh.build(RING3_VERTICES, tets)
        before = audit(mesh)
        self.assertEqual(before.negative_tets, (1,))
        self.assertEqual(before.total_signed_volume, Fraction(1, 2))  # 手算 3/6

        res = split_internal_edge(mesh, (0, 1), Fraction(1, 4))
        m2 = res.mesh
        ia, ib = res.child_tet_indices[1]
        # 父 6 倍体积 -6 → 子单元 -3/2 与 -9/2（独立参考核对）。
        self.assertEqual(_vol6(m2.vertices, m2.tets[ia]), Fraction(-3, 2))
        self.assertEqual(_vol6(m2.vertices, m2.tets[ib]), Fraction(-9, 2))
        rep = audit(m2)
        self.assertEqual(rep.negative_tets, (ia, ib))
        self.assertEqual(rep.total_signed_volume, Fraction(1, 2))
        # 分裂前冲突侧面是 (0,1,2)（奇数环长固有）与 (0,1,3)（翻转 T1
        # 新增）；每个冲突侧面的两个半边分裂后仍冲突，一致的 (0,1,4)
        # 两个半边仍一致——分裂不制造也不修复冲突。
        inconsistent = {f.key for f in rep.inconsistent_faces}
        self.assertEqual(
            inconsistent,
            {(0, 2, 9), (1, 2, 9), (0, 3, 9), (1, 3, 9)},
        )
        consistent = {
            f.key for f in rep.interior_faces if f.opposite
        }
        self.assertIn((0, 4, 9), consistent)
        self.assertIn((1, 4, 9), consistent)

    def test_chained_splits(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        first = split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        second = split_internal_edge(first.mesh, (0, 8), Fraction(1, 2))
        m3 = second.mesh
        self.assertEqual(second.new_vertex_index, 9)
        self.assertEqual(
            m3.vertices[9],
            (Fraction(1, 4), Fraction(1, 4), Fraction(1, 4)),
        )
        rep = audit(m3)
        self.assertEqual(rep.n_tets, 18)
        self.assertEqual(rep.total_signed_volume, Fraction(1))
        self.assertEqual(len(rep.boundary_faces), 12)
        self.assertTrue(rep.consistently_oriented)
        self.assertTrue(rep.is_manifold)
        self.assertEqual(
            ref_flux_volume(m3.vertices, [f.ordered for f in rep.boundary_faces]),
            Fraction(1),
        )

    def test_split_result_mesh_supports_repair_and_audit(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        m2 = split_internal_edge(mesh, (0, 6), Fraction(1, 2)).mesh
        fixed = repair_orientation(m2)  # 已全正，应原样返回
        self.assertEqual(fixed.tets, m2.tets)
        self.assertEqual(fixed.materials, m2.materials)


class TestSplitInputValidation(unittest.TestCase):
    """参数形状/类型/范围非法 → MeshInputError（可定位）。"""

    def setUp(self):
        self.mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)

    def test_parameter_out_of_open_interval(self):
        for bad in (0, 1, Fraction(-1, 2), Fraction(3, 2), "0", "5/2"):
            with self.subTest(bad=bad):
                with self.assertRaises(MeshInputError) as cm:
                    split_internal_edge(self.mesh, (0, 6), bad)
                self.assertEqual(cm.exception.location, "t")

    def test_parameter_rejects_float_bool_nan(self):
        for bad in (0.5, True, "nan", "inf", "1/0"):
            with self.subTest(bad=bad):
                with self.assertRaises(MeshInputError) as cm:
                    split_internal_edge(self.mesh, (0, 6), bad)
                self.assertEqual(cm.exception.location, "t")

    def test_edge_shape_and_index_errors(self):
        cases = (
            ((0,), "edge"),
            ((0, 1, 2), "edge"),
            ((0, 0), "edge"),
            ((0, True), "edge[1]"),
            ((0, 99), "edge[1]"),
            ((Fraction(1, 2), 6), "edge[0]"),
            ("06", "edge"),
        )
        for edge, loc in cases:
            with self.subTest(edge=edge):
                with self.assertRaises(MeshInputError) as cm:
                    split_internal_edge(self.mesh, edge, Fraction(1, 2))
                self.assertEqual(cm.exception.location, loc)

    def test_mesh_type_checked(self):
        with self.assertRaises(TypeError):
            split_internal_edge("not a mesh", (0, 6), Fraction(1, 2))


class TestSplitRejections(unittest.TestCase):
    """拓扑/保护/退化拒绝：EdgeSplitError，且原网格状态不变（原子）。"""

    def assert_mesh_unchanged(self, mesh, snapshot):
        self.assertEqual(
            (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces),
            snapshot,
        )

    def test_edge_not_found(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (1, 7), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "edge_not_found")
        self.assert_mesh_unchanged(mesh, snap)

    def test_boundary_edge_rejected(self):
        mesh = TetMesh.build(CUBE_VERTICES, KUHN_TETS)
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (0, 1), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "edge_on_boundary")
        self.assertTrue(cm.exception.faces)
        self.assert_mesh_unchanged(mesh, snap)

    def test_protected_face_rejected(self):
        mesh = TetMesh.build(
            CUBE_VERTICES, KUHN_TETS, protected_faces=[(0, 1, 6)]
        )
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "edge_on_protected_face")
        self.assertEqual(cm.exception.faces, ((0, 1, 6),))
        self.assert_mesh_unchanged(mesh, snap)

    def test_non_manifold_edge_rejected(self):
        # Kuhn 环外加两个单元，使侧面 (0,2,6)、(0,3,6) 被 3 个单元引用，
        # 同时新侧面 (0,6,8) 仍保持重数 2（不引入边界侧面）。
        verts = (*CUBE_VERTICES, (2, 0, 0))
        tets = [*KUHN_TETS, (0, 2, 6, 8), (0, 6, 3, 8)]
        mesh = TetMesh.build(verts, tets)
        self.assertEqual(len(audit(mesh).non_manifold_faces), 2)
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "non_manifold_edge")
        self.assertEqual(cm.exception.faces, ((0, 2, 6), (0, 3, 6)))
        self.assert_mesh_unchanged(mesh, snap)

    def test_disconnected_ring_rejected(self):
        # Kuhn 6 环之外，再绕同一边 (0,6) 挂一个独立的 3 环：
        # 所有侧面重数都是 2，但邻接图是两个不相交的环。
        verts = (*CUBE_VERTICES, (2, 0, 0), (2, 2, 0), (1, 2, 0))
        tets = [*KUHN_TETS, (0, 6, 8, 9), (0, 6, 9, 10), (0, 6, 10, 8)]
        mesh = TetMesh.build(verts, tets)
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "ring_disconnected")
        self.assert_mesh_unchanged(mesh, snap)

    def test_zero_volume_star_rejected(self):
        # 绕边 (0,1) 的 3 环，其中 T2 四点共面（都在平面 y=z 上，手算 det=0）。
        verts = ((0, 0, 0), (1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1))
        tets = ((0, 1, 2, 3), (0, 1, 3, 4), (0, 1, 2, 4))
        mesh = TetMesh.build(verts, tets)
        self.assertEqual(audit(mesh).zero_tets, (2,))
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (0, 1), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "zero_volume_star")
        self.assert_mesh_unchanged(mesh, snap)

    def test_insertion_point_collides_with_existing_vertex(self):
        # 孤立顶点 8 恰好在对角线中点；t=1/2 的插入点与其重合。
        verts = (*CUBE_VERTICES, (Fraction(1, 2),) * 3)
        mesh = TetMesh.build(verts, KUHN_TETS)
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (0, 6), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "point_collides")
        self.assert_mesh_unchanged(mesh, snap)

    def test_too_many_tets_rejected(self):
        mesh = _ring_mesh(251)  # 251 个单元，分裂后 502 > 500
        snap = (mesh.vertices, mesh.tets, mesh.materials, mesh.protected_faces)
        with self.assertRaises(EdgeSplitError) as cm:
            split_internal_edge(mesh, (0, 1), Fraction(1, 2))
        self.assertEqual(cm.exception.reason, "too_many_tets")
        self.assert_mesh_unchanged(mesh, snap)

    def test_split_up_to_exactly_max_tets_accepted(self):
        mesh = _ring_mesh(250)  # 分裂后恰 500 = MAX_TETS
        res = split_internal_edge(mesh, (0, 1), Fraction(1, 2))
        rep = audit(res.mesh)
        self.assertEqual(rep.n_tets, 500)
        # 总 6 倍有向体积手算闭式：sum_{j=1}^{249} j(j+1) + (250-250^2)
        # = 249*250*251/3 - 62250 = 5146000。
        self.assertEqual(rep.total_signed_volume6, Fraction(5146000))
        self.assertEqual(rep.total_signed_volume, Fraction(2573000, 3))


class TestSplitEnumeration(unittest.TestCase):
    """确定性小规模枚举：随机顶点置换（覆盖槽位组合）× 有理 t。

    全部期望值由独立 Bareiss 参考（tests/reference.py）给出，不调用
    被测核心的行列式。
    """

    def test_random_permutations_and_parameters(self):
        rng = random.Random(20260924)
        perms = list(itertools.permutations(range(4)))
        ts = (
            Fraction(1, 2), Fraction(1, 3), Fraction(2, 5),
            Fraction(3, 7), Fraction(5, 11), Fraction(7, 13),
        )
        seen_slot_pairs = set()
        for _ in range(120):
            chosen = [rng.choice(perms) for _ in KUHN_TETS]
            tets = [tuple(t[p] for p in chosen[ti]) for ti, t in enumerate(KUHN_TETS)]
            t = rng.choice(ts)
            mesh = TetMesh.build(CUBE_VERTICES, tets)
            res = split_internal_edge(mesh, (0, 6), t)
            m2 = res.mesh
            for ti, (ia, ib) in enumerate(res.child_tet_indices):
                seen_slot_pairs.add(
                    (tets[ti].index(0), tets[ti].index(6))
                )
                pv = _vol6(CUBE_VERTICES, tets[ti])
                va = _vol6(m2.vertices, m2.tets[ia])
                vb = _vol6(m2.vertices, m2.tets[ib])
                self.assertEqual(va, t * pv)
                self.assertEqual(vb, (1 - t) * pv)
                self.assertNotEqual(va, 0)
                self.assertEqual(va > 0, pv > 0)
            old_total = sum(_vol6(CUBE_VERTICES, q) for q in tets)
            new_total = sum(_vol6(m2.vertices, q) for q in m2.tets)
            self.assertEqual(old_total, new_total)
            # 边界有序面多重集不变。
            old_b = sorted(
                o.ordered
                for f in build_face_incidence(mesh).values()
                if f.multiplicity == 1
                for o in (f.occurrences[0],)
            )
            new_b = sorted(
                o.ordered
                for f in build_face_incidence(m2).values()
                if f.multiplicity == 1
                for o in (f.occurrences[0],)
            )
            self.assertEqual(old_b, new_b)
            # 冲突/一致状态继承。
            old_inc = build_face_incidence(mesh)
            new_inc = build_face_incidence(m2)
            c = res.new_vertex_index
            for f in old_inc.values():
                if f.multiplicity == 2 and 0 in f.key and 6 in f.key:
                    x = next(v for v in f.key if v not in (0, 6))
                    for half in ((0, c, x), (6, c, x)):
                        hk = tuple(sorted(half))
                        self.assertEqual(
                            new_inc[hk].is_manifold_interior,
                            f.is_manifold_interior,
                        )
            # 全正父网格（偶置换）时结果必须整体一致定向。
            if all(_vol6(CUBE_VERTICES, q) > 0 for q in tets):
                self.assertTrue(audit(m2).consistently_oriented)
        # 枚举确实覆盖了父单元中 (a, b) 槽位的全部 6 种有序组合。
        self.assertEqual(
            seen_slot_pairs, set(itertools.permutations(range(4), 2))
        )


class TestMaterialsAndProtectedValidation(unittest.TestCase):
    """build 新选项的校验：材料标签与保护面。"""

    def test_materials_length_mismatch(self):
        with self.assertRaises(MeshInputError) as cm:
            TetMesh.build(CUBE_VERTICES, KUHN_TETS, materials=(1, 2))
        self.assertEqual(cm.exception.location, "materials")

    def test_materials_rejects_bool_and_float(self):
        for bad in (True, 1.5, "x"):
            with self.subTest(bad=bad):
                with self.assertRaises(MeshInputError) as cm:
                    TetMesh.build(
                        CUBE_VERTICES, KUHN_TETS,
                        materials=(0, 0, 0, 0, 0, bad),
                    )
                self.assertEqual(cm.exception.location, "materials[5]")

    def test_protected_face_validation(self):
        cases = (
            ([(0, 1)], "protected_faces[0]"),
            ([(0, 1, 2, 3)], "protected_faces[0]"),
            ([(0, 1, 1)], "protected_faces[0]"),
            ([(0, 1, True)], "protected_faces[0][2]"),
            ([(0, 1, 99)], "protected_faces[0][2]"),
        )
        for faces, loc in cases:
            with self.subTest(faces=faces):
                with self.assertRaises(MeshInputError) as cm:
                    TetMesh.build(
                        CUBE_VERTICES, KUHN_TETS, protected_faces=faces
                    )
                self.assertEqual(cm.exception.location, loc)

    def test_protected_faces_deduplicated_and_sorted(self):
        mesh = TetMesh.build(
            CUBE_VERTICES, KUHN_TETS,
            protected_faces=[(2, 1, 0), (0, 1, 2), (0, 6, 7)],
        )
        self.assertEqual(mesh.protected_faces, ((0, 1, 2), (0, 6, 7)))

    def test_repair_preserves_materials_and_protected(self):
        mats = (7, 6, 5, 4, 3, 2)
        tets = [list(t) for t in KUHN_TETS]
        tets[0] = [tets[0][0], tets[0][1], tets[0][3], tets[0][2]]  # 翻负
        mesh = TetMesh.build(
            CUBE_VERTICES, tets, materials=mats, protected_faces=[(0, 1, 2)]
        )
        fixed = repair_orientation(mesh)
        self.assertEqual(fixed.materials, mats)
        self.assertEqual(fixed.protected_faces, ((0, 1, 2),))
        self.assertEqual(fixed.tets[0], KUHN_TETS[0])
        self.assertEqual(audit(fixed).negative_tets, ())


if __name__ == "__main__":
    unittest.main()
