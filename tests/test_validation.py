"""输入校验：NaN/Infinity、bool、越界、重复单元、上限等错误语义。"""

import unittest
from fractions import Fraction

from tetmesh_audit import (
    DegenerateTetrahedronError,
    DuplicateTetrahedronError,
    InvalidMeshError,
    MAX_TETRAHEDRA,
    TetMesh,
)

CUBE_VERTS = [
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
]
CUBE_TETS = [
    (0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
    (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6),
]


class ValidationTests(unittest.TestCase):
    def test_nan_coordinate_rejected_with_location(self):
        verts = [list(v) for v in CUBE_VERTS]
        verts[2][1] = float("nan")
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(verts, CUBE_TETS)
        self.assertIn("vertices[2][1]", str(cm.exception))
        self.assertEqual(cm.exception.path, "vertices[2][1]")

    def test_infinity_coordinates_rejected(self):
        for bad in (float("inf"), float("-inf")):
            verts = [list(v) for v in CUBE_VERTS]
            verts[0][2] = bad
            with self.assertRaises(InvalidMeshError) as cm:
                TetMesh(verts, CUBE_TETS)
            self.assertEqual(cm.exception.path, "vertices[0][2]")

    def test_bool_as_vertex_index_rejected(self):
        tets = [list(t) for t in CUBE_TETS]
        tets[1][2] = True
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(CUBE_VERTS, tets)
        self.assertEqual(cm.exception.path, "tetrahedra[1][2]")

    def test_float_index_rejected(self):
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(CUBE_VERTS, [(0, 1.0, 2, 6)])
        self.assertEqual(cm.exception.path, "tetrahedra[0][1]")

    def test_bool_coordinate_rejected(self):
        verts = [list(v) for v in CUBE_VERTS]
        verts[3][0] = True
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(verts, CUBE_TETS)
        self.assertEqual(cm.exception.path, "vertices[3][0]")

    def test_out_of_range_and_negative_index(self):
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(CUBE_VERTS, [(0, 1, 2, 8)])
        self.assertEqual(cm.exception.path, "tetrahedra[0][3]")
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(CUBE_VERTS, [(0, -1, 2, 6)])
        self.assertEqual(cm.exception.path, "tetrahedra[0][1]")

    def test_repeated_vertices_in_tet_rejected(self):
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(CUBE_VERTS, [(0, 0, 2, 6)])
        self.assertEqual(cm.exception.path, "tetrahedra[0]")

    def test_wrong_lengths(self):
        with self.assertRaises(InvalidMeshError):
            TetMesh([(0, 0)], CUBE_TETS)
        with self.assertRaises(InvalidMeshError):
            TetMesh(CUBE_VERTS, [(0, 1, 2)])

    def test_empty_inputs_rejected(self):
        with self.assertRaises(InvalidMeshError):
            TetMesh([], CUBE_TETS)
        with self.assertRaises(InvalidMeshError):
            TetMesh(CUBE_VERTS, [])

    def test_non_numeric_coordinate_rejected(self):
        verts = [list(v) for v in CUBE_VERTS]
        verts[0] = ("x", 0, 0)
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(verts, CUBE_TETS)
        self.assertEqual(cm.exception.path, "vertices[0][0]")

    def test_decimal_string_and_fraction_accepted_exactly(self):
        verts = [(Fraction(0), "0.5", "1/3"), (1, 0, 0), (1, 1, 0), (0, 1, 1)]
        mesh = TetMesh(verts, [(0, 1, 2, 3)])
        self.assertEqual(mesh.vertices[0][1], Fraction(1, 2))
        self.assertEqual(mesh.vertices[0][2], Fraction(1, 3))

    def test_generators_accepted(self):
        mesh = TetMesh((v for v in CUBE_VERTS), (t for t in CUBE_TETS))
        self.assertEqual(mesh.num_tetrahedra, 6)

    def test_limit_500_accepted_501_rejected(self):
        verts = []
        tets = []
        for k in range(MAX_TETRAHEDRA):
            base = len(verts)
            verts += [
                (10 * k, 0, 0), (10 * k + 1, 0, 0),
                (10 * k, 1, 0), (10 * k, 0, 1),
            ]
            tets.append((base, base + 1, base + 2, base + 3))
        mesh = TetMesh(verts, tets)
        self.assertEqual(mesh.num_tetrahedra, 500)

        verts.append((0, 100, 0))
        tets.append((0, 1, 2, len(verts) - 1))
        with self.assertRaises(InvalidMeshError) as cm:
            TetMesh(verts, tets)
        self.assertIn("500", str(cm.exception))
        self.assertEqual(cm.exception.path, "tetrahedra[500]")


class DuplicateTests(unittest.TestCase):
    def test_identical_duplicate_reported_separately(self):
        tets = CUBE_TETS + [CUBE_TETS[0]]
        with self.assertRaises(DuplicateTetrahedronError) as cm:
            TetMesh(CUBE_VERTS, tets)
        self.assertEqual(cm.exception.duplicates[0][:2], (0, 6))

    def test_permuted_duplicate_is_still_duplicate(self):
        # 同一四点集合，奇排列（定向相反）也算重复单元。
        tets = CUBE_TETS + [(6, 2, 1, 0)]
        with self.assertRaises(DuplicateTetrahedronError) as cm:
            TetMesh(CUBE_VERTS, tets)
        self.assertEqual(cm.exception.duplicates[0][2], (0, 1, 2, 6))

    def test_degenerate_error_type_distinct(self):
        verts = CUBE_VERTS + [(2, 0, 0), (2, 1, 0)]
        tets = CUBE_TETS + [(1, 2, 8, 9)]
        mesh = TetMesh(verts, tets)  # 构造期不报体积问题
        from tetmesh_audit import audit, repair_orientation

        result = audit(mesh)
        self.assertEqual(result.zero_volume_tets, (6,))
        with self.assertRaises(DegenerateTetrahedronError) as cm:
            repair_orientation(mesh)
        self.assertEqual(cm.exception.zero_volume_tets, [6])


if __name__ == "__main__":
    unittest.main()
