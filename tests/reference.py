"""测试专用的*独立*参考实现。

纪律：本文件绝不 import orient_audit，期望值由下列与被测核心不同的
算法独立给出：

1. ``ref_signed_volume6``：4x4 齐次坐标行列式，Bareiss 无分数高斯消元
   （被测核心用的是 3x3 余子式展开，算法路径不同）；
2. ``ref_flux_volume``：闭合定向三角面的散度/面通量公式
   ``V = (1/6) * sum det(p, q, r)``，直接对边界面求和，
   不经过任何单元体积；
3. ``ref_outward``：用面心 - 单元体形心与叉积法向做点积，独立判断
   诱导面朝向是否朝外。
"""

from __future__ import annotations

from fractions import Fraction

Vec3 = tuple[Fraction, Fraction, Fraction]


def _as_vec(p) -> Vec3:
    return (Fraction(p[0]), Fraction(p[1]), Fraction(p[2]))


def ref_signed_volume6(verts, tet) -> Fraction:
    """4x4 齐次行列式（行 [1, x, y, z]），Bareiss 精确消元。"""

    a = [
        [Fraction(1), *_as_vec(verts[tet[i]])]
        for i in range(4)
    ]
    n = 4
    prev = Fraction(1)
    swap_sign = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            for i in range(k + 1, n):
                if a[i][k] != 0:
                    a[k], a[i] = a[i], a[k]
                    swap_sign = -swap_sign
                    break
            else:
                return Fraction(0)
        pivot = a[k][k]
        b = [row[:] for row in a]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                b[i][j] = (a[i][j] * pivot - a[i][k] * a[k][j]) / prev
        a = b
        prev = pivot
    return swap_sign * a[n - 1][n - 1]


def _det3(p: Vec3, q: Vec3, r: Vec3) -> Fraction:
    return (
        p[0] * (q[1] * r[2] - q[2] * r[1])
        - p[1] * (q[0] * r[2] - q[2] * r[0])
        + p[2] * (q[0] * r[1] - q[1] * r[0])
    )


def ref_flux_volume(verts, boundary_triples) -> Fraction:
    """面通量公式：对闭合定向三角网格 V = (1/6) sum det(p,q,r)。"""

    total = Fraction(0)
    for i, j, k in boundary_triples:
        total += _det3(_as_vec(verts[i]), _as_vec(verts[j]), _as_vec(verts[k]))
    return total / 6


def ref_outward(verts, tet, ordered_face) -> bool:
    """面诱导朝向是否朝外：叉积法向·(面心-体形心) > 0（精确）。"""

    pts = [_as_vec(verts[t]) for t in tet]
    p, q, r = (_as_vec(verts[t]) for t in ordered_face)
    u = (q[0] - p[0], q[1] - p[1], q[2] - p[2])
    v = (r[0] - p[0], r[1] - p[1], r[2] - p[2])
    normal = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    face_c = tuple((p[k] + q[k] + r[k]) / 3 for k in range(3))
    tet_c = tuple(sum(pt[k] for pt in pts) / 4 for k in range(3))
    outward = tuple(face_c[k] - tet_c[k] for k in range(3))
    return sum(normal[k] * outward[k] for k in range(3)) > 0


def ref_barycentric(verts, tet, point):
    """点相对四面体的重心坐标，Cramer 法则精确求解（独立于被测核心）。

    返回 4 个 Fraction；退化单元返回 None。坐标全部落在 [0,1]
    即表示点（含边界）在单元内。
    """

    p = [_as_vec(verts[t]) for t in tet]
    x = _as_vec(point)
    cols = [tuple(p[i][k] - p[0][k] for k in range(3)) for i in (1, 2, 3)]
    rhs = tuple(x[k] - p[0][k] for k in range(3))
    d = _det3(cols[0], cols[1], cols[2])
    if d == 0:
        return None
    l1 = _det3(rhs, cols[1], cols[2]) / d
    l2 = _det3(cols[0], rhs, cols[2]) / d
    l3 = _det3(cols[0], cols[1], rhs) / d
    return (1 - l1 - l2 - l3, l1, l2, l3)
