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
4. ``ref_clip_segment``：独立的逐四面体裁线（Cyrus–Beck 风格：
   对每个面用叉积外法向 + 点积半空间裁剪），与被测核心的
   "顶点替换行列式 f_i(t)" 是两条不同的代数路径；
5. ``ref_affine_integral``：仿射场在线段参数区间上的解析积分，
   不经过任何网格结构。
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


_FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (2, 1, 0))


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(u: Vec3, v: Vec3) -> Vec3:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def _dot(u: Vec3, v: Vec3) -> Fraction:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def ref_clip_segment(verts, tet, p0, p1):
    """独立裁线参考：返回线段在四面体内的参数区间 (lo, hi) 或 None。

    方法（与被测核心不同）：对每个诱导面取叉积法向，按"与相对顶点
    反向"翻成外法向，再用半空间 n·(P(t) - p) <= 0 做 Cyrus–Beck 裁剪。
    """

    pts = [_as_vec(verts[k]) for k in tet]
    a = _as_vec(p0)
    b = _as_vec(p1)
    d = _sub(b, a)
    lo, hi = Fraction(0), Fraction(1)
    for local in range(4):
        face = _FACES[local]
        p, q, r = (pts[face[0]], pts[face[1]], pts[face[2]])
        n = _cross(_sub(q, p), _sub(r, p))
        # 诱导朝向的法向朝外（正体积单元），直接用即可；
        # 为稳妥仍按相对顶点校正符号。
        opposite_idx = ({0, 1, 2, 3} - set(face)).pop()
        if _dot(n, _sub(pts[opposite_idx], p)) > 0:
            n = (-n[0], -n[1], -n[2])
        # 半空间：n·(a + t d - p) <= 0
        c0 = _dot(n, _sub(a, p))
        c1 = _dot(n, d)
        if c1 > 0:
            bound = -c0 / c1
            if bound < hi:
                hi = bound
        elif c1 < 0:
            bound = -c0 / c1
            if bound > lo:
                lo = bound
        elif c0 > 0:
            return None
        if lo > hi:
            return None
    return (lo, hi)


def ref_affine_integral(coeffs, p0, p1, t_a, t_b) -> Fraction:
    """仿射场 f(x,y,z) = c0*x + c1*y + c2*z + c3 在 P(t) 上的解析积分。

    ∫_{t_a}^{t_b} f(p0 + t*(p1-p0)) dt，直接对一次多项式求原函数。
    """

    c0, c1, c2, c3 = (Fraction(c) for c in coeffs)
    a = _as_vec(p0)
    b = _as_vec(p1)
    d = _sub(b, a)
    # f(t) = (c0*d0 + c1*d1 + c2*d2) * t + (c0*a0 + c1*a1 + c2*a2 + c3)
    slope = c0 * d[0] + c1 * d[1] + c2 * d[2]
    intercept = c0 * a[0] + c1 * a[1] + c2 * a[2] + c3
    ta, tb = Fraction(t_a), Fraction(t_b)
    return slope * (tb * tb - ta * ta) / 2 + intercept * (tb - ta)


def ref_affine_value(coeffs, point) -> Fraction:
    """仿射场在一点的精确值。"""

    c0, c1, c2, c3 = (Fraction(c) for c in coeffs)
    p = _as_vec(point)
    return c0 * p[0] + c1 * p[1] + c2 * p[2] + c3
