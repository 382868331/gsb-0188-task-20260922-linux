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


# ============================================================
# 平面剖切的独立参考（同样不 import orient_audit）
# ============================================================


def ref_dot(u, v) -> Fraction:
    return sum(Fraction(u[k]) * Fraction(v[k]) for k in range(3))


def ref_cross(u, v):
    u = tuple(Fraction(x) for x in u)
    v = tuple(Fraction(x) for x in v)
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def ref_on_plane(coord, normal, offset) -> bool:
    """点精确满足平面方程 n·x = d。"""

    return ref_dot(coord, normal) == Fraction(offset)


def _vol6(a, b, c, d) -> Fraction:
    """det(b-a, c-a, d-a)，直接展开（独立于被测核心的实现路径）。"""

    a, b, c, d = (_as_vec(p) for p in (a, b, c, d))
    u = tuple(b[k] - a[k] for k in range(3))
    v = tuple(c[k] - a[k] for k in range(3))
    w = tuple(d[k] - a[k] for k in range(3))
    mat = [
        [Fraction(1), *a],
        [Fraction(1), *b],
        [Fraction(1), *c],
        [Fraction(1), *d],
    ]
    # 4x4 行列式（与 ref_signed_volume6 相同的 Bareiss 路径，输入为任意点）。
    n = 4
    prev = Fraction(1)
    swap = 1
    for k in range(n - 1):
        if mat[k][k] == 0:
            for i in range(k + 1, n):
                if mat[i][k] != 0:
                    mat[k], mat[i] = mat[i], mat[k]
                    swap = -swap
                    break
            else:
                return Fraction(0)
        pivot = mat[k][k]
        bm = [row[:] for row in mat]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                bm[i][j] = (mat[i][j] * pivot - mat[i][k] * mat[k][j]) / prev
        mat = bm
        prev = pivot
    return swap * mat[n - 1][n - 1]


def ref_barycentric(coord, tet_points):
    """点相对四面体的 4 个重心坐标（精确）。

    返回 (λ0..λ3, 总体积6)；点在闭单体内 iff 全部 λ >= 0 且和为 1。
    实现独立于被测核心：用 4x4 Bareiss 体积比。
    """

    p = _as_vec(coord)
    a, b, c, d = (_as_vec(q) for q in tet_points)
    total = _vol6(a, b, c, d)
    la = _vol6(p, b, c, d) / total
    lb = _vol6(a, p, c, d) / total
    lc = _vol6(a, b, p, d) / total
    ld = _vol6(a, b, c, p) / total
    return (la, lb, lc, ld), total


def ref_polygon_alpha(cyclic_coords, normal) -> Fraction:
    """平面多边形面积系数 α 的独立算法：边界叉积和。

        2·(面积/|n|) = |n·Σ p_i × p_{i+1}| / (n·n)

    与被测核心的“逐三角形 (q-p)×(r-p) 扇形和”是不同的公式
    （边界和 vs 三角剖分和），二者恒等可作为独立核对。
    """

    n = tuple(Fraction(x) for x in normal)
    nn = ref_dot(n, n)
    total = [Fraction(0), Fraction(0), Fraction(0)]
    m = len(cyclic_coords)
    for i in range(m):
        cr = ref_cross(cyclic_coords[i], cyclic_coords[(i + 1) % m])
        for k in range(3):
            total[k] += cr[k]
    return abs(ref_dot(n, total)) / (2 * nn)


def ref_pyramid_alpha(base_cyclic, apex, normal, offset) -> Fraction:
    """解析棱锥核对面积系数：棱锥体积 V = 底面积·高/3。

    底为截面三角形，顶点为被切下的角顶点：
        面积 = 3V/h，  h = |n·apex − d| / |n|
        => α = 面积/|n| = |det6| / (2|n·apex − d|)
    det6 用 4x4 Bareiss 独立给出。仅适用于三角形底。
    """

    assert len(base_cyclic) == 3
    n = tuple(Fraction(x) for x in normal)
    p, q, r = (_as_vec(x) for x in base_cyclic)
    a = _as_vec(apex)
    det6 = abs(_vol6(a, p, q, r))
    height_num = abs(ref_dot(n, a) - Fraction(offset))
    return det6 / (2 * height_num)
