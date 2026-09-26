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


# ---------------------------------------------------------------------------
# 线段穿越的独立参考（第二轮新增）。
# 与 orient_audit.traverse 完全独立：不使用面邻接行走，也不调用库中任何
# 函数，而是对每个四面体独立解半空间不等式组裁出参数区间。

Vec = tuple[Fraction, Fraction, Fraction]


def _sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(u: Vec, v: Vec) -> Vec:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def _dot(u: Vec, v: Vec) -> Fraction:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def ref_clip_intervals(verts, tets, start, end):
    """逐四面体独立裁线：返回 [(tet_index, lo, hi), ...]（裁剪到 [0,1]）。

    方法：对每个面用叉积法向建立半空间不等式（内侧符号由对顶点确定），
    沿线段参数 t 解线性不等式组。与库中行走实现不共享任何代码路径。
    """

    a = _as_vec(start)
    d = _sub(_as_vec(end), a)
    out = []
    for ti, tet in enumerate(tets):
        p = [_as_vec(verts[k]) for k in tet]
        lo, hi = Fraction(0), Fraction(1)
        for i in range(4):
            q = [p[j] for j in range(4) if j != i]
            n = _cross(_sub(q[1], q[0]), _sub(q[2], q[0]))
            # 内侧方向：对顶点一侧。
            side = _dot(n, _sub(p[i], q[0]))
            if side == 0:
                continue  # 退化四面体；参考只用于合法网格
            # 不等式：side * (n.(X - q0)) >= 0，X = a + t d
            c0 = side * _dot(n, _sub(a, q[0]))
            c1 = side * _dot(n, d)
            if c1 == 0:
                if c0 < 0:
                    lo, hi = Fraction(1), Fraction(0)
                    break
                continue
            root = -c0 / c1
            if c1 > 0:
                lo = max(lo, root)
            else:
                hi = min(hi, root)
        if lo < hi:
            out.append((ti, lo, hi))
    return out


def ref_merge_intervals(intervals):
    """把 [(lo, hi), ...] 合并为互不重叠的区间并集（排序后扫描）。"""

    merged = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1]:
            if hi > merged[-1][1]:
                merged[-1] = (merged[-1][0], hi)
        else:
            merged.append((lo, hi))
    return merged


def ref_affine_integral(alpha, beta, start, end, intervals):
    """仿射场 f(X) = alpha.X + beta 在参数区间集上的解析积分。

    X(t) = a + t d，积分为 (alpha.a + beta)(t1-t0) + (alpha.d)(t1^2-t0^2)/2。
    不经过任何单元/重心计算，是独立的解析参考。
    """

    a = _as_vec(start)
    d = _sub(_as_vec(end), a)
    base = _dot(alpha, a) + beta
    slope = _dot(alpha, d)
    total = Fraction(0)
    for lo, hi in intervals:
        total += base * (hi - lo) + slope * (hi * hi - lo * lo) / 2
    return total
