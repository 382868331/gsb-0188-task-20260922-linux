"""测试专用独立参考实现。

约束：本模块不得 import tetmesh_audit，所有期望值均由独立推导的公式
（Leibniz 排列展开、散度定理通量公式、解析体积）直接算出。
"""

from fractions import Fraction
from itertools import permutations


def _mat3_det(a, b, c):
    """3x3 行列式，Leibniz 六项带符号排列展开（与核心的余子式展开独立）。"""
    cols = (a, b, c)
    total = 0
    for perm in permutations(range(3)):
        # 排列逆序数奇偶决定 Leibniz 符号。
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        sign = 1 if inv % 2 == 0 else -1
        total += sign * cols[0][perm[0]] * cols[1][perm[1]] * cols[2][perm[2]]
    return total


def tet_signed_six_volume(p0, p1, p2, p3):
    u = tuple(p1[k] - p0[k] for k in range(3))
    v = tuple(p2[k] - p0[k] for k in range(3))
    w = tuple(p3[k] - p0[k] for k in range(3))
    return _mat3_det(u, v, w)


def tet_signed_volume(p0, p1, p2, p3):
    return Fraction(tet_signed_six_volume(p0, p1, p2, p3), 6)


def triangle_det_about_origin(p, q, r):
    """det(p,q,r)：原点到有向三角面的 6 倍有向通量贡献，独立计算。"""
    return _mat3_det(tuple(p), tuple(q), tuple(r))


def flux_volume(vertices, oriented_triangles):
    """散度定理独立核对：对有向闭合三角面网求 1/6 * Σ det(p,q,r)。

    由 F=(x,y,z)/3 推出体积 V = 1/6 Σ det(p,q,r)；此处逐面用参考
    行列式实现求和，不经过任何被测函数。
    """
    total = 0
    for tri in oriented_triangles:
        p, q, r = (vertices[i] for i in tri)
        total += triangle_det_about_origin(p, q, r)
    return Fraction(total, 6)


def face_incidence(tetrahedra):
    """从原始单元表独立统计无序三点面的关联（单元下标, 删除位）集合。"""
    removal = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    inc = {}
    for ti, t in enumerate(tetrahedra):
        for lf, rem in enumerate(removal):
            key = tuple(sorted((t[rem[0]], t[rem[1]], t[rem[2]])))
            inc.setdefault(key, []).append((ti, lf))
    return inc


def box_volume_analytic(side):
    """轴对齐立方体解析体积 s^3。"""
    s = Fraction(side)
    return s * s * s


def triangle_normal_sign(vertices, tri, axis, direction):
    """独立判断三角面外法向沿某坐标轴的符号（+1/-1/0）。

    用于边界定向抽查：例如立方体 z=1 面上的面外法向应为 +z。
    """
    p, q, r = (vertices[i] for i in tri)
    u = tuple(q[k] - p[k] for k in range(3))
    v = tuple(r[k] - p[k] for k in range(3))
    # 叉积分量（独立写出）
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    n = cross[axis]
    return (1 if n > 0 else -1 if n < 0 else 0) * (1 if direction > 0 else -1)
