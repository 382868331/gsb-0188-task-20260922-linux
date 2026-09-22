"""精确有理几何原语。"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence, Tuple

# 单元 (a,b,c,d) 中，与第 i 个顶点相对的面；三元组按“删除第 i 个顶点后
# 剩余顶点在原四元组中的顺序”排列（单纯形边界算子的原始顺序）。
FACE_REMOVAL_ORDER: Tuple[Tuple[int, int, int], ...] = (
    (1, 2, 3),  # 面 0，与顶点 0 相对
    (0, 2, 3),  # 面 1，与顶点 1 相对
    (0, 1, 3),  # 面 2，与顶点 2 相对
    (0, 1, 2),  # 面 3，与顶点 3 相对
)


def signed_six_volume(
    a: Sequence[Fraction],
    b: Sequence[Fraction],
    c: Sequence[Fraction],
    d: Sequence[Fraction],
):
    """6 倍有向体积 det(b-a, c-a, d-a)，精确整数或 Fraction。

    正值表示 (a,b,c,d) 为右手定向；零表示四点共面；负值表示定向翻转。
    """
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    wx, wy, wz = d[0] - a[0], d[1] - a[1], d[2] - a[2]
    return (
        ux * (vy * wz - vz * wy)
        - uy * (vx * wz - vz * wx)
        + uz * (vx * wy - vy * wx)
    )


def chain_triangle(tet: Sequence[int], local_face: int) -> Tuple[int, int, int]:
    """单元在其面（local_face）上的**边界链诱导方向**：(-1)^i · 删除序。

    这是纯粹的单纯形边界算子，与单元体积符号无关：两个一致定向的
    相邻单元在共享面上的诱导方向必然相反。
    """
    rem = FACE_REMOVAL_ORDER[local_face]
    triple = (tet[rem[0]], tet[rem[1]], tet[rem[2]])
    if local_face % 2 == 1:
        return (triple[2], triple[1], triple[0])
    return triple


def outward_triangle(tet: Sequence[int], local_face: int, six_volume) -> Tuple[int, int, int]:
    """返回该面**几何外法向**定向的顶点三元组（边界输出用）。

    正体积单元的边界链方向即朝外；负体积单元整体反向，需再取反；
    零体积时方向无定义，按删除顺序原样返回（调用方以
    ``orientation_valid=False`` 另行标注）。
    """
    if six_volume == 0:
        rem = FACE_REMOVAL_ORDER[local_face]
        return (tet[rem[0]], tet[rem[1]], tet[rem[2]])
    triple = chain_triangle(tet, local_face)
    if six_volume < 0:
        return (triple[2], triple[1], triple[0])
    return triple


def _permutation_parity(triple: Tuple[int, int, int]) -> int:
    """三元组相对其升序排列的逆序数奇偶（0 偶 / 1 奇）。"""
    inv = 0
    for i in range(3):
        for j in range(i + 1, 3):
            if triple[i] > triple[j]:
                inv += 1
    return inv % 2


def induced_orientations_opposite(
    key: Tuple[int, int, int], tr0: Tuple[int, int, int], tr1: Tuple[int, int, int]
) -> bool:
    """两个三点诱导方向是否相反。

    循环重排（如 (6,3,0) 与 (3,6,0)）表示同一方向，不能按元组字面值
    比较；正确判据是两者相对无序面键 ``key`` 的排列奇偶性不同。
    """
    rank0 = tuple(key.index(v) for v in tr0)
    rank1 = tuple(key.index(v) for v in tr1)
    return _permutation_parity(rank0) != _permutation_parity(rank1)
