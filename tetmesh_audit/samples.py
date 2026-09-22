"""仓库内确定性示例网格（无外部数据、无随机源）。

单元立方体的 6 四面体扇形剖分：以主对角线 0→6 为轴，每个四面体
6 倍有向体积均为 1，总体积解析值为 1。

顶点编号（单位立方体）::

    0:(0,0,0)  1:(1,0,0)  2:(1,1,0)  3:(0,1,0)
    4:(0,0,1)  5:(1,0,1)  6:(1,1,1)  7:(0,1,1)
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

from .mesh import TetMesh

CUBE_VERTICES: List[Tuple[int, int, int]] = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
]

# 每个单元 6 倍体积 = +1，共享面诱导方向相反。
CUBE_TETS: List[Tuple[int, int, int, int]] = [
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
]


def cube_tetmesh() -> TetMesh:
    """定向完全一致的单位立方体（6 四面体扇形剖分）。"""
    return TetMesh(CUBE_VERTICES, CUBE_TETS)


def cube_with_flipped_tet(flip_index: int = 2) -> TetMesh:
    """单位立方体，但第 ``flip_index`` 个单元被翻转（奇置换 1↔3 号位）。"""
    tets = [list(t) for t in CUBE_TETS]
    t = tets[flip_index]
    t[1], t[3] = t[3], t[1]
    return TetMesh(CUBE_VERTICES, tets)


def nonmanifold_mesh():
    """立方体底面三角形 (0,1,2) 被三个单元引用的非流形网格。

    在立方体下方追加两个不同顶点 8、9，各自与 (0,1,2) 组成单元，
    使该面关联单元数为 3（T1 与两个下方单元）。两个追加单元体积为负，
    可被修复规则翻转，但非流形拓扑不因此消失。
    """
    verts = list(CUBE_VERTICES) + [(0, 0, -1), (1, 1, -1)]
    tets = list(CUBE_TETS) + [(0, 1, 2, 8), (0, 1, 2, 9)]
    return TetMesh(verts, tets)


def coplanar_degenerate_mesh() -> TetMesh:
    """含一个零体积单元的网格：四点全部落在 z=0 平面上。"""
    verts = list(CUBE_VERTICES) + [(2, 0, 0), (2, 1, 0)]
    # (1,2,8,9) 全在 z=0 → 零体积；另附一个正常单元保证其余部分有效。
    tets = list(CUBE_TETS) + [(1, 2, 8, 9)]
    return TetMesh(verts, tets)


def overlapping_same_side_mesh() -> TetMesh:
    """两个正体积单元共享面且第四顶点在面的同一侧（诱导方向相同）。

    面 (0,1,2) 关联恰两个单元，但两单元互相“穿刺”重叠：内面诱导方向
    不相反，构成 inconsistent_interior_face 缺陷（修复规则无法处理）。
    """
    verts = list(CUBE_VERTICES) + [(0, 0, 2)]
    tets = [(0, 1, 2, 6), (0, 1, 2, 8)]
    return TetMesh(verts, tets)


def scaled_cube_tetmesh(scale) -> TetMesh:
    """坐标为分数的缩放立方体，用于验证有理输入与精确体积 s³。"""
    from fractions import Fraction

    s = Fraction(scale)
    verts = [(s * x, s * y, s * z) for (x, y, z) in CUBE_VERTICES]
    return TetMesh(verts, CUBE_TETS)


def two_disjoint_cubes() -> TetMesh:
    """两个互不接触的单位立方体，用于连通分量测试（总体积 2）。"""
    verts = list(CUBE_VERTICES) + [
        (x + 3, y, z) for (x, y, z) in CUBE_VERTICES
    ]
    tets = list(CUBE_TETS) + [(a + 8, b + 8, c + 8, d + 8) for (a, b, c, d) in CUBE_TETS]
    return TetMesh(verts, tets)
