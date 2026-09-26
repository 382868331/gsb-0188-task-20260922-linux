"""守恒重映射真实离线演示（约 1 秒，远小于 90 秒预算）。

全部数据由本文件确定性构造，真实调用 orient_audit.remap 计算，
不依赖账号、密钥、网络、sleep 或预录结果。

演示三个场景：
1. 同一立方体的两种四面体剖分之间转常量场（手算可核验）；
2. 跨分界面 x=1/2 的分段常量场，源网格贴合、目标不贴合；
3. 目标区域大于源区域：缺口保留，不归一。
"""

from __future__ import annotations

import itertools
import time
from fractions import Fraction as F

from orient_audit import TetMesh, remap


def idx(x, y, z):
    return x * 4 + y * 2 + z


T6 = [
    (idx(0, 0, 0), idx(1, 0, 0), idx(1, 1, 0), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(1, 0, 0), idx(1, 0, 1), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 1, 0), idx(1, 1, 0), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 1, 0), idx(0, 1, 1), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 0, 1), idx(1, 0, 1), idx(1, 1, 1)),
    (idx(0, 0, 0), idx(0, 0, 1), idx(0, 1, 1), idx(1, 1, 1)),
]

T5 = [
    (idx(0, 0, 0), idx(1, 0, 0), idx(0, 1, 0), idx(0, 0, 1)),
    (idx(1, 1, 0), idx(1, 0, 0), idx(0, 1, 0), idx(1, 1, 1)),
    (idx(1, 0, 1), idx(1, 0, 0), idx(0, 0, 1), idx(1, 1, 1)),
    (idx(0, 1, 1), idx(0, 1, 0), idx(0, 0, 1), idx(1, 1, 1)),
    (idx(1, 0, 0), idx(0, 1, 0), idx(0, 0, 1), idx(1, 1, 1)),
]


def unit_cube_vertices():
    return [
        (F(x), F(y), F(z)) for x in (0, 1) for y in (0, 1) for z in (0, 1)
    ]


def grid_mesh(n, scale):
    """[0, n*scale]^3 的规则网格，每个小立方体 6 四面体剖分。"""
    verts = {}
    tets = []

    def vid(i, j, k):
        key = (i, j, k)
        if key not in verts:
            verts[key] = len(verts)
        return verts[key]

    t6_offsets = [
        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)),
        ((0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1)),
        ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)),
        ((0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)),
        ((0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)),
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)),
    ]
    for c in itertools.product(range(n), repeat=3):
        for t in t6_offsets:
            tets.append(tuple(vid(c[0] + d[0], c[1] + d[1], c[2] + d[2]) for d in t))
    vlist = [None] * len(verts)
    for (i, j, k), v in verts.items():
        vlist[v] = (F(i) * scale, F(j) * scale, F(k) * scale)
    return TetMesh.build(vlist, tets)


def show(name, rep, n_target):
    print(f"  稀疏权重非零项: {len(rep.weights)}")
    print(f"  源总质量:       {rep.source_mass}")
    print(f"  源被覆盖质量:   {rep.covered_source_mass}")
    print(f"  目标总质量:     {rep.target_mass}")
    print(f"  总量误差:       {rep.mass_error}")
    print(f"  未覆盖体积合计: {sum(rep.uncovered_volumes)} / 目标总体积 {sum(rep.target_volumes)}")


def main():
    t0 = time.perf_counter()

    print("场景 1：同一立方体，6 四面体剖分(源) -> 5 四面体剖分(目标)，常量场 7/3")
    src = TetMesh.build(unit_cube_vertices(), T6)
    tgt = TetMesh.build(unit_cube_vertices(), T5)
    rep = remap(src, [F(7, 3)] * 6, tgt)
    show("常量场", rep, 5)
    print(f"  目标单元平均密度: {list(rep.target_averages)}（手算：全为 7/3）")
    assert rep.mass_error == 0 and all(a == F(7, 3) for a in rep.target_averages)

    print()
    print("场景 2：[0,2]^3 上 2x2x2 源网格 -> 3x3x3 目标网格，分段常量场")
    src2 = grid_mesh(2, F(1))
    tgt2 = grid_mesh(3, F(2, 3))
    # 每个源小立方体的 6 个四面体密度相同；x<1 为 1，x>=1 为 3（分界面不贴合目标）
    rho = []
    for tet in src2.tets:
        cx = sum(src2.vertices[i][0] for i in tet) / 4
        rho.append(F(1) if cx < 1 else F(3))
    rep2 = remap(src2, rho, tgt2)
    show("分段常量", rep2, len(tgt2.tets))
    print(f"  手算总质量: {1 * 4 + 3 * 4}（左半 4 立方体密度 1，右半密度 3）")
    assert rep2.mass_error == 0 and rep2.target_mass == 16

    print()
    print("场景 3：目标 [0,3]^3 大于源 [0,2]^3，缺口保留不归一")
    tgt3 = grid_mesh(2, F(3, 2))
    rep3 = remap(src2, rho, tgt3)
    show("部分覆盖", rep3, len(tgt3.tets))
    print(f"  手算未覆盖体积: {3**3 - 2**3}（目标 27 减去源 8）")
    assert rep3.mass_error == 0
    assert sum(rep3.uncovered_volumes) == 19
    assert rep3.target_mass == rep3.covered_source_mass == 16

    dt = time.perf_counter() - t0
    print()
    print(f"全部场景守恒核验通过，真实计算耗时 {dt:.2f} 秒（预算 90 秒）")


if __name__ == "__main__":
    main()
