#!/usr/bin/env python3
"""守恒重映射演示：两张四面体网格之间的单元平均量转移。

真实离线计算，无 sleep、无预录结果、无网络。三部分：

A. 同一单位立方体的两种四面体剖分（5 剖分 -> 绕对角线 6 剖分），
   常量密度场：每个目标单元的平均密度必须精确等于该常量，
   且总质量精确守恒（mass_error == 0）。

B. 跨界分段常量场：源为 [0,1/2]（密度 1）与 [1/2,1]（密度 3）
   两个半立方体，目标为中间条 [1/4,3/4]x[0,1]x[0,1]。
   总质量手算 = 1/4*1 + 1/4*3 = 1，与程序精确值对照。

C. 目标网格超出源区域（[0,1]x[0,1]x[0,2]）：上半无覆盖，
   缺口保留在 target_uncovered_volume 中，不做满覆盖归一。

运行：python demo_iteration.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import TetMesh, remap

# 单位立方体顶点（与 tests/fixtures.py 相同编号）。
CUBE_VERTICES = (
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
)
# 标准 5 四面体剖分：中心正四面体 1/3 + 四角各 1/6。
CUBE5_TETS = (
    (0, 1, 2, 5), (0, 2, 3, 7), (0, 4, 5, 7), (2, 5, 6, 7), (0, 2, 7, 5),
)
# 绕主对角线 (0,0,0)-(1,1,1) 的 6 四面体剖分，每单元 1/6。
CUBE6_TETS = (
    (0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
    (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6),
)


def box_vertices(x0, x1, y0=0, y1=1, z0=0, z1=1):
    return (
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    )


def main() -> None:
    start = time.perf_counter()

    print("=" * 68)
    print("A. 同一立方体：5 剖分(源) -> 6 剖分(目标)，常量密度 7/2")
    print("=" * 68)
    source = TetMesh.build(CUBE_VERTICES, CUBE5_TETS)
    target = TetMesh.build(CUBE_VERTICES, CUBE6_TETS)
    result = remap(source, target, [Fraction(7, 2)] * 5)
    print(f"稀疏权重条数: {len(result.weights)}（仅交叠体积>0 的源/目标对）")
    for w in result.weights:
        print(f"  源[{w.source_index}] -> 目标[{w.target_index}]  交叠体积 {w.volume}")
    print("目标单元平均密度:", [str(v) for v in result.target_average_density])
    print("未覆盖体积:", [str(v) for v in result.target_uncovered_volume])
    print(f"目标总质量 {result.total_target_mass}  "
          f"覆盖部分源质量 {result.covered_source_mass}  "
          f"总量误差 {result.mass_error}")
    assert all(v == Fraction(7, 2) for v in result.target_average_density)
    assert result.mass_error == 0

    print()
    print("=" * 68)
    print("B. 跨界分段常量：左半密度 1、右半密度 3，目标为中间条")
    print("=" * 68)
    left = TetMesh.build(box_vertices(0, Fraction(1, 2)), CUBE5_TETS)
    right = TetMesh.build(box_vertices(Fraction(1, 2), 1), CUBE5_TETS)
    src2 = TetMesh.build(
        list(left.vertices) + list(right.vertices),
        list(left.tets) + [tuple(k + 8 for k in t) for t in right.tets],
    )
    tgt2 = TetMesh.build(
        box_vertices(Fraction(1, 4), Fraction(3, 4)), CUBE5_TETS
    )
    result2 = remap(src2, tgt2, [1] * 5 + [3] * 5)
    print("目标单元平均密度:", [str(v) for v in result2.target_average_density])
    print(f"目标总质量 {result2.total_target_mass}（手算期望 1 = 1/4*1 + 1/4*3）")
    print(f"覆盖部分源质量 {result2.covered_source_mass}  总量误差 {result2.mass_error}")
    assert result2.total_target_mass == 1
    assert result2.mass_error == 0

    print()
    print("=" * 68)
    print("C. 目标超出源区域（z 到 2）：缺口保留，不归一")
    print("=" * 68)
    tgt3 = TetMesh.build(box_vertices(0, 1, 0, 1, 0, 2), CUBE5_TETS)
    result3 = remap(source, tgt3, [2] * 5)
    print("目标单元覆盖体积:  ", [str(v) for v in result3.target_covered_volume])
    print("目标单元未覆盖体积:", [str(v) for v in result3.target_uncovered_volume])
    print("目标单元平均密度:  ", [str(v) for v in result3.target_average_density])
    print(f"目标总质量 {result3.total_target_mass}（源总质量 2，未放大到满覆盖）")
    print(f"总缺口 {sum(result3.target_uncovered_volume, Fraction(0))}（上半单位立方体）")
    assert result3.total_target_mass == 2
    assert sum(result3.target_uncovered_volume, Fraction(0)) == 1

    elapsed = time.perf_counter() - start
    print()
    print(f"全部核验通过，真实计算耗时 {elapsed:.3f} 秒（预算 90 秒）。")


if __name__ == "__main__":
    main()
