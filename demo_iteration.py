#!/usr/bin/env python3
"""固定输入演示：线段穿越四面体网格的采样与积分（第二轮）。

真实计算，无 sleep、无预录结果、无网络。三部分：

A. 正常结果：单位立方体 5 四面体剖分，每顶点标量取仿射场
   f(x,y,z) = 2x + 3y - z + 5 的精确值，一条斜线段从网格外穿入再穿出。
   打印依次落入的单元、各段参数区间、线性插值端点值、网格外间隙与
   总积分；总积分另用*独立*的仿射场解析积分公式核对（不经过任何
   单元/重心计算）。

B. 边界结果：线段恰过共享顶点 0（进入点同时属于 4 个单元，按确定
   规则唯一归属、长度不计两次），以及一段完全在网格外（整段为间隙，
   不补 0）。

C. 实际触发的拒绝：网格含负定向单元时 sample_segment 抛
   TraversalError，不输出看起来完整的穿越序列；按原审计规则
   repair_orientation 修复后同一调用成功。

运行：python demo_iteration.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import (
    TetMesh,
    TraversalError,
    repair_orientation,
    sample_segment,
)

# ---------------------------------------------------------------- 固定数据
VERTICES = (
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
)
TETS = (
    (0, 1, 2, 5),
    (0, 2, 3, 7),
    (0, 4, 5, 7),
    (2, 5, 6, 7),
    (0, 2, 7, 5),
)
ALPHA = (Fraction(2), Fraction(3), Fraction(-1))
BETA = Fraction(5)
VALUES = tuple(
    ALPHA[0] * x + ALPHA[1] * y + ALPHA[2] * z + BETA for x, y, z in VERTICES
)


def affine_integral_reference(start, end, intervals):
    """独立解析参考：∫ (α·X + β) dt，X(t) = a + t·d（不调用本库）。"""

    a = tuple(Fraction(c) for c in start)
    d = tuple(Fraction(e) - Fraction(s) for s, e in zip(start, end, strict=True))
    base = sum(al * ai for al, ai in zip(ALPHA, a, strict=True)) + BETA
    slope = sum(al * di for al, di in zip(ALPHA, d, strict=True))
    return sum(
        base * (hi - lo) + slope * (hi * hi - lo * lo) / 2 for lo, hi in intervals
    )


def show_traversal(title, mesh, start, end):
    print(f"--- {title} ---")
    print(f"线段: {start} -> {end}")
    tr = sample_segment(mesh, VALUES, start, end)
    for h in tr.hits:
        print(
            f"  命中 tets[{h.tet_index}]  t ∈ [{h.t0}, {h.t1}]"
            f"  端点值 ({h.value0}, {h.value1})  段积分 {h.integral}"
        )
    for g in tr.gaps:
        print(f"  网格外间隙  t ∈ [{g.t0}, {g.t1}]（不补 0，不计入积分）")
    print(f"  总积分（精确）: {tr.total_integral}")
    print(f"  总积分（float 视图）: {tr.total_integral_float()!r}")
    return tr


def main() -> None:
    t_start = time.perf_counter()
    mesh = TetMesh.build(VERTICES, TETS)

    # A. 正常穿越：网格外 -> 中心四面体 -> 角四面体 -> 网格外。
    start, end = (Fraction(-1),) * 3, (Fraction(2),) * 3
    tr = show_traversal("A. 斜线段穿越立方体（仿射场）", mesh, start, end)
    expect = affine_integral_reference(
        start, end, [(h.t0, h.t1) for h in tr.hits]
    )
    assert tr.total_integral == expect, "总积分与独立解析参考不一致"
    print(f"  独立解析积分核对: {expect}  ✓ 一致")

    # B1. 边界：线段恰过共享顶点 0（4 个单元共点，确定规则唯一归属）。
    show_traversal(
        "B1. 线段恰过共享顶点 0",
        mesh,
        (Fraction(-1), Fraction(0), Fraction(0)),
        (Fraction(1), Fraction(1), Fraction(1)),
    )

    # B2. 边界：整段在网格外，只有间隙。
    show_traversal(
        "B2. 整段在网格外",
        mesh,
        (Fraction(2), Fraction(2), Fraction(2)),
        (Fraction(3), Fraction(3), Fraction(3)),
    )

    # C. 拒绝路径：负定向单元 -> TraversalError；修复后成功。
    print("--- C. 负定向网格拒绝穿越 ---")
    bad_tets = TETS[:4] + ((0, 2, 5, 7),)  # 中心单元交换最后两顶点
    bad_mesh = TetMesh.build(VERTICES, bad_tets)
    try:
        sample_segment(bad_mesh, VALUES, (Fraction(0),) * 3, (Fraction(1),) * 3)
    except TraversalError as exc:
        print(f"  已拒绝（符合契约）: {exc}")
    repaired = repair_orientation(bad_mesh)
    tr = sample_segment(repaired, VALUES, (Fraction(0),) * 3, (Fraction(1),) * 3)
    print(f"  repair_orientation 后同一调用成功：{len(tr.hits)} 段命中，"
          f"总积分 {tr.total_integral}")

    elapsed = time.perf_counter() - t_start
    print(f"全部演示真实计算完成，用时 {elapsed:.3f} 秒（预算 90 秒）")


if __name__ == "__main__":
    main()

