#!/usr/bin/env python3
"""线段穿越采样与积分演示（本轮新增入口）。

真实计算，无 sleep、无预录结果、无网络。两部分：

A. 正常结果：单位立方体 5 四面体剖分 + 仿射场 f = 2x - y + 3z + 1/2。
   打印穿越序列、各段参数区间与插值端点值；总积分用两条独立路径核对：
   (1) 仿射场的解析原函数积分；(2) 逐四面体独立裁线（叉积外法向
   Cyrus–Beck，不调用被测核心）核对区间并集。

B. 边界结果：对角线恰过两个共享顶点，网外间隙单独列出、过点不计长；
   负定向单元网格被拒绝（InvalidTraversalMeshError），按原审计规则
   repair_orientation 修复后再穿越成功。

运行：python demo_iteration.py
"""

from __future__ import annotations

import time
from fractions import Fraction

from orient_audit import (
    InvalidTraversalMeshError,
    TetMesh,
    repair_orientation,
    traverse_segment,
)

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

# 仿射场 f(x,y,z) = 2x - y + 3z + 1/2。
AFFINE = (Fraction(2), Fraction(-1), Fraction(3), Fraction(1, 2))


def affine_at(p):
    return AFFINE[0] * p[0] + AFFINE[1] * p[1] + AFFINE[2] * p[2] + AFFINE[3]


def affine_integral(p0, p1, ta, tb):
    """解析参考：∫ f(p0 + t(p1-p0)) dt（一次多项式原函数）。"""

    d = tuple(p1[k] - p0[k] for k in range(3))
    slope = sum(AFFINE[k] * d[k] for k in range(3))
    intercept = affine_at(p0)
    return slope * (tb * tb - ta * ta) / 2 + intercept * (tb - ta)


def independent_clip(vertices, tet, p0, p1):
    """独立裁线参考：叉积外法向 + 点积半空间裁剪（不调用被测核心）。"""

    def sub(a, b):
        return tuple(a[k] - b[k] for k in range(3))

    def cross(u, v):
        return (u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0])

    def dot(u, v):
        return sum(u[k] * v[k] for k in range(3))

    pts = [tuple(Fraction(c) for c in vertices[k]) for k in tet]
    a = tuple(Fraction(c) for c in p0)
    d = sub(tuple(Fraction(c) for c in p1), a)
    lo, hi = Fraction(0), Fraction(1)
    for face in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (2, 1, 0)):
        p, q, r = (pts[face[0]], pts[face[1]], pts[face[2]])
        n = cross(sub(q, p), sub(r, p))
        opp = pts[({0, 1, 2, 3} - set(face)).pop()]
        if dot(n, sub(opp, p)) > 0:
            n = tuple(-x for x in n)
        c0, c1 = dot(n, sub(a, p)), dot(n, d)
        if c1 > 0:
            hi = min(hi, -c0 / c1)
        elif c1 < 0:
            lo = max(lo, -c0 / c1)
        elif c0 > 0:
            return None
        if lo > hi:
            return None
    return (lo, hi)


def independent_union(mesh, p0, p1):
    ivs = []
    for tet in mesh.tets:
        r = independent_clip(mesh.vertices, tet, p0, p1)
        if r is not None and r[0] < r[1]:
            ivs.append(r)
    ivs.sort()
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [tuple(x) for x in merged]


def show(result):
    for p in result.pieces:
        print(f"  单元 {p.tet_index}: t ∈ [{p.t0}, {p.t1}]"
              f"  端点值 ({p.value0}, {p.value1})  段积分 {p.integral}")
    for g in result.gaps:
        print(f"  网外间隙: t ∈ [{g.t0}, {g.t1}]（不计入积分，不补 0）")
    print(f"  穿越序列: {list(result.tet_sequence)}")
    print(f"  总积分（精确）: {result.total_integral}"
          f"  ≈ {result.total_integral_float():.6f}")


def main() -> None:
    t0 = time.perf_counter()
    mesh = TetMesh.build(VERTICES, TETS)
    values = [affine_at(v) for v in VERTICES]

    print("=" * 68)
    print("A. 立方体内线段 + 仿射场：穿越采样与积分")
    print("=" * 68)
    p0 = (Fraction(1, 8), Fraction(1, 8), Fraction(1, 8))
    p1 = (Fraction(7, 8), Fraction(7, 8), Fraction(7, 8))
    result = traverse_segment(mesh, values, p0, p1)
    show(result)
    analytic = affine_integral(p0, p1, Fraction(0), Fraction(1))
    union = independent_union(mesh, p0, p1)
    print(f"  解析积分（独立参考）: {analytic}")
    print(f"  独立逐四面体裁线并集: {union}")
    assert result.total_integral == analytic
    assert [p.t0 for p in result.pieces][0] == 0
    assert [p.t1 for p in result.pieces][-1] == 1
    assert all(a.t1 == b.t0 for a, b in zip(result.pieces, result.pieces[1:]))
    assert union == [(Fraction(0), Fraction(1))]
    print("  总积分与解析参考一致；区间并集与独立裁线一致。")

    print()
    print("=" * 68)
    print("B. 边界：恰过共享顶点的对角线 / 无效网格拒绝")
    print("=" * 68)
    q0, q1 = (-1, -1, -1), (2, 2, 2)
    diag = traverse_segment(mesh, values, q0, q1)
    print("[B1] 对角线 (-1,-1,-1)->(2,2,2)，恰过顶点 (0,0,0) 与 (1,1,1)：")
    show(diag)
    assert diag.total_integral == affine_integral(
        q0, q1, Fraction(1, 3), Fraction(2, 3)
    )
    assert [(g.t0, g.t1) for g in diag.gaps] == [
        (Fraction(0), Fraction(1, 3)), (Fraction(2, 3), Fraction(1))
    ]
    print("  顶点处只计一次长度；间隙单独列出；积分与解析参考一致。")

    print("[B2] 负定向单元网格：")
    broken = [list(t) for t in TETS]
    broken[2] = [TETS[2][0], TETS[2][1], TETS[2][3], TETS[2][2]]
    bad_mesh = TetMesh.build(VERTICES, broken)
    try:
        traverse_segment(bad_mesh, values, p0, p1)
    except InvalidTraversalMeshError as exc:
        print(f"  traverse_segment 拒绝 -> {type(exc).__name__}:")
        print(f"    {exc}")
    fixed = repair_orientation(bad_mesh)
    ok = traverse_segment(fixed, values, p0, p1)
    print(f"  repair_orientation 修复后穿越成功: 序列 {list(ok.tet_sequence)}，"
          f"总积分 {ok.total_integral}")
    assert ok.total_integral == analytic

    elapsed = time.perf_counter() - t0
    print()
    print(f"演示完成，实际耗时 {elapsed:.3f} 秒（预算 90 秒）。")


if __name__ == "__main__":
    main()
