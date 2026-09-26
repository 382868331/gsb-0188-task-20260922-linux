# orient_audit — 四面体网格的定向与边界审计

离线、纯标准库的可复用题库：输入有理坐标顶点与最多 500 个四面体，
精确判定并修复单元定向，建立面关联，报告边界面、非流形、连通分量、
总有向体积与缺陷位置。第二轮新增：给定每顶点标量场，沿线段精确
采样与积分——返回线段依次落入的四面体、各段参数区间、线性插值
端点值与总积分，网格外间隙单独列出。

- **Python 版本**：CPython 3.14.7（仅用到 3.10+ 即有的语法特性）。
- **依赖**：无第三方依赖，只用标准库（`fractions`、`dataclasses`、
  `unittest` 等）。装好 Python 后无需安装任何其他东西。
- **离线**：核心、测试、演示全部离线运行，不访问网络、文件、数据库、
  云 API 或任何中间件；示例数据由仓库内确定性代码直接给出。

## 准确命令

# venv 准备（可选；本库无第三方依赖，venv 仅为隔离解释器）
python3.14 -m venv .venv && . .venv/bin/activate

```bash
# 必要测试（在仓库根目录执行）
python -m unittest discover -s tests -v

# 固定输入演示（约 0.05 秒，远小于 8 秒预算；真实计算，无 sleep/预录）
python demo.py

# 第二轮演示：线段穿越采样与积分（真实计算，远小于 90 秒预算）
python demo_iteration.py
```

## 公开接口与数据结构

```python
from orient_audit import (
    TetMesh, audit, repair_orientation,      # 核心入口
    sample_segment,                          # 线段穿越采样与积分
    AuditReport, BoundaryFace, InteriorFace, NonManifoldFace,  # 结果类型
    SegmentTraversal, SegmentHit, SegmentGap,                  # 穿越结果类型
    MeshInputError, DuplicateTetrahedraError, ZeroVolumeError, # 异常
    TraversalError,                          # 网格不可行走
    MAX_TETS,                                # 500
)
```

### 构造

`TetMesh.build(vertices, tets) -> TetMesh`（不可变对象；构造失败不产生
部分结果）：

- `vertices`：三维顶点序列，每个顶点是 3 个**有理**坐标。接受 `int`、
  `fractions.Fraction` 与有理字符串（如 `"1/3"`）。**拒绝 `float`**，
  从源头消除浮点误差；因此 NaN/Infinity 在 float 入口即被拒，字符串
  `"nan"`/`"inf"`/`"1/0"` 同样拒绝，错误带定位（如 `verts[2][1]`）。
- `tets`：单元序列，每个单元是恰好 4 个**互不相同**的顶点下标（`int`，
  **拒绝 `bool`**，越界可定位，如 `tets[3][0]`）。
- **重复单元单独报错**：四个顶点集合相同即视为重复（与排列顺序无关），
  抛 `DuplicateTetrahedraError`，其 `duplicates` 列出全部重复组
  `(排序后的顶点元组, (单元下标...))`。

### `audit(mesh) -> AuditReport`（只读，不修改输入）

字段：

| 字段 | 含义 |
|---|---|
| `n_tets` | 单元数 |
| `negative_tets` | 精确行列式 < 0 的单元下标（升序） |
| `zero_tets` | 精确体积 = 0（共面退化）的单元下标（升序） |
| `total_signed_volume6` | 全部单元 6 倍有向体积的精确和（`Fraction`） |
| `total_signed_volume` | 精确有向体积和（属性，`Fraction`） |
| `boundary_faces` | 边界面 `BoundaryFace(ordered, tet_index)`，`ordered` 为单元诱导朝向 |
| `interior_faces` | 恰被两单元引用的面 `InteriorFace(key, tet_indices, opposite)` |
| `non_manifold_faces` | 被 >2 单元引用的面 `NonManifoldFace(key, tet_indices)` |
| `components` | 共享任一三点面的单元连通分量（并查集），按下标排序 |

派生判定属性：`consistently_oriented`（无负向单元且内面全部反向）、
`is_manifold`、`has_zero_volume`、`has_non_manifold`、
`inconsistent_faces`（两单元共享但诱导同向的内面）。

### `repair_orientation(mesh) -> TetMesh`

- 负体积单元**只交换最后两个顶点**（局部位 2 与 3，奇置换使行列式变号）；
  正体积单元（含偶置换排列）保持原样；
- **零体积单元拒绝修复**：存在任一零体积单元时抛 `ZeroVolumeError`
  （`tet_indices` 列出全部），不返回部分修改，原网格保持不变；
- **修复不移动顶点**：返回新网格，顶点坐标对象恒等保留（有断言守护）；
- 非流形/定向冲突属于拓扑状态，由 `audit` 报告，不阻止逐单元定向修复。

## 面关联语义（无序三点键）

每个面取排序后的三个顶点下标作为键。一个四面体的四个诱导面（外法向
朝向）固定为：

```
(1,2,3), (0,3,2), (0,1,3), (2,1,0)
```

- 键出现 **1 次**：边界面；
- 出现 **2 次且两诱导朝向相反**（排列奇偶性不同）：合法内面；
- 出现 **2 次但同向**：定向冲突内面（`opposite=False`）；
- 出现 **>2 次**：非流形面。

## 浮点结果误差口径

## 线段穿越采样与积分（sample_segment）

```python
tr = sample_segment(mesh, values, start, end)
```

- `values`：每顶点一个标量（int / Fraction / 有理字符串；拒绝 float/bool，
  长度必须等于顶点数，错误带定位如 `values[2]`）。
- `start` / `end`：线段端点，各 3 个有理坐标；参数 `t=0` 在 start、
  `t=1` 在 end。两端点重合抛 `MeshInputError`。
- 返回 `SegmentTraversal`：
  - `hits`：`SegmentHit(tet_index, t0, t1, value0, value1)` 按 t 升序、
    互不重叠；`value0/value1` 为分段线性场在段端点的精确插值，
    `hit.integral` 为该段积分（线性场梯形公式即精确值）；
  - `gaps`：`SegmentGap(t0, t1)` 网格外区间，单独列出，**不自动补 0
    混进结果**；`hits` 与 `gaps` 互补地覆盖 `[0, 1]`；
  - `total_integral`：各命中段积分之和（`Fraction` 精确值）；
    `total_integral_float()` 为唯一浮点出口（一次正确舍入，≤0.5 ULP）。

### 行走策略与退化规则

- **跨面行走不重扫全网格**：从当前单元经面邻接表跨到相邻单元；全网格
  扫描只发生在初始定位与离开网格后的下一次进入定位（且用的是构造期
  一次性预计算的区间表）。
- **线恰过共享棱/顶点**：同一时刻多个面同时到期时，候选下一单元取所有
  到期面的邻居中"区间覆盖 t 之后"者，按**单元下标最小**唯一选定；
  邻居表无合法候选（线穿顶点进入仅共点的单元）时回退到预计算区间表，
  同样取下标最小者。长度为零的接触（擦过顶点/棱）不产生命中段，
  同一长度不会计两次。
- **沿边界/共享面行走的区间归属唯一**：每个参数 t 至多属于一个命中段，
  由"下标最小合法单元"规则确定归属。
- **无效网格拒绝**：存在零体积单元、负定向单元（需先
  `repair_orientation`）或非流形面时抛 `TraversalError`
  （`problems` 列出全部可定位问题），不输出看起来完整的穿越序列。
  负定向单元仍按原审计规则处理：`audit` 报告、`repair_orientation`
  修复后再穿越。

### 独立核对

- **仿射场解析积分**：测试与 demo 用 ∫(α·X+β)dt 的解析公式
  （不经过任何单元/重心计算）核对 `total_integral`。
- **独立逐四面体裁线**：`tests/reference.py` 的 `ref_clip_intervals`
  对每个四面体独立解半空间不等式（叉积法向，不用面邻接行走、不 import
  本库），核对命中区间的并集精确一致、间隙互补覆盖 `[0,1]`、
  每个命中段确属其声明的单元。


所有几何/拓扑判定（体积正负、零体积、面朝向、非流形）全部用
`fractions.Fraction` 精确完成，不存在容差，也**不存在搜索预算耗尽**
这类不确定结果——本算法是确定性的一次遍历构造。

唯一的浮点出口是 `AuditReport.total_signed_volume_float()`，它对精确
`Fraction` 做一次 IEEE-754 binary64 正确舍入：误差 ≤ 0.5 ULP；
对量级正常的非零值相对误差 ≤ 2⁻⁵³（约 1.11×10⁻¹⁶）。该舍入误差不
参与、也不可能掩盖任何边界判定。默认建议直接使用精确的
`total_signed_volume`（`Fraction`）。

## 设计取舍

- **拒绝 float 输入而非内部加容差**：题面要求精确行列式判体积、零体积
  必须可证，浮点配阈值会使“恰好共面”依赖阈值。有理输入 + Fraction
  让正负与为零成为精确事实。
- **修复策略固定为“交换最后两顶点”**：不做全局定向传播、不改变顶点
  坐标。逐单元独立翻转可在输入存在非流形/互不接触分量时依然定义良好；
  内面是否因此一致由 `audit` 如实报告。
- **重复单元在构造期拒绝**：重复体会让面的关联计数虚增（例如把合法内面
  伪装成非流形），因此必须作为独立错误先于几何处理报出，而不是混进审计。
- **`audit` 宽容、`repair` 严格**：审计如实报告零体积单元；修复对零体积
  直接拒绝（交换顶点不可能让零行列式非零）。
- 连通分量按“共享完整三点面”的并查集计算；仅共享棱或点不算相连。

## 验收要点与独立核对

- 单位立方体标准 5 四面体剖分（中心正四面体 1/3 + 四角各 1/6）：
  精确总体积 = 1，12 个边界面，4 个内面全部诱导反向，单一连通分量。
- **面通量独立核对**：测试与 demo 另用散度定理公式
  `V = (1/6)·Σ det(p,q,r)` 对边界面直接求和，该实现不调用被测核心、
  不经过单元体积，独立得到 1。
- **体积独立核对**：测试参考实现用 4×4 齐次坐标行列式 + Bareiss 无分数
  消元（`tests/reference.py`，不 import 本库），与核心的 3×3 余子式
  展开是两条独立算法路径；并独立用叉积法向·(面心−体形心) 检查边界面朝外。
- 已覆盖：翻转单元修复、重复体（含排列重复）、共面退化拒绝、
  三单元共面的非流形、共享面同向冲突、多个连通分量、500 上限边界、
  bool/float/NaN/Infinity/越界等非法输入定位。

## 已知限制（本题有限协议，不宣称行业标准兼容）

- 坐标必须有理；不支持浮点网格、无理数或度量单位换算。
- 只处理四面体（4 顶点）单元；不支持混合单元、棱柱或多面体。
- 非流形的判定只按“同一三点面被 >2 个单元引用”；不检测悬点、
  自相交、孔洞规范、网格标签/材料属性等更广义的网格质量问题。
- 修复只做逐单元体积极性，不做全局一致性传播（例如彼此孤立、无法靠
  几何符号定向的零体积簇保持可报告但不可修复——这正是契约要求）。
- 规模上限 500 个四面体（接口边界，非性能极限）。
- 线段穿越只支持四面体网格上的分段线性标量场；曲线、高阶场、
  非流形/退化网格不在契约内（后者直接拒绝而非近似穿越）。

## 目录

```
orient_audit/
  __init__.py   # 公开接口
  errors.py     # MeshInputError / DuplicateTetrahedraError / ZeroVolumeError
  exact.py      # 有理坐标解析、精确 3x3 行列式、输入上限
  mesh.py       # TetMesh、无序三点面关联
  audit.py      # audit / repair_orientation / AuditReport
  traverse.py   # sample_segment：线段穿越采样与积分
tests/
  fixtures.py   # 确定性立方体剖分数据
  reference.py  # 独立参考（4x4 Bareiss 行列式、面通量、外法向检查、
                #              逐四面体裁线、仿射场解析积分）
  test_audit.py # 必要测试
  test_traverse.py # 线段穿越必要测试
demo.py         # 固定输入演示
demo_iteration.py # 第二轮演示：线段穿越采样与积分
```
