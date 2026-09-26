# orient_audit — 四面体网格的定向与边界审计

离线、纯标准库的可复用题库：输入有理坐标顶点与最多 500 个四面体，
精确判定并修复单元定向，建立面关联，报告边界面、非流形、连通分量、
总有向体积与缺陷位置。

增量能力（守恒重映射）：输入两张覆盖（部分）同一区域的四面体网格与
源单元常量密度，用半空间交精确计算每对源/目标四面体的真实交叠体积，
输出稀疏传递权重、每个目标单元的平均密度、未覆盖体积（缺口保留，
不归一）与总量守恒误差。

- **Python 版本**：CPython 3.14.7（仅用到 3.10+ 即有的语法特性）。
- **依赖**：无第三方依赖，只用标准库（`fractions`、`dataclasses`、
  `unittest` 等）。装好 Python 后无需安装任何其他东西。
- **离线**：核心、测试、演示全部离线运行，不访问网络、文件、数据库、
  云 API 或任何中间件；示例数据由仓库内确定性代码直接给出。

## 环境准备（隔离 venv）

```bash
python3.14 -m venv .venv          # 隔离环境；仅用标准库，无需 pip install
source .venv/bin/activate
python --version                   # 应为 3.14.7
```

## 准确命令

```bash
# 必要测试（在仓库根目录执行；基线回归 + 重映射新测试）
python -m unittest discover -s tests -v

# 基线演示：定向审计与修复（约 0.05 秒）
python demo.py

# 增量演示：两张网格间的守恒重映射（约 0.3 秒，远小于 90 秒预算；
# 真实离线计算，无账号/密钥/sleep/预录结果）
python demo_iteration.py
```

## 公开接口与数据结构

```python
from orient_audit import (
    TetMesh, audit, repair_orientation,      # 核心入口
    AuditReport, BoundaryFace, InteriorFace, NonManifoldFace,  # 结果类型
    MeshInputError, DuplicateTetrahedraError, ZeroVolumeError, # 异常
    MAX_TETS,                                # 500
    # 守恒重映射（增量）
    remap, RemapResult, TransferWeight,
    intersection_volume, tet_volume,
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

## 守恒重映射：`remap(source, target, source_density) -> RemapResult`

- **输入**：两个已校验的 :class:`TetMesh`（只读，不被修改；仅凸四面体
  单元，任一网格含零体积单元会带定位报错）与逐源单元常量密度
  （有理数序列，拒绝 float/bool，长度必须等于源单元数）。
- **交叠体积**：对每对源/目标四面体取共 8 个外法向半空间，枚举平面
  三元组求交点、保留落在全部半空间内的点作为凸交体顶点，再对每张
  平面上的顶点按外法向排序、扇形三角化，用散度定理公式
  ``V = (1/6)·Σ det(p,q,r)`` 精确求和。同平面异侧的半空间不视为重复；
  退化交集自然得到 0。全部为 Fraction 精确运算，无容差。
- **输出**（`RemapResult`，全部精确 Fraction）：
  - `weights`：稀疏传递权重 `TransferWeight(source_index, target_index,
    volume)`，只存交叠体积 > 0 的项；
  - `target_mass` / `target_average_density`：每个目标单元的质量与
    平均密度，平均按**整个**目标单元体积计算（未覆盖部分按 0 贡献）；
  - `target_covered_volume` / `target_uncovered_volume`：覆盖与缺口
    体积，目标未被源完全覆盖时缺口原样保留，**不归一到满覆盖**；
  - `total_target_mass`、`covered_source_mass`（按源侧独立汇总）、
    `mass_error`（两者之差，精确算术中恒为 0，作守恒核验）。
- **不是**形心采样：权重来自真实交叠体积，目标中心落在哪个源单元
  与结果无关。

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
- 重映射核验：同一立方体 5 剖分 -> 6 剖分的常量场精确保持；
  跨界分段常量（左半密度 1、右半密度 3）总质量手算 = 1 精确吻合，
  单个跨界四面体的左右体积 7/96 与 1/96 为手算值；目标超出源区域时
  缺口精确等于上半单位立方体体积 1，且总质量不被放大。
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
- 重映射只支持凸四面体单元的真实交叠；不实现一般非凸单元、
  高阶场或网格随时间运动。两网格不必覆盖同一区域，未覆盖部分
  以缺口形式报告。

## 目录

```
orient_audit/
  __init__.py   # 公开接口
  errors.py     # MeshInputError / DuplicateTetrahedraError / ZeroVolumeError
  exact.py      # 有理坐标解析、精确 3x3 行列式、输入上限
  mesh.py       # TetMesh、无序三点面关联
  audit.py      # audit / repair_orientation / AuditReport
  remap.py      # 守恒重映射：半空间交叠体积、稀疏权重、缺口与守恒误差
tests/
  fixtures.py   # 确定性立方体剖分数据
  reference.py  # 独立参考（4x4 Bareiss 行列式、面通量、外法向检查）
  test_audit.py # 基线回归测试
  test_remap.py # 重映射测试（手算核对、守恒、缺口、非法输入）
demo.py         # 基线固定输入演示
demo_iteration.py  # 重映射固定输入演示
```
