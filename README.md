# orient_audit — 四面体网格的定向与边界审计 + 有理平面截面

离线、纯标准库的可复用题库：输入有理坐标顶点与最多 500 个四面体，
精确判定并修复单元定向，建立面关联，报告边界面、非流形、连通分量、
总有向体积与缺陷位置；并可对网格做**有理平面截面**，输出精确有理
多边形片（三角网与平面的交）。

- **Python 版本**：CPython 3.14.7（仅用到 3.10+ 即有的语法特性）。
- **依赖**：无第三方依赖，只用标准库（`fractions`、`dataclasses`、
  `unittest` 等）。装好 Python 后无需安装任何其他东西。
- **离线**：核心、测试、演示全部离线运行，不访问网络、文件、数据库、
  云 API 或任何中间件；示例数据由仓库内确定性代码直接给出。

## 准确命令

```bash
# 必要测试（在仓库根目录执行；回归入口）
./.venv/bin/python -m unittest discover -s tests -v

# 固定输入演示（约 0.1 秒，远小于 90 秒预算；真实计算，无 sleep/预录）
./.venv/bin/python demo.py
```

`.venv` 是仓库旁的 Python 3.14.7 隔离环境；用其它解释器时把
`./.venv/bin/python` 换成对应 `python` 即可，无额外安装步骤。

## 公开接口与数据结构

```python
from orient_audit import (
    TetMesh, audit, repair_orientation,      # 核心入口
    section,                                 # 平面截面入口
    AuditReport, BoundaryFace, InteriorFace, NonManifoldFace,  # 审计结果类型
    SectionReport, SectionPatch, SectionPoint,                 # 截面结果类型
    PointContact, EdgeContact,
    MeshInputError, DuplicateTetrahedraError, ZeroVolumeError, # 异常
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

## 平面截面 `section(mesh, normal, offset) -> SectionReport`（只读，不修改输入）

计算四面体网格与有理平面的交，输出精确有理多边形片。

### 参数与错误状态

- `mesh`：必须是 `TetMesh`（否则 `TypeError`）。
- `normal`：平面法向，恰好 3 个有理分量（`int` / `Fraction` / 有理字符串，
  **拒绝 `float`**）。零向量抛 `MeshInputError`（`location="plane.normal"`）；
  非法分量可定位（如 `plane.normal[2]`）。
- `offset`：平面偏置（同上有理类型），平面方程为 **`normal·x = offset`**；
  非法值定位 `plane.offset`。

### 结果语义

- **统一分类**：每个顶点按精确符号 `sign(normal·v − offset)` 分为
  正 / 负 / 恰好在平面上，全程 `Fraction`，无容差。
- **点身份不靠浮点焊接**：截面点 `SectionPoint(pid, coords, origin)`
  的身份是 `("vertex", i)`（就是原顶点 `i`，坐标对象恒等保留）或
  `("edge", i, j)`（原网格棱 `(i,j)`，`i<j` 的交点）。两个单元共享
  一条棱时，交点天然是同一个全局 ID。
- **面片 `SectionPatch`**：`kind` 为 `"triangle"` / `"quad"` /
  `"coplanar_face"`；`polygon` 为有序边界点全局 ID；`triangles` 为
  三角剖分；`source_tets` 为来源单元；`area2` 为**精确面积平方**
  （`Fraction`），`area()` 是浮点便利出口。
  - 每个非退化四面体贡献一个三角形（1-3 分布）或四边形（2-2 分布）；
  - **四边形按全局交点 ID 选确定对角线**：含四个全局 ID 中最小者的
    那条对角线，规则与遍历顺序无关；
  - **完全共面的原面只保留一份**（`kind="coplanar_face"`），
    `source_tets` 记录所有来源单元。
- **零面积接触单独列出**：只在原顶点处碰平面的单元进
  `point_contacts`，只沿原棱贴在平面上的进 `edge_contacts`，
  均不混进 `patches`。
- **朝向统一**：所有面片 / 三角形的法向与给定 `normal` 点积为正。
- **负向单元**：先按原修正规则（只交换最后两个顶点，与
  `repair_orientation` 相同）在内部修正再截面，修正的单元列入
  `repaired_tets`；输入网格本身不被修改。
- **退化单元**：零体积四面体不参与截面，列入 `degenerate_tets`。
- **不改网格、不分裂边**：截面点是独立新建的有理点，原网格的顶点
  与单元对象恒等保留。

### 截面面积的误差口径

精确出口是每片的 `area2`（面积平方，`Fraction`，精确）与
`SectionReport.plane_residuals()`（应全为精确 0）。浮点便利出口
`SectionPatch.area()` / `SectionReport.total_area()` 是对精确值做
IEEE-754 binary64 开方 / 求和，单次舍入误差 ≤ 0.5 ULP 量级，
不参与任何几何判定。

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
- 已覆盖：翻转单元修复、重复体（含排列重复）、共面退化拒绝、
  三单元共面的非流形、共享面同向冲突、多个连通分量、500 上限边界、
  bool/float/NaN/Infinity/越界等非法输入定位。
- **截面独立核对**：测试用 Cramer 法则重心坐标（`tests/reference.py`
  的 `ref_barycentric`，不调用被测核心）验证每片面片质心落在全部
  来源单元内；逐点验证精确平面方程残差为 0；并用**解析棱锥**
  （方底棱锥在 z=1/2 处截面是边长 1 的正方形）核对截面总面积；
  另验证共享棱交点全局 ID 一致、共面面去重、四边形对角线规则、
  负向单元修正前后结果一致、输入网格对象恒等不变。

## 已知限制（本题有限协议，不宣称行业标准兼容）

- 坐标必须有理；不支持浮点网格、无理数或度量单位换算。
- 只处理四面体（4 顶点）单元；不支持混合单元、棱柱或多面体。
- 非流形的判定只按“同一三点面被 >2 个单元引用”；不检测悬点、
  自相交、孔洞规范、网格标签/材料属性等更广义的网格质量问题。
- 修复只做逐单元体积极性，不做全局一致性传播（例如彼此孤立、无法靠
  几何符号定向的零体积簇保持可报告但不可修复——这正是契约要求）。
- 截面只输出每个单元与平面的交（三角形 / 四边形 / 共面面 / 零面积
  接触），不做跨单元的多边形合并、洞检测或拓扑缝合；共享棱交点靠
  身份键天然对齐，已足够消费方自行拼接。
- 截面面积的精确出口是面积平方（有理点的多边形面积一般不是有理数）；
  浮点 `area()` 仅为便利。
- 规模上限 500 个四面体（接口边界，非性能极限）。

## 目录

```
orient_audit/
  __init__.py   # 公开接口
  errors.py     # MeshInputError / DuplicateTetrahedraError / ZeroVolumeError
  exact.py      # 有理坐标解析、精确 3x3 行列式、输入上限
  mesh.py       # TetMesh、无序三点面关联
  audit.py      # audit / repair_orientation / AuditReport
  section.py    # section / SectionReport（有理平面截面）
tests/
  fixtures.py   # 确定性立方体剖分数据
  reference.py  # 独立参考（4x4 Bareiss 行列式、面通量、外法向、重心坐标）
  test_audit.py # 审计必要测试
  test_section.py # 截面必要测试
demo.py         # 固定输入演示（审计 + 截面）
```
