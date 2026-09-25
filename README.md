# orient_audit — 四面体网格的定向审计与有理平面剖切

离线、纯标准库的可复用题库：输入有理坐标顶点与最多 500 个四面体，
精确判定并修复单元定向，建立面关联，报告边界面、非流形、连通分量、
总有向体积与缺陷位置；并以有理平面 ``n·x = d`` 剖切网格，输出三角网
与平面相交的**有理**多边形片（三角/四边形）、共面原面与零面积接触。

- **Python 版本**：CPython 3.14.7（仅用到 3.10+ 即有的语法特性）。
- **依赖**：无第三方依赖，只用标准库（`fractions`、`dataclasses`、
  `unittest` 等）。装好 Python 后无需安装任何其他东西。
- **离线**：核心、测试、演示全部离线运行，不访问网络、文件、数据库、
  云 API 或任何中间件；示例数据由仓库内确定性代码直接给出。

## 准确命令

```bash
# 隔离环境（仓库已附带指向内置 Python 3.14.7 的 .venv 符号链接）
./.venv/bin/python --version        # Python 3.14.7

# 必要测试（在仓库根目录执行；保留的回归入口）
./.venv/bin/python -m unittest discover -s tests -v

# 固定输入演示（约 0.01 秒，远小于 90 秒预算；真实计算，无 sleep/预录）
./.venv/bin/python demo.py
```

## 公开接口与数据结构

```python
from orient_audit import (
    TetMesh, audit, repair_orientation,      # 核心入口
    section, Plane,                          # 平面剖切入口
    AuditReport, BoundaryFace, InteriorFace, NonManifoldFace,  # 审计结果
    SectionReport, SectionPoint, SectionPolygon,              # 剖切结果
    CoplanarFacePatch, ContactPoint, ContactEdge,
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

## 平面剖切：`section(mesh, plane) -> SectionReport`（只读）

`Plane.build(normal, offset=0)` 校验平面（结构失败不产生部分结果）：

- `normal`：3 个**有理**分量（`int`/`Fraction`/有理字符串），拒绝
  `float`、`bool`，错误定位 `normal[k]`；**零向量拒绝**（位置 `normal`）；
- `offset`：有理标量（默认 0），拒绝 `float`，错误定位 `offset`。

`section` 的语义（全部判定走 `Fraction`，无容差）：

1. **负向输入单元先按原修正规则处理**：内部调用与
   `repair_orientation` 相同的逐单元规则（只交换最后两个顶点）；存在
   零体积单元时抛 `ZeroVolumeError`，不产生部分结果。**原网格对象不被
   修改，也不执行边分裂**——交点全部是派生出的截面点。
2. **顶点统一分类**：每个原顶点相对平面的符号
   （`n·v−d` 的正/负/零）全局只算一次。
3. **交点身份**（绝不靠浮点坐标焊接）：
   - 平面上的原顶点：`("vertex", i)`；
   - 原棱 `(i,j)` 与平面的交点：`("edge", (min(i,j), max(i,j)))`，
     坐标为精确凸插值 `(s_b·v_i − s_i·v_j)/(s_b−s_i)`（`s_k=n·v_k−d`）。
     任意多个相邻单元共享同一条棱时得到的是**同一个点对象**。
4. **逐单元组合分类**（对四面体四个三角面分别求与平面的交线段，拼成
   截面边图），每个非退化单元给出：
   - 1 个正侧 / 3 个另一侧（或反之）：**三角形片**；
   - 2 个正侧 / 2 个负侧：**四边形片**，按全局交点 ID 最小顶点与其
     环上对点的对角线**确定性**拆成两个三角形（输入排列改变不影响结果）；
   - 三个顶点在平面上：该**共面原面**登记到全局表；
   - 只触到一个顶点或一条棱：列入零面积结果，**不混进面片**；
   - 四顶点严格同侧：跳过。
5. **共面原面全局去重**：同一无序三点键只输出一份
   `CoplanarFacePatch`，`tet_indices` 记录引用它的全部来源单元
   （例如平面沿两单元的共享内面时，一份片记录两个单元）。
6. **统一朝向**：所有片的绕序满足 `n·(q−p)×(r−p) > 0`（含共面原面）。
   平面写作 `(−n,−d)` 时几何不变、绕序相反。

`SectionReport` 字段：

| 字段 | 含义 |
|---|---|
| `plane` | 输入平面 |
| `points` | 本截面用到的全部全局点（`SectionPoint(id, coord)`），按全局 ID 排序，相同身份唯一 |
| `polygons` | 每个非退化单元的三角/四边形片 `SectionPolygon` |
| `coplanar_faces` | 共面原面（每三点键一份），含全部来源单元 |
| `contact_points` | 零面积点接触 `ContactPoint(id, tet_indices)` |
| `contact_edges` | 零面积棱接触 `ContactEdge(endpoints, tet_indices)` |
| `total_area_coefficient` | 全部面积片的系数和 α（`Fraction`） |

`SectionPolygon` 含 `tet_index`、`vertices`（3 或 4 个全局 id 的绕序）、
`triangles`（三角形 1 个、四边形 2 个）、`area_coefficient` 与 `is_quad`。
**真实面积 = `area_coefficient · ‖n‖`**，其中系数
`α = n·(q−p)×(r−p)/(2 n·n)` 恒为非负有理数；`n` 为单位法向时 α 即面积。
`report.total_area_float()` 是唯一浮点出口（一次分数舍入 + 一次 `sqrt`）。

### section 错误状态

| 情形 | 行为 |
|---|---|
| `normal` 不是 3 个有理分量 / 含 float、bool | `MeshInputError`，定位 `normal[k]` |
| `normal` 为零向量 | `MeshInputError`，定位 `normal` |
| `offset` 非法（float/bool/非有理字符串） | `MeshInputError`，定位 `offset` |
| 网格含零体积单元 | `ZeroVolumeError`，`tet_indices` 列出全部，无部分结果 |
| 首参不是 `TetMesh` / 次参不是 `Plane` | `TypeError` |

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

平面剖切同样全部精确：交点坐标、平面方程、重心包含、四边形对角线
均为 `Fraction` 事实。唯一浮点出口是
`SectionReport.total_area_float()`（`float(α)·sqrt(float(n·n))`，
一次分数舍入 + 一次正确舍入 `sqrt`，误差为个位数 ULP），仅用于打印
估算，精确面积应使用 `total_area_coefficient · ‖n‖`。

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

### 平面剖切的独立核对（新增）

- 交点全部**精确满足** `n·x = d`；片内点（三角形重心）用
  **4×4 Bareiss 体积比的重心坐标**独立验证位于来源单元闭单体内
  （`tests/reference.py`，不 import 本库）。
- 面积两条独立路径：多边形**边界叉积和**公式
  `α = |n·Σ p_i×p_{i+1}|/(2 n·n)`（区别于核心的扇形三角剖分和），
  以及**解析棱锥** `V = 底面积·高/3 ⇒ α = |det6|/(2|n·apex−d|)`。
- 四个题设场景均有用例：两个共享面四面体（含平面同时切两单元、共享
  棱交点同一全局身份）、穿顶点（点/棱零面积接触单列）、平面沿共享
  内面（共面原面只一份且记录两来源单元）、穿四边形截面（确定对角线，
  单位正方形总 α=1）；另覆盖 `(0,0,+,-)` 三角截面、反向法向、负向
  输入单元的原规则修正、零体积拒绝剖切、原网格不被修改/不边分裂。

## 已知限制（本题有限协议，不宣称行业标准兼容）

- 坐标必须有理；不支持浮点网格、无理数或度量单位换算。
- 只处理四面体（4 顶点）单元；不支持混合单元、棱柱或多面体。
- 非流形的判定只按“同一三点面被 >2 个单元引用”；不检测悬点、
  自相交、孔洞规范、网格标签/材料属性等更广义的网格质量问题。
- 修复只做逐单元体积极性，不做全局一致性传播（例如彼此孤立、无法靠
  几何符号定向的零体积簇保持可报告但不可修复——这正是契约要求）。
- 规模上限 500 个四面体（接口边界，非性能极限）。
- 剖切只处理单个平面；不做多平面构造实体几何（CSG）、不输出新的持久
  网格、不执行边分裂（截面点只在 `SectionReport` 中按原顶点对身份
  派生存在）。共面判定就是精确的 `n·v−d = 0`，不提供容差模式。

## 目录

```
orient_audit/
  __init__.py   # 公开接口
  errors.py     # MeshInputError / DuplicateTetrahedraError / ZeroVolumeError
  exact.py      # 有理坐标解析、精确 3x3 行列式、输入上限
  mesh.py       # TetMesh、无序三点面关联
  audit.py      # audit / repair_orientation / AuditReport
  section.py    # Plane / section / SectionReport 及面片、接触结果类型
tests/
  fixtures.py   # 确定性立方体剖分数据
  reference.py  # 独立参考（4x4 Bareiss 行列式、面通量、外法向检查；
                # 剖切：平面方程、重心包含、边界叉积面积、解析棱锥）
  test_audit.py   # 定向审计必要测试
  test_section.py # 平面剖切必要测试
demo.py         # 固定输入演示（A 审计 / B 缺陷边界 / C 剖切四场景）
```
