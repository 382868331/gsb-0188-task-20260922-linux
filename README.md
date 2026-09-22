# tetmesh_audit — 四面体网格的定向与边界审计

离线、可复用的四面体网格质量审计库：精确判定单元定向、修复负向单元、
建立面关联、检测非流形与定向冲突、输出定向边界三角面、连通分量、
总有向体积与缺陷位置。**Python 3.14.7 标准库实现，无任何第三方依赖。**

## 运行环境与命令

- Python：**3.14.7**（仅使用 `fractions`、`math`、`typing` 等标准库；
  安装 Python 后无需安装其他依赖）。
- 全部计算离线进行：不需要网络、账号、密钥、数据库、中间件、云 API 或
  Docker。示例数据由仓库内 `tetmesh_audit/samples.py` 的确定性代码生成。

```bash
# 必要测试（39 个用例）
python -m unittest discover -s tests -v

# 固定输入演示（实测约 0.03 秒，远小于 8 秒预算；无 sleep、无压力测试）
python demo.py
```

## 公开接口

```python
from tetmesh_audit import TetMesh, audit, repair_orientation
from tetmesh_audit import (
    InvalidMeshError, DuplicateTetrahedronError, DegenerateTetrahedronError,
)
```

### 1. 数据结构 `TetMesh`

```python
TetMesh(vertices, tetrahedra)
```

- `vertices`：顶点序列，每点恰 3 个**有理坐标**。接受 `int`、
  `fractions.Fraction`、有限 `float`（按其二进制值精确转换为 Fraction，
  不做十进制舍入）或十进制字符串（如 `"0.5"`、`"1/3"`）。解析后可通过
  `mesh.vertices` 取到 `tuple[Fraction, Fraction, Fraction]`。
- `tetrahedra`：单元序列，每单元恰 4 个 `int` 顶点下标，四点互不相同。
- 构造后对象只读；`audit` / `repair_orientation` 均为纯函数，不修改入参。

### 2. `audit(mesh) -> AuditResult`

字段：

| 字段 | 含义 |
|---|---|
| `six_volumes` | 各单元 6 倍有向体积（精确 `int`/`Fraction` 的元组） |
| `total_signed_six_volume` | 非零体积单元之和（精确） |
| `total_signed_volume` | 总有向体积（`Fraction`，= 前者 / 6） |
| `boundary_det_sum` | 边界定向三角面 `det(p,q,r)` 之和 |
| `negative_tets` / `zero_volume_tets` | 负向 / 零体积单元下标 |
| `issues` | 逐单元结构化问题（`Issue(kind="negative"/"zero_volume")`） |
| `boundary_triangles` | `BoundaryTriangle(vertices, tet_index, local_face, component, orientation_valid)`，按几何外法向定向 |
| `interior_faces` | `InteriorFace(vertices, tet_indices, consistent, induced)`，恰两单元共享 |
| `defects` | `Defect(kind, face, tet_indices, detail)`，精确定位 |
| `components` | 共享面连通的四面体分量（含各分量有向体积与边界面数） |
| `is_consistent` | 无负向、无零体积、无任何缺陷 |

缺陷类型（`Defect.kind`）：

- `non_manifold_face`：同一无序三点面被 **>2** 个单元引用；
- `inconsistent_interior_face`：内面恰两单元，但边界链诱导方向相同
  （两第四顶点在面的同一侧，单纯翻转符号无法形成一致定向）；
- `zero_volume_tet`：四点共面，体积精确为零，定向无定义。

### 3. `repair_orientation(mesh) -> RepairResult`

- 精确行列式判号；**负向单元只交换其最后两个顶点（位置 2、3）**；
- **零体积单元拒绝修复**：只要存在一个即整体抛
  `DegenerateTetrahedronError`（先检查后构造，原子失败，无部分变更，
  原网格不被修改）；
- **顶点永不移动、不增删**；修复结果是一个新 `TetMesh`
  （`RepairResult.repaired_mesh`），并附 `flipped_tets`、修复前后审计；
- 非流形等结构缺陷**不阻止**符号修复，但仍保留在修复后审计中——
  本题修复规则只处理单元符号，不做全局定向传播或拓扑改动。

### 4. 错误语义

| 异常 | 触发条件 |
|---|---|
| `InvalidMeshError` | 坐标为 NaN/Infinity 或非数值；下标非 `int`（**`bool` 被拒绝**，即使它是 `int` 子类）、越界、重复；单元非 4 点；顶点数/单元数为 0 或超上限 |
| `DuplicateTetrahedronError` | 两个单元顶点集合相同（与排列无关），在网格构造阶段单独报错；`.duplicates` 给出 `(先现下标, 重复下标, 四点键)` |
| `DegenerateTetrahedronError` | 修复时存在零体积单元；`.zero_volume_tets` 给出全部下标 |

`InvalidMeshError` 带 `path` 属性与 `"(at vertices[i][j])"` 形式的可定位
信息。非法输入一律报错，绝不静默纠正。

### 输入上限（接口边界）

- 四面体数：**≤ 500**（题面契约），第 501 个在 `tetrahedra[500]` 处报错；
- 顶点数：≤ 100 000（防御性上限）；每点 3 坐标、每单元 4 下标。

## 数学口径与误差说明（重要）

- **体积判号与面关联零容差**：所有坐标解析为 `fractions.Fraction`，
  6 倍有向体积 `det(b-a, c-a, d-a)` 为精确整数/分数。`0` 就是严格的零，
  不存在“接近共面”被误判；例如第四点高度为 `10⁻⁹` 时体积精确为
  `1/10⁹ · 底面积`，不会被任何容差吞掉。因此**不用宽容差掩盖边界错误**。
- 浮点输入仅在入口按 IEEE-754 二进制值**精确**转为分数；若需要十进制
  语义，请直接传 `Fraction("0.1")` 或字符串 `"0.1"`。
- 所有输出体积为 `Fraction`；需要展示为浮点时由调用方自行
  `float(...)`（常规二进制舍入，约 15–17 位有效数字），库不以浮点值
  做任何判定。
- 面通量恒等式：对**一致定向的闭合流形**（允许多个互不接触分量），
  `boundary_det_sum == total_signed_six_volume` 精确成立
  （散度定理用于 `F=(x,y,z)/3`：`V = ⅟₆ Σ det(p,q,r)`）。
  存在翻转单元、非流形面或零体积单元时该式**可以不成立**，这本身就是
  审计信号，库不会把差异“凑”平。

## 设计取舍

- **两种面方向严格区分**：边界链诱导方向 `(-1)^i·删除序`（仅取决于单元
  顶点排列，用于判断两单元定向是否相反）与几何外法向方向（用于输出
  边界三角面，负体积单元取反）。判“相反”使用相对排序面键的排列奇偶性，
  而非三元组字面值比较（循环重排是同一方向）。
- **检测是穷举枚举**：每个面的关联单元与诱导方向都被完整列出，不存在
  搜索过程，因此没有“预算耗尽”这一状态；非流形、定向冲突、零体积都
  是被**证明存在**的缺陷，而非超时放弃。
- 修复采用题面指定的最小局部规则（仅交换最后两顶点），不实现全局
  定向传播；因此“同侧重叠”这类拓扑会如实保留为
  `inconsistent_interior_face`。
- 连通分量按“共享任意面（含非流形面）即相连”的并查集计算。
- 重复单元在构造期拒绝，审计与修复永远看不到重复数据。

## 已知限制（本题数学子集，非行业标准实现）

- 仅支持 3D 四面体单元与 0..3 局部面编号约定，不读取 Abaqus/ANSYS/Gmsh
  等网格文件格式，不做网格生成、重剖分、光滑化或全局重定向求解。
- 不检测自相交曲面、不保证修复后网格“可定向全局一致”（局部符号规则
  对莫比乌斯式/同侧重叠拓扑无能为力，此类情况以缺陷形式报告）。
- 重复单元按顶点集合定义；坐标相同但下标不同的顶点视为不同点
  （不做顶点焊接），这是精确语义下的有意选择。
- 有理坐标数值很大时分数字长会增长；500 单元规模下无实际压力。

## 目录

```
tetmesh_audit/        本题库（errors/mesh/geometry/audit/samples）
tests/                unittest 用例 + 独立参考实现 tests/reference.py
demo.py               固定输入演示（正常、修复、非流形、两类拒绝）
README.md
```

`tests/reference.py` 是**独立参考**：不 import 被测核心，期望值由
Leibniz 排列展开行列式、散度定理通量公式与立方体解析体积 `s³` 直接
算出；“审计值与通量值相等”只在其作为数学定理成立的场景用作补充核对。
