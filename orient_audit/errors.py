"""错误类型与输入校验语义。

所有非法输入都抛出带有 ``location`` 的异常，调用方可以据此定位到
具体是哪个顶点 / 单元 / 下标出了问题，不做静默纠正。
"""

from __future__ import annotations


class MeshInputError(ValueError):
    """网格数据本身非法（结构问题，而非几何缺陷）。

    属性:
        location: 形如 ``"tets[3][0]"`` 的定位字符串，可能为 ``None``。
    """

    def __init__(self, message: str, location: str | None = None) -> None:
        self.location = location
        if location:
            full = f"{location}: {message}"
        else:
            full = message
        super().__init__(full)


class DuplicateTetrahedraError(MeshInputError):
    """存在重复四面体（四个顶点集合相同，与顶点排列顺序无关）。

    属性:
        duplicates: 每项为 ``(顶点下标排序元组, (首次位置, 重复位置, ...))``，
                    列出全部重复组，而不是只报第一个。
    """

    def __init__(
        self,
        duplicates: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
    ) -> None:
        self.duplicates = tuple(
            (tuple(sorted(verts)), tuple(positions))
            for verts, positions in duplicates
        )
        lines = ["检测到重复四面体（顶点集合相同即视为同一单元）:"]
        for verts, positions in self.duplicates:
            lines.append(f"  单元 {positions} 顶点集合均为 {verts}")
        super().__init__("\n".join(lines))


class EdgeSplitError(ValueError):
    """内部边分裂被拒绝或分裂后校验失败。

    属性:
        reason: 机器可读的原因串，取值为

            - ``"edge_not_found"``：网格中不存在该边；
            - ``"edge_on_boundary"``：边关联边界面（不是内部边）；
            - ``"edge_on_protected_face"``：边触碰保护面；
            - ``"non_manifold_edge"``：星内某侧面被引用次数 != 2；
            - ``"ring_disconnected"``：星形邻接环不连通（非单一闭环）；
            - ``"zero_volume_star"``：星内存在零体积单元；
            - ``"point_collides"``：插入点与现有顶点重合；
            - ``"verification_failed"``：分裂后自检（体积/朝向/邻接/
              边界/保护面守恒）未通过。

        edge: 引发拒绝的边 ``(a, b)``（按调用方给定顺序）。
        faces: 与原因相关的面键（仅 ``edge_on_protected_face`` 时给出）。
    """

    def __init__(
        self,
        reason: str,
        message: str,
        edge: tuple[int, int] | None = None,
        faces: tuple[tuple[int, int, int], ...] = (),
    ) -> None:
        self.reason = reason
        self.edge = edge
        self.faces = tuple(faces)
        super().__init__(f"[{reason}] {message}")


class ZeroVolumeError(ValueError):
    """要求修复方向时存在零体积（共面 / 退化）四面体，按契约拒绝修复。

    属性:
        tet_indices: 全部零体积单元的下标（元组，按升序）。
    """

    def __init__(self, tet_indices: tuple[int, ...]) -> None:
        self.tet_indices = tuple(sorted(tet_indices))
        if len(self.tet_indices) == 1:
            detail = f"四面体 tets[{self.tet_indices[0]}]"
        else:
            detail = "四面体 " + ", ".join(f"tets[{i}]" for i in self.tet_indices)
        super().__init__(
            f"{detail} 体积精确为零，无法通过交换顶点修复，已拒绝（顶点未被移动）"
        )
