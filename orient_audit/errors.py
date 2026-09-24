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


class EdgeSplitError(ValueError):
    """内部边分裂在几何 / 拓扑 / 规模上无法执行。

    与 :class:`MeshInputError` 的区别：输入结构本身合法（网格已成功构造），
    但指定边不满足分裂前提（边不存在、邻接非流形、星内零体积、超规模）。
    拒绝时不产生任何新网格，调用方网格保持原状（原子语义）。

    属性:
        edge: 被请求分裂的边 ``(a, b)``（可能为 ``None``）。
        tet_indices: 涉及的单元下标（如零体积星内单元），升序元组。
        reason: 机器可读的短原因码。
    """

    def __init__(
        self,
        message: str,
        *,
        edge: tuple[int, int] | None = None,
        tet_indices: tuple[int, ...] = (),
        reason: str = "edge_split_rejected",
    ) -> None:
        self.edge = edge
        self.tet_indices = tuple(sorted(tet_indices))
        self.reason = reason
        super().__init__(message)


class ProtectedFaceError(EdgeSplitError):
    """分裂会破坏受保护面（保护边界），按契约拒绝。

    受保护面有两个来源：

    - 全部边界面（边界默认受保护：分裂边界边会把边界三角面一分为二）；
    - 构造网格时经 ``protected_faces`` 显式声明的面（如材料分界面）。

    属性:
        boundary_faces: 被冒犯的边界面三点键（升序）。
        explicit_faces: 被冒犯的显式保护面三点键（升序）。
    """

    def __init__(
        self,
        *,
        edge: tuple[int, int],
        boundary_faces: tuple[tuple[int, int, int], ...] = (),
        explicit_faces: tuple[tuple[int, int, int], ...] = (),
    ) -> None:
        self.boundary_faces = tuple(sorted(boundary_faces))
        self.explicit_faces = tuple(sorted(explicit_faces))
        offending = self.boundary_faces + self.explicit_faces
        lines = [
            f"边 {edge} 位于受保护面上，分裂会细分并破坏保护边界，已拒绝"
            "（网格未被修改）:"
        ]
        for key in self.boundary_faces:
            lines.append(f"  边界面 {key}（边界默认受保护）")
        for key in self.explicit_faces:
            lines.append(f"  显式保护面 {key}")
        super().__init__(
            "\n".join(lines),
            edge=edge,
            reason="protected_face",
        )
        # offending 仅用于便于调用方一次性取全部键。
        self.violating_faces = offending
