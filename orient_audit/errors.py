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


class TraversalError(ValueError):
    """网格不满足线段穿越采样的前置条件，拒绝输出穿越序列。

    线段穿越依赖逐单元跨面行走，要求网格几何上可行走：无零体积单元、
    无负定向单元（需先 repair_orientation 修复）、无非流形面。
    无效网格直接拒绝，不输出看起来完整实则不可信的穿越序列。

    属性:
        problems: 问题描述字符串元组（可定位，如负向/零体积单元下标、
                  非流形面的顶点键）。
    """

    def __init__(self, problems: tuple[str, ...]) -> None:
        self.problems = tuple(problems)
        lines = ["网格不满足线段穿越条件，已拒绝（不产生穿越序列）:"]
        lines.extend(f"  {p}" for p in self.problems)
        super().__init__("\n".join(lines))
