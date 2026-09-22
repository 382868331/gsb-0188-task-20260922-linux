"""错误类型与结构化问题记录。"""

from __future__ import annotations

from typing import List


class TetMeshAuditError(Exception):
    """本库所有异常的基类。"""


class InvalidMeshError(TetMeshAuditError):
    """输入不合法（NaN/无穷、索引越界、重复顶点等）。

    :attr path: 可定位的字段路径，例如 ``"vertices[3][1]"``、
    ``"tetrahedra[2]"``；顶层问题路径为空字符串。
    """

    def __init__(self, message: str, path: str = "") -> None:
        self.path = path
        if path:
            message = f"{message} (at {path})"
        super().__init__(message)


class DuplicateTetrahedronError(TetMeshAuditError):
    """存在重复单元（顶点集合相同，与排列方向无关）。

    :attr duplicates: 每项为 ``(先出现的单元下标, 重复单元下标, 排序后的四点索引元组)``。
    """

    def __init__(self, duplicates: List[tuple]) -> None:
        self.duplicates = duplicates
        n = len(duplicates)
        first = duplicates[0]
        detail = (
            f"tetrahedron #{first[1]} duplicates tetrahedron #{first[0]} "
            f"(vertex set {first[2]})"
        )
        if n > 1:
            detail += f"; {n - 1} more duplicate pair(s)"
        super().__init__(f"{n} duplicate tetrahedron pair(s): {detail}")


class DegenerateTetrahedronError(TetMeshAuditError):
    """存在零体积（共面退化）单元，拒绝修复。

    :attr zero_volume_tets: 零体积单元下标列表（升序、去重）。
    """

    def __init__(self, zero_volume_tets: List[int]) -> None:
        self.zero_volume_tets = sorted(set(zero_volume_tets))
        super().__init__(
            f"{len(self.zero_volume_tets)} zero-volume tetrahedron(s) at "
            f"indices {self.zero_volume_tets}; repair refused"
        )


class Issue:
    """单个单元的定向问题（结构化记录，非异常）。

    属性：
        tet_index: 单元下标。
        kind: ``"negative"``（体积为负）或 ``"zero_volume"``（零体积）。
        signed_six_volume: 6 倍有向体积（精确整数，坐标为整数时；否则为 Fraction）。
    """

    __slots__ = ("tet_index", "kind", "signed_six_volume")

    def __init__(self, tet_index: int, kind: str, signed_six_volume) -> None:
        self.tet_index = tet_index
        self.kind = kind
        self.signed_six_volume = signed_six_volume

    def __repr__(self) -> str:  # pragma: no cover - 调试辅助
        return f"Issue(tet_index={self.tet_index!r}, kind={self.kind!r})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Issue):
            return NotImplemented
        return (
            self.tet_index == other.tet_index
            and self.kind == other.kind
            and self.signed_six_volume == other.signed_six_volume
        )
