"""网格数据结构与严格输入校验。"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Iterable, List, Sequence, Tuple, Union

from .errors import InvalidMeshError

MAX_VERTICES = 100_000
"""顶点数上限（本题规模远低于此，作为接口边界存在）。"""

MAX_TETRAHEDRA = 500
"""四面体数量上限（题面契约）。"""

# 接受的有理坐标输入类型：int / Fraction / 有限 float / 十进制字符串。
RationalInput = Union[int, float, Fraction, str]
Vertex = Tuple[Fraction, Fraction, Fraction]
"""解析后的顶点：三个 :class:`~fractions.Fraction`。"""


def _parse_coordinate(value: object, path: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, RationalInput):
        raise InvalidMeshError(
            f"coordinate must be int, float, Fraction or decimal str, got {type(value).__name__}",
            path,
        )
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            raise InvalidMeshError("coordinate must not be NaN or Infinity", path)
        return Fraction(value)
    if isinstance(value, str):
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise InvalidMeshError(f"invalid rational coordinate {value!r}", path) from exc
    return Fraction(value)  # int / Fraction


class TetMesh:
    """不可变的四面体网格：有理顶点 + 四点索引单元。

    构造时执行完整校验（失败抛 :class:`InvalidMeshError`，错误信息带可定位
    路径），因此构造成功即保证：

    * 顶点数 1..MAX_VERTICES，坐标全部有限且被精确转为 Fraction；
    * 单元数 1..MAX_TETRAHEDRA，每个单元恰 4 个非 bool 整数下标，
      互不相同且在顶点范围内；
    * 不存在顶点集合相同的重复单元（否则抛
      :class:`~tetmesh_audit.errors.DuplicateTetrahedronError`）。
    """

    __slots__ = ("_vertices", "_tetrahedra", "_fraction_coordinates")

    def __init__(
        self,
        vertices: Iterable[Iterable[RationalInput]],
        tetrahedra: Iterable[Iterable[int]],
    ) -> None:
        verts: List[Vertex] = []
        try:
            vert_iter = iter(vertices)
        except TypeError as exc:
            raise InvalidMeshError("vertices must be an iterable of 3-tuples", "vertices") from exc
        for i, raw in enumerate(vert_iter):
            if i >= MAX_VERTICES:
                raise InvalidMeshError(
                    f"too many vertices: limit is {MAX_VERTICES}", f"vertices[{i}]"
                )
            try:
                coords = tuple(raw)
            except TypeError as exc:
                raise InvalidMeshError("vertex must be a 3-coordinate sequence", f"vertices[{i}]") from exc
            if len(coords) != 3:
                raise InvalidMeshError(
                    f"vertex must have exactly 3 coordinates, got {len(coords)}",
                    f"vertices[{i}]",
                )
            verts.append(
                tuple(
                    _parse_coordinate(c, f"vertices[{i}][{j}]") for j, c in enumerate(coords)
                )  # type: ignore[arg-type]
            )
        if not verts:
            raise InvalidMeshError("mesh must contain at least one vertex", "vertices")

        tets: List[Tuple[int, int, int, int]] = []
        seen: dict = {}
        duplicate_pairs: List[tuple] = []
        try:
            tet_iter = iter(tetrahedra)
        except TypeError as exc:
            raise InvalidMeshError("tetrahedra must be an iterable of 4-tuples", "tetrahedra") from exc
        for i, raw in enumerate(tet_iter):
            if i >= MAX_TETRAHEDRA:
                raise InvalidMeshError(
                    f"too many tetrahedra: limit is {MAX_TETRAHEDRA}", f"tetrahedra[{i}]"
                )
            try:
                ids = tuple(raw)
            except TypeError as exc:
                raise InvalidMeshError(
                    "tetrahedron must be a 4-index sequence", f"tetrahedra[{i}]"
                ) from exc
            if len(ids) != 4:
                raise InvalidMeshError(
                    f"tetrahedron must have exactly 4 vertex indices, got {len(ids)}",
                    f"tetrahedra[{i}]",
                )
            for j, idx in enumerate(ids):
                if isinstance(idx, bool) or not isinstance(idx, int):
                    raise InvalidMeshError(
                        f"vertex index must be int, got {type(idx).__name__}",
                        f"tetrahedra[{i}][{j}]",
                    )
                if not 0 <= idx < len(verts):
                    raise InvalidMeshError(
                        f"vertex index {idx} out of range [0, {len(verts) - 1}]",
                        f"tetrahedra[{i}][{j}]",
                    )
            if len(set(ids)) != 4:
                raise InvalidMeshError(
                    "tetrahedron references repeated vertices (zero-volume topology)",
                    f"tetrahedra[{i}]",
                )
            key = tuple(sorted(ids))
            if key in seen:
                duplicate_pairs.append((seen[key], i, key))
            else:
                seen[key] = i
            tets.append(ids)  # type: ignore[arg-type]
        if not tets:
            raise InvalidMeshError("mesh must contain at least one tetrahedron", "tetrahedra")
        if duplicate_pairs:
            from .errors import DuplicateTetrahedronError

            raise DuplicateTetrahedronError(duplicate_pairs)

        self._vertices = verts
        self._tetrahedra = tets

    @property
    def vertices(self) -> Tuple[Vertex, ...]:
        return tuple(self._vertices)

    @property
    def tetrahedra(self) -> Tuple[Tuple[int, int, int, int], ...]:
        return tuple(self._tetrahedra)

    @property
    def num_vertices(self) -> int:
        return len(self._vertices)

    @property
    def num_tetrahedra(self) -> int:
        return len(self._tetrahedra)

    def tet_vertices(self, i: int) -> Tuple[Vertex, Vertex, Vertex, Vertex]:
        t = self._tetrahedra[i]
        return (
            self._vertices[t[0]],
            self._vertices[t[1]],
            self._vertices[t[2]],
            self._vertices[t[3]],
        )

    def __repr__(self) -> str:  # pragma: no cover - 调试辅助
        return f"TetMesh(vertices={len(self._vertices)}, tetrahedra={len(self._tetrahedra)})"
