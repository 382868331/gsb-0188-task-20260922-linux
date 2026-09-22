"""精确有理运算层。

顶点坐标接受有理数：``int`` 或 ``fractions.Fraction``（含会精确转成
Fraction 的 ``str``，如 ``"1/3"``）。拒绝 ``float``，从源头避免浮点
误差；坐标与行列式全部走 :class:`fractions.Fraction`，因此体积极性
与零体积判定是精确的。
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from .errors import MeshInputError

from collections.abc import Iterable

#: 接口上限：四面体数量。
MAX_TETS = 500
#: 接口上限：单个顶点坐标分量个数。
DIMS = 3
#: 接口上限：四面体顶点数。
TET_SIZE = 4

# 有理标量的允许类型（bool 必须显式排除）。
_Rational = int | Fraction


def parse_coord(value: Any, location: str) -> Fraction:
    """把单个坐标分量解析成精确的 :class:`Fraction`。

    接受 ``int``（但拒绝 ``bool``）和 :class:`Fraction`，以及能被
    Fraction 精确解释的 ``str``（如 ``"1"``、``"-3/7"``）。
    拒绝 ``float`` 与非有限数——NaN/Infinity 在 float 入口即被拒。
    """

    if isinstance(value, bool) or not isinstance(value, int | Fraction | str):
        raise MeshInputError(
            "坐标分量必须是 int / fractions.Fraction / 有理字符串，"
            f"得到 {type(value).__name__}（不接受 float，请用 Fraction）",
            location,
        )
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    try:
        return Fraction(value)
    except (ZeroDivisionError, ValueError, ArithmeticError):
        raise MeshInputError(
            f"坐标分量 {value!r} 不是有限有理数（拒绝 NaN/Infinity）",
            location,
        )


def parse_vertex(raw: Any, index: int) -> tuple[Fraction, Fraction, Fraction]:
    """解析并校验一个三维顶点。"""

    loc = f"verts[{index}]"
    if isinstance(raw, str) or not isinstance(raw, Iterable):
        raise MeshInputError("顶点必须是长度为 3 的坐标序列", loc)
    coords = list(raw)
    if len(coords) != DIMS:
        raise MeshInputError(
            f"顶点必须恰好有 3 个坐标分量，得到 {len(coords)} 个", loc
        )
    return tuple(  # type: ignore[return-value]
        parse_coord(c, f"{loc}[{k}]") for k, c in enumerate(coords)
    )


def parse_index(value: Any, location: str, n_verts: int) -> int:
    """校验一个顶点下标：必须是真正的 int（拒绝 bool）且在范围内。"""

    if isinstance(value, bool) or not isinstance(value, int):
        raise MeshInputError(
            f"顶点下标必须是 int，得到 {type(value).__name__}", location
        )
    if not 0 <= value < n_verts:
        raise MeshInputError(
            f"顶点下标 {value} 越界（当前顶点数 {n_verts}，合法区间 [0, {n_verts - 1}]）",
            location,
        )
    return value


def signed_volume6(
    a: tuple[Fraction, ...],
    b: tuple[Fraction, ...],
    c: tuple[Fraction, ...],
    d: tuple[Fraction, ...],
) -> Fraction:
    """精确返回 6 倍有向体积 det(b-a, c-a, d-a)。

    用分数运算展开，不经过任何浮点。
    """

    u0, u1, u2 = (b[k] - a[k] for k in range(3))
    v0, v1, v2 = (c[k] - a[k] for k in range(3))
    w0, w1, w2 = (d[k] - a[k] for k in range(3))
    return (
        u0 * (v1 * w2 - v2 * w1)
        - u1 * (v0 * w2 - v2 * w0)
        + u2 * (v0 * w1 - v1 * w0)
    )
