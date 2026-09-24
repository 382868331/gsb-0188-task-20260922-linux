"""四面体网格数据结构：解析、校验、面关联。

公开数据结构（见模块 docstring 与 README）：

``TetMesh(vertices, tets)``
    不可变网格。``vertices`` 为三元组序列（内部存 :class:`Fraction`），
    ``tets`` 为 4 个顶点下标组成的序列。

``FaceIncidence(key, occurrences)``
    面关联记录。``key`` 是排序后的三个顶点下标（无序三点键）；
    ``occurrences`` 是 ``(单元下标, 有序三元组)`` 列表，有序三元组为
    该单元中该面的诱导朝向（按 (v1,v2,v3),(v0,v3,v2),(v0,v1,v3),(v2,v1,v0)）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .errors import DuplicateTetrahedraError, MeshInputError
from .exact import MAX_TETS, TET_SIZE, parse_index, parse_vertex
from fractions import Fraction

#: 四面体四个面的局部顶点位（由右手定则与外法向一致的诱导顺序）。
#: 当面在单元中以此序出现时，其外法向指向单元外侧。
FACE_LOCAL: tuple[tuple[int, int, int], ...] = (
    (1, 2, 3),  # 与 v0 相对
    (0, 3, 2),  # 与 v1 相对
    (0, 1, 3),  # 与 v2 相对
    (2, 1, 0),  # 与 v3 相对
)

Vertex = tuple[Fraction, Fraction, Fraction]
#: 材料标签：任意可哈希值（默认 0）。分裂时子单元继承所属单元标签。
TetLabel = Any
#: 受保护面以无序三点键表示（排序后的顶点下标）。
FaceKey = tuple[int, int, int]


def _normalize_face_key(raw: Any, location: str, n_verts: int) -> FaceKey:
    """校验单个保护面声明：恰好 3 个互不相同且不越界的 int 下标。

    返回排序后的三点键。
    """

    if isinstance(raw, str) or not isinstance(raw, (list, tuple)):
        raise MeshInputError("保护面必须是 3 个顶点下标的序列", location)
    if len(raw) != 3:
        raise MeshInputError(
            f"保护面必须恰好含 3 个顶点下标，得到 {len(raw)} 个", location
        )
    idx = tuple(parse_index(x, f"{location}[{k}]", n_verts) for k, x in enumerate(raw))
    if len(set(idx)) != 3:
        raise MeshInputError(
            f"保护面的 3 个顶点下标必须互不相同，得到 {list(idx)}", location
        )
    return tuple(sorted(idx))  # type: ignore[return-value]


@dataclass(frozen=True)
class TetMesh:
    """已校验的四面体网格（不可变）。

    属性:
        vertices: ``tuple[Vertex, ...]``，坐标均为精确 :class:`Fraction`。
        tets: ``tuple[tuple[int, int, int, int], ...]``。
        tet_labels: 每个四面体的材料标签元组；未指定时全为 0。
                    内部边分裂产生的子单元继承原单元标签。
        protected_faces: 显式声明的受保护面三点键元组（排序三点）。
                         边界面默认始终受保护，无需在此重复声明；
                         此处用于内部材料分界面等额外保护面。
    """

    vertices: tuple[Vertex, ...]
    tets: tuple[tuple[int, int, int, int], ...]
    tet_labels: tuple[TetLabel, ...] = ()
    protected_faces: tuple[FaceKey, ...] = ()

    def __post_init__(self) -> None:
        # 兼容旧的两字段位置构造 TetMesh(vertices, tets)：标签缺省补全为 0。
        if len(self.tet_labels) != len(self.tets):
            if self.tet_labels:
                raise ValueError(
                    "tet_labels 长度必须与 tets 一致（请用 TetMesh.build 构造）"
                )
            object.__setattr__(self, "tet_labels", tuple(0 for _ in self.tets))

    @classmethod
    def build(
        cls,
        vertices,
        tets,
        *,
        tet_labels=None,
        protected_faces=None,
    ) -> "TetMesh":
        """从类序列输入构造网格，做完整结构校验。

        校验顺序固定（先顶点后面，最后附属数据），所有错误带 ``location``：

        1. 顶点数 / 坐标合法性（拒绝 bool、float、NaN/Infinity）；
        2. 单元数不超过 :data:`~orient_audit.exact.MAX_TETS`；
        3. 每个单元恰好 4 个下标，下标为 int（拒绝 bool）且不越界，
           单元内下标互不相同；
        4. 无重复单元（四点集合相同即重复，与排列无关）；
        5. ``tet_labels``（若给）长度与单元数一致，每个标签可哈希；
        6. ``protected_faces``（若给）每个为 3 个合法、互异下标，去重排序。

        任何一步失败都抛出 :class:`MeshInputError` 或其子类，
        不产生部分构造结果（本方法要么返回完整对象要么抛异常）。
        """

        # ---- 顶点 ----
        if isinstance(vertices, str) or not isinstance(vertices, (list, tuple)):
            raise MeshInputError(
                "vertices 必须是 list/tuple 形式的顶点序列", "verts"
            )
        parsed_verts = tuple(
            parse_vertex(v, k) for k, v in enumerate(vertices)
        )
        n_verts = len(parsed_verts)

        # ---- 单元列表 ----
        if isinstance(tets, str) or not isinstance(tets, (list, tuple)):
            raise MeshInputError("tets 必须是 list/tuple 形式的单元序列", "tets")
        if len(tets) > MAX_TETS:
            raise MeshInputError(
                f"四面体数量 {len(tets)} 超过上限 {MAX_TETS}", "tets"
            )
        parsed_tets: list[tuple[int, int, int, int]] = []
        for ti, raw_tet in enumerate(tets):
            loc = f"tets[{ti}]"
            if isinstance(raw_tet, str) or not isinstance(raw_tet, (list, tuple)):
                raise MeshInputError("四面体必须是 4 个顶点下标的序列", loc)
            if len(raw_tet) != TET_SIZE:
                raise MeshInputError(
                    f"四面体必须恰好含 {TET_SIZE} 个顶点下标，得到 {len(raw_tet)} 个",
                    loc,
                )
            idx = tuple(
                parse_index(x, f"{loc}[{k}]", n_verts) for k, x in enumerate(raw_tet)
            )
            if len(set(idx)) != TET_SIZE:
                raise MeshInputError(
                    f"四面体的 4 个顶点下标必须互不相同，得到 {list(idx)}", loc
                )
            parsed_tets.append(idx)  # type: ignore[arg-type]

        # ---- 重复单元：四点集合相同即重复，与排列无关 ----
        seen: dict[tuple[int, ...], list[int]] = {}
        for ti, idx in enumerate(parsed_tets):
            seen.setdefault(tuple(sorted(idx)), []).append(ti)
        dup = tuple((key, tuple(pos)) for key, pos in seen.items() if len(pos) > 1)
        if dup:
            raise DuplicateTetrahedraError(dup)

        # ---- 材料标签 ----
        n_tets = len(parsed_tets)
        if tet_labels is None:
            labels: tuple[TetLabel, ...] = tuple(0 for _ in parsed_tets)
        else:
            if isinstance(tet_labels, str) or not isinstance(
                tet_labels, (list, tuple)
            ):
                raise MeshInputError(
                    "tet_labels 必须是 list/tuple 且长度与单元数一致", "tet_labels"
                )
            if len(tet_labels) != n_tets:
                raise MeshInputError(
                    f"tet_labels 长度 {len(tet_labels)} 与单元数 {n_tets} 不一致",
                    "tet_labels",
                )
            checked: list[TetLabel] = []
            for li, label in enumerate(tet_labels):
                try:
                    hash(label)
                except TypeError:
                    raise MeshInputError(
                        f"材料标签必须可哈希，得到 {type(label).__name__}",
                        f"tet_labels[{li}]",
                    )
                checked.append(label)
            labels = tuple(checked)

        # ---- 显式保护面 ----
        if protected_faces is None:
            prot: tuple[FaceKey, ...] = ()
        else:
            if isinstance(protected_faces, str) or not isinstance(
                protected_faces, (list, tuple)
            ):
                raise MeshInputError(
                    "protected_faces 必须是 list/tuple 形式的三点面序列",
                    "protected_faces",
                )
            keys: set[FaceKey] = set()
            for fi, raw_face in enumerate(protected_faces):
                keys.add(
                    _normalize_face_key(raw_face, f"protected_faces[{fi}]", n_verts)
                )
            prot = tuple(sorted(keys))

        return cls(parsed_verts, tuple(parsed_tets), labels, prot)

    def with_repaired_tets(
        self, new_tets: tuple[tuple[int, int, int, int], ...]
    ) -> "TetMesh":
        """返回更换单元顺序后的新网格；顶点对象原样保留（修复不移动顶点）。

        材料标签与保护面随单元位置原样保留（逐单元修复不增删单元）。
        """

        return TetMesh(
            self.vertices, new_tets, self.tet_labels, self.protected_faces
        )

    def boundary_face_keys(self) -> frozenset[FaceKey]:
        """当前网格全部边界面的无序三点键（只被一个单元引用的面）。"""

        counts: dict[FaceKey, int] = {}
        for tet in self.tets:
            for local in FACE_LOCAL:
                key = tuple(sorted((tet[local[0]], tet[local[1]], tet[local[2]])))
                counts[key] = counts.get(key, 0) + 1  # type: ignore[index]
        return frozenset(key for key, n in counts.items() if n == 1)


@dataclass(frozen=True)
class FaceOccurrence:
    """面在某个单元中的一次出现。"""

    tet_index: int
    ordered: tuple[int, int, int]


@dataclass(frozen=True)
class FaceIncidence:
    """单个无序三点面的全部关联信息。"""

    key: tuple[int, int, int]
    occurrences: tuple[FaceOccurrence, ...]

    @property
    def multiplicity(self) -> int:
        return len(self.occurrences)

    @property
    def is_boundary(self) -> bool:
        """恰被一个单元引用 → 边界面。"""

        return self.multiplicity == 1

    @property
    def is_manifold_interior(self) -> bool:
        """恰被两个单元引用且诱导方向相反 → 合法内面。"""

        return self.multiplicity == 2 and self._opposite_orientation()

    def _opposite_orientation(self) -> bool:
        if self.multiplicity != 2:
            return False
        (o1, o2) = self.occurrences
        # 同一三点集，排列奇偶性相反即朝向相反。
        return _permutation_parity(o1.ordered, self.key) != _permutation_parity(
            o2.ordered, self.key
        )


def _permutation_parity(ordered: tuple[int, int, int], key: tuple[int, int, int]) -> int:
    """返回 ordered 相对排序 key 的排列奇偶性（偶=0，奇=1）。"""

    pos = {v: i for i, v in enumerate(key)}
    p = [pos[v] for v in ordered]
    inversions = 0
    for i in range(3):
        for j in range(i + 1, 3):
            if p[i] > p[j]:
                inversions += 1
    return inversions % 2


def build_face_incidence(
    mesh: TetMesh,
) -> dict[tuple[int, int, int], FaceIncidence]:
    """按无序三点键汇总每个面被哪些单元以何种诱导朝向引用。"""

    table: dict[tuple[int, int, int], list[FaceOccurrence]] = {}
    for ti, tet in enumerate(mesh.tets):
        for local in FACE_LOCAL:
            ordered = (tet[local[0]], tet[local[1]], tet[local[2]])
            key = tuple(sorted(ordered))
            table.setdefault(key, []).append(FaceOccurrence(ti, ordered))  # type: ignore[index]
    return {
        key: FaceIncidence(key, tuple(occ))
        for key, occ in sorted(table.items())
    }
