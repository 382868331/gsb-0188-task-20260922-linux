"""orient_audit —— 四面体网格的定向与边界审计（离线、精确有理）。

公开接口：

- :class:`TetMesh` -- 已校验的不可变四面体网格（``TetMesh.build`` 构造）；
- :func:`audit` -- 定向与边界审计，返回 :class:`AuditReport`；
- :func:`repair_orientation` -- 修复负体积单元定向（只交换最后两顶点）；
- 报告类型：:class:`AuditReport`、:class:`BoundaryFace`、
  :class:`InteriorFace`、:class:`NonManifoldFace`；
- 线段穿越：:func:`traverse_segment` 及结果类型 :class:`SegmentTraversal`、
  :class:`SegmentPiece`、:class:`SegmentGap`；
- 异常：:class:`MeshInputError`、:class:`DuplicateTetrahedraError`、
  :class:`ZeroVolumeError`、:class:`InvalidTraversalMeshError`；
- 接口常量：:data:`MAX_TETS`（500）。
"""

from __future__ import annotations

from .audit import (
    AuditReport,
    BoundaryFace,
    InteriorFace,
    NonManifoldFace,
    audit,
    repair_orientation,
)
from .errors import (
    DuplicateTetrahedraError,
    InvalidTraversalMeshError,
    MeshInputError,
    ZeroVolumeError,
)
from .exact import MAX_TETS
from .mesh import FaceIncidence, FaceOccurrence, TetMesh
from .traverse import (
    SegmentGap,
    SegmentPiece,
    SegmentTraversal,
    traverse_segment,
)

__all__ = [
    "TetMesh",
    "FaceIncidence",
    "FaceOccurrence",
    "audit",
    "repair_orientation",
    "AuditReport",
    "BoundaryFace",
    "InteriorFace",
    "NonManifoldFace",
    "MeshInputError",
    "DuplicateTetrahedraError",
    "ZeroVolumeError",
    "InvalidTraversalMeshError",
    "traverse_segment",
    "SegmentTraversal",
    "SegmentPiece",
    "SegmentGap",
    "MAX_TETS",
]
