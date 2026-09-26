"""orient_audit —— 四面体网格的定向与边界审计（离线、精确有理）。

公开接口：

- :class:`TetMesh` -- 已校验的不可变四面体网格（``TetMesh.build`` 构造）；
- :func:`audit` -- 定向与边界审计，返回 :class:`AuditReport`；
- :func:`repair_orientation` -- 修复负体积单元定向（只交换最后两顶点）；
- 报告类型：:class:`AuditReport`、:class:`BoundaryFace`、
  :class:`InteriorFace`、:class:`NonManifoldFace`；
- 异常：:class:`MeshInputError`、:class:`DuplicateTetrahedraError`、
  :class:`ZeroVolumeError`、:class:`TraversalError`；
- :func:`sample_segment` -- 线段穿越采样与积分（分段线性场），返回
  :class:`SegmentTraversal`（命中段 :class:`SegmentHit`、网格外间隙
  :class:`SegmentGap` 分开列出）；
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
    MeshInputError,
    TraversalError,
    ZeroVolumeError,
)
from .exact import MAX_TETS
from .mesh import FaceIncidence, FaceOccurrence, TetMesh
from .traverse import SegmentGap, SegmentHit, SegmentTraversal, sample_segment

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
    "TraversalError",
    "sample_segment",
    "SegmentTraversal",
    "SegmentHit",
    "SegmentGap",
    "MAX_TETS",
]
