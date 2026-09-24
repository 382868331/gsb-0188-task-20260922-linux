"""orient_audit —— 四面体网格的定向与边界审计（离线、精确有理）。

公开接口：

- :class:`TetMesh` -- 已校验的不可变四面体网格（``TetMesh.build`` 构造，
  支持逐单元材料标签 ``materials`` 与 ``protected_faces`` 保护面）；
- :func:`audit` -- 定向与边界审计，返回 :class:`AuditReport`；
- :func:`repair_orientation` -- 修复负体积单元定向（只交换最后两顶点）；
- :func:`split_internal_edge` -- 内部边有理位置插点分裂，返回
  :class:`SplitResult`（保护材料标签与保护面，失败原子拒绝）；
- 报告类型：:class:`AuditReport`、:class:`BoundaryFace`、
  :class:`InteriorFace`、:class:`NonManifoldFace`、:class:`SplitResult`；
- 异常：:class:`MeshInputError`、:class:`DuplicateTetrahedraError`、
  :class:`ZeroVolumeError`、:class:`EdgeSplitError`；
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
    EdgeSplitError,
    MeshInputError,
    ZeroVolumeError,
)
from .exact import MAX_TETS
from .mesh import FaceIncidence, FaceOccurrence, TetMesh
from .split import SplitResult, split_internal_edge

__all__ = [
    "TetMesh",
    "FaceIncidence",
    "FaceOccurrence",
    "audit",
    "repair_orientation",
    "split_internal_edge",
    "AuditReport",
    "BoundaryFace",
    "InteriorFace",
    "NonManifoldFace",
    "SplitResult",
    "MeshInputError",
    "DuplicateTetrahedraError",
    "ZeroVolumeError",
    "EdgeSplitError",
    "MAX_TETS",
]
