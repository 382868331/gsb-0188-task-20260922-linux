"""orient_audit —— 四面体网格的定向与边界审计（离线、精确有理）。

公开接口：

- :class:`TetMesh` -- 已校验的不可变四面体网格（``TetMesh.build`` 构造）；
- :func:`audit` -- 定向与边界审计，返回 :class:`AuditReport`；
- :func:`repair_orientation` -- 修复负体积单元定向（只交换最后两顶点）；
- :func:`split_interior_edge` -- 内部边分裂（有理位置插点重建相邻四面体，
  保留材料标签与保护面，失败原子回滚），返回 :class:`EdgeSplitResult`；
- 报告类型：:class:`AuditReport`、:class:`BoundaryFace`、
  :class:`InteriorFace`、:class:`NonManifoldFace`；
- 异常：:class:`MeshInputError`、:class:`DuplicateTetrahedraError`、
  :class:`ZeroVolumeError`、:class:`EdgeSplitError`、
  :class:`ProtectedFaceError`；
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
    ProtectedFaceError,
    ZeroVolumeError,
)
from .exact import MAX_TETS
from .mesh import FaceIncidence, FaceOccurrence, TetMesh
from .split import EdgeSplitResult, split_interior_edge

__all__ = [
    "TetMesh",
    "FaceIncidence",
    "FaceOccurrence",
    "audit",
    "repair_orientation",
    "split_interior_edge",
    "EdgeSplitResult",
    "AuditReport",
    "BoundaryFace",
    "InteriorFace",
    "NonManifoldFace",
    "MeshInputError",
    "DuplicateTetrahedraError",
    "ZeroVolumeError",
    "EdgeSplitError",
    "ProtectedFaceError",
    "MAX_TETS",
]
