"""orient_audit —— 四面体网格的定向与边界审计（离线、精确有理）。

公开接口：

- :class:`TetMesh` -- 已校验的不可变四面体网格（``TetMesh.build`` 构造）；
- :func:`audit` -- 定向与边界审计，返回 :class:`AuditReport`；
- :func:`repair_orientation` -- 修复负体积单元定向（只交换最后两顶点）；
- :func:`merge_meshes` -- 按调用方显式指定的重合点对拼接两份网格
  （不自动吸附），返回 :class:`MergeResult`（含新旧编号关系与复检报告）；
- 报告类型：:class:`AuditReport`、:class:`BoundaryFace`、
  :class:`InteriorFace`、:class:`NonManifoldFace`；
- 异常：:class:`MeshInputError`、:class:`DuplicateTetrahedraError`、
  :class:`ZeroVolumeError`、:class:`MergePointError`；
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
    MergePointError,
    ZeroVolumeError,
)
from .exact import MAX_TETS
from .merge import MergeResult, merge_meshes
from .mesh import FaceIncidence, FaceOccurrence, TetMesh

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
    "MergePointError",
    "ZeroVolumeError",
    "MergeResult",
    "merge_meshes",
    "MAX_TETS",
]
