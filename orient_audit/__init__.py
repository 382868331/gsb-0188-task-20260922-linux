"""orient_audit —— 四面体网格的定向与边界审计（离线、精确有理）。

公开接口：

- :class:`TetMesh` -- 已校验的不可变四面体网格（``TetMesh.build`` 构造）；
- :func:`audit` -- 定向与边界审计，返回 :class:`AuditReport`；
- :func:`repair_orientation` -- 修复负体积单元定向（只交换最后两顶点）；
- :func:`section` -- 有理平面剖切，输出三角网与平面相交的有理多边形片；
- 平面/报告类型：:class:`Plane`、:class:`SectionReport`、
  :class:`SectionPoint`、:class:`SectionPolygon`、:class:`CoplanarFacePatch`、
  :class:`ContactPoint`、:class:`ContactEdge`；
- 报告类型：:class:`AuditReport`、:class:`BoundaryFace`、
  :class:`InteriorFace`、:class:`NonManifoldFace`；
- 异常：:class:`MeshInputError`、:class:`DuplicateTetrahedraError`、
  :class:`ZeroVolumeError`；
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
    ZeroVolumeError,
)
from .exact import MAX_TETS
from .mesh import FaceIncidence, FaceOccurrence, TetMesh
from .section import (
    ContactEdge,
    ContactPoint,
    CoplanarFacePatch,
    Plane,
    SectionPoint,
    SectionPolygon,
    SectionReport,
    section,
)

__all__ = [
    "TetMesh",
    "FaceIncidence",
    "FaceOccurrence",
    "audit",
    "repair_orientation",
    "section",
    "Plane",
    "SectionReport",
    "SectionPoint",
    "SectionPolygon",
    "CoplanarFacePatch",
    "ContactPoint",
    "ContactEdge",
    "AuditReport",
    "BoundaryFace",
    "InteriorFace",
    "NonManifoldFace",
    "MeshInputError",
    "DuplicateTetrahedraError",
    "ZeroVolumeError",
    "MAX_TETS",
]
