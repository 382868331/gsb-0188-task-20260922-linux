"""tetmesh_audit — 四面体网格定向与边界审计（纯标准库）。

公开接口见模块文档字符串与 README。所有几何符号判定均使用
:mod:`fractions` 精确有理运算，不使用任何浮点容差。
"""

from .errors import (
    DegenerateTetrahedronError,
    DuplicateTetrahedronError,
    InvalidMeshError,
    Issue,
    TetMeshAuditError,
)
from .mesh import MAX_TETRAHEDRA, MAX_VERTICES, TetMesh, Vertex
from .audit import (
    AuditResult,
    BoundaryTriangle,
    Component,
    Defect,
    InteriorFace,
    RepairResult,
    audit,
    repair_orientation,
)
from . import samples

__version__ = "1.0.0"

__all__ = [
    "TetMesh",
    "Vertex",
    "audit",
    "repair_orientation",
    "AuditResult",
    "RepairResult",
    "BoundaryTriangle",
    "InteriorFace",
    "Component",
    "Defect",
    "Issue",
    "TetMeshAuditError",
    "InvalidMeshError",
    "DuplicateTetrahedronError",
    "DegenerateTetrahedronError",
    "MAX_TETRAHEDRA",
    "MAX_VERTICES",
    "samples",
    "__version__",
]
