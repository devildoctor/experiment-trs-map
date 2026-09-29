"""串联重复候选 read 选择所需的核心数据模型与算法。"""

from .catalog import Locus, parse_repeat_catalog
from .evidence import AlignmentEvidence, SelectionConfig, evaluate_alignment

__all__ = [
    "AlignmentEvidence",
    "Locus",
    "SelectionConfig",
    "evaluate_alignment",
    "parse_repeat_catalog",
]
