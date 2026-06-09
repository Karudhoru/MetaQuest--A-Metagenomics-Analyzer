"""
MetaQuest Pipeline Architecture
================================
Decouples orchestration (sequencing of stages, error handling, progress reporting)
from computation (the actual bioinformatics logic in core/).
"""

from .runner import PipelineRunner
from .context import PipelineContext

__all__ = ["PipelineRunner", "PipelineContext"]
