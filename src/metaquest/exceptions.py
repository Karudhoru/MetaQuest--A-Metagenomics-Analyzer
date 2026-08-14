"""
MetaQuest Exception Hierarchy
"""


class MetaQuestError(Exception):
    """Base exception for all MetaQuest errors."""

    def __init__(self, message: str, suggestions: list[str] | None = None):
        self.suggestions = suggestions or []
        super().__init__(message)


class ConfigError(MetaQuestError):
    """Configuration loading or validation error."""


class DatabaseNotFoundError(ConfigError):
    """Required database file not found."""

    def __init__(self, db_name: str, path: str):
        super().__init__(
            f"Database '{db_name}' not found at: {path}",
            suggestions=[
                "Run 'metaquest setup-db' to download and configure databases",
                "Set METAQUEST_DB_DIR environment variable to your database directory",
                "Check databases.base_dir in your config YAML",
            ],
        )


class InputValidationError(MetaQuestError):
    """Invalid input file or parameters."""


class PipelineStageError(MetaQuestError):
    """Error during a specific pipeline stage."""

    def __init__(self, stage: str, message: str, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        full_msg = f"[{stage}] {message}"
        if cause:
            full_msg += f" (caused by: {cause})"
        super().__init__(full_msg)


class ClassificationError(PipelineStageError):
    """Error during taxonomic classification (Kraken2/Bracken)."""

    def __init__(self, message: str, cause: Exception | None = None):
        super().__init__("classification", message, cause)


class AssemblyError(PipelineStageError):
    """Error during metagenomic assembly."""

    def __init__(self, message: str, cause: Exception | None = None):
        super().__init__("assembly", message, cause)


class AnnotationError(PipelineStageError):
    """Error during gene prediction or functional annotation."""

    def __init__(self, message: str, cause: Exception | None = None):
        super().__init__("annotation", message, cause)


class ReportingError(MetaQuestError):
    """Error during report generation."""
