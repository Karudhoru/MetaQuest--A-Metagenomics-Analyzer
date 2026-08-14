"""Load and validate tabular outputs used by the stable pipeline."""

from pathlib import Path

import pandas as pd


def load_bracken_report(bracken_file: Path) -> pd.DataFrame:
    """Load a Bracken abundance report with its required stable columns."""
    bracken_file = Path(bracken_file)
    if not bracken_file.is_file():
        raise FileNotFoundError(f"Bracken report not found: {bracken_file}")

    table = pd.read_csv(bracken_file, sep="\t")
    required = {
        "name",
        "taxonomy_id",
        "new_est_reads",
        "fraction_total_reads",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(
            "Bracken report is missing required columns: "
            + ", ".join(sorted(missing))
        )
    return table
