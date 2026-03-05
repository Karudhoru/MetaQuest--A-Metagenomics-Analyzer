# Changelog

All notable changes to MetaQuest are documented here.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [5.0.0] — 2026-03-03

### Added
- **`reporting/validation_engine.py`** — Multi-evidence pathogen validation engine. Vectorized gene-to-species linkage (pandas merge, 5–10× faster). Bracken confidence scoring: assigns HIGH/MODERATE/LOW confidence to species based on ratio of direct Kraken hits vs. inferred reads. Reduces risk scores by 30% for low-confidence dominant taxa.
- **`io/data_loaders.py`** — Centralized loaders for annotation files, ML predictions, and protein sequences. Streaming FASTA loader pulls only sequences with database hits (~99% memory reduction).
- **`io/text_parsers.py`** — Shared text parsing utilities (`extract_organism_name`, Prokka GFF parsers) extracted from multiple modules.
- **`core/assembly.py`** — Dedicated assembly module wrapping MEGAHIT/SPAdes with contig filtering, coverage extraction, and quality metrics.
- **`config/pathogen_config.json`** — Pathogen configuration externalized from `config.py` into a versioned config directory.

### Changed
- **`core/comparative_analysis.py`**
  - Rarefaction deprecated: warns with McMurdie & Holmes (2014) citation, falls back to TSS
  - CLR validation added: rejects data with >95% sparsity or negative values
  - Power analysis corrected using Noether (1987) nonparametric formula
  - Reproducibility: `random_seed` parameter, saves `analysis_metadata.json`
  - Memory: explicit `del` + cleanup of intermediate DataFrames
  - OutputFormatter integration replaces bare `print()` calls
  - Raw counts preserved separately for alpha diversity (which requires counts, not normalized values)
- **`reporting/risk_scoring.py`**
  - Scientific model: z-score baseline comparison replaced with absolute functional density scoring (sample-agnostic)
  - `DENSITY_THRESHOLDS`: virulence 2%, AMR 3%, transposase 5% of proteome
  - Fixed divide-by-zero when `total_cds` is `None`, `0`, or non-numeric
  - `_analyze_pathogen_db_hits()` vectorized: pandas masks replace `iterrows()` (10–50× faster)
- **`reporting/taxonomy_reporter.py`**
  - Classification rate now derived from Kraken2 `unclassified` row (Bracken estimates always showed ~100%)
  - BSL lookup reads from `pathogen_config.json` (was hardcoded dict)
  - Confidence column (HIGH/MODERATE/LOW) added using `ValidationEngine` data
  - Non-bacterial reads section added with accurate fractions
- **`reporting/functional_reporter.py`**
  - Annotation counts clarified: DIAMOND reports multiple hits per protein
  - Section 2.1 now shows "Annotated CDS" (unique) vs "Total annotations"
  - Added `_generate_methodology_note()` explaining DIAMOND multi-match behavior
- **`reporting/pathogen_risk_reporter.py`**
  - Section 5 interpretation made sample-agnostic (valid for clinical, environmental, food samples)
- **`reporting/base_reporter.py`**
  - `VERSION` constant corrected: `"5.0.0"` (was `"5.0.1"`, embedded in all report footers)
- **`core/pathogen_analysis.py`**
  - 7-factor fragment-aware confidence scoring replaces simple identity/coverage filter
  - Factors: fragment completeness (sigmoid), identity×coverage, organism abundance (log₁₀), contig depth, critical motifs, Prokka cross-validation, signal peptide penalty
  - Tiers: HIGH ≥70, MEDIUM 45–69, LOW <45
  - Outputs: `pathogen_detections_validated.json` (full) + `pathogen_results.tsv` (legacy)
- **`cli.py`** — Version display reads from `config.__version__` dynamically; `--normalization` help updated to note rarefaction deprecation
- **`config.py`** — `__version__ = "5.0.0"`, `__release_date__ = "2026-03-03"`

### Fixed
- `setup.py`: `check_external_tools()` moved inside `if __name__ == "__main__"` guard (was running at import time)
- `setup.py`: `package_data` path corrected to `config/pathogen_config.json`
- `README.md`: All Table of Contents links fixed (were pointing to `google.com/search` URLs)
- `ml/feature_extractor.py`: Removed `✅` symbols from `print()` statements being emitted to stdout

### Removed
- Development comment markers (`✅`) cleaned from all source files

---

## [4.1.0] — 2025-11-01

### Added
- **Comparative Analysis v2.0**: Full compositional data normalization (TSS, CLR, rarefaction)
- **Statistical Power Assessment**: Cohort size analysis with detectable effect size estimates
- **Multi-group Alpha Diversity**: Kruskal-Wallis with automated pairwise post-hoc tests
- **Effect Size Metrics**: Cliff's Delta (differential abundance), epsilon-squared (alpha diversity), R² (PERMANOVA)
- **ML Biomarker Discovery**: SelectKBest feature selection + GridSearchCV hyperparameter tuning
- **Prevalence Filtering**: Filters rare species before differential abundance testing

### Changed
- **Pathogen DB v4.1**: Rebuilt to filter housekeeping genes (*rpoB*, *gyrA*, etc.)
  - Sensitivity: 33.3% → 100% | Specificity: 33.3% → 100%
  - Database reduced from 70,752 to 21,152 sequences (70% smaller, 2.6× faster)
- PERMANOVA results now include R² effect size
- ANOSIM results include group separability interpretation

---

## [4.0.0] — 2025-10-01

### Added
- Enhanced COG + SwissProt dual-database annotation
- ML model v2.1 with updated training data
- Advanced gene prediction controls (contig filtering, tbl2asn timeout)
- Professional logging system (`OutputFormatter`) with `--debug` mode

---

## [3.6.0] — 2025-10-01

### Added
- SPAdes metagenomic assembly support
- Modular visualization system with specialized plotters (`taxonomic_visualizer`, `functional_visualizer`, `pathogenic_visualizer`)

---

## [3.5.0] — 2025-08-01

### Added
- Enhanced statistical testing for comparative analysis
- Advanced ML biomarker discovery pipeline
- Comprehensive comparative visualizations

---

## [3.3.0] — 2025-08-01

### Added
- Comparative analysis pipeline (initial)
- Beta diversity analysis and PCoA visualization
- Alpha diversity visualization
