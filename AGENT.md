# AGENT.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Snakemake-based differential gene expression analysis pipeline that supports multiple quantification methods (RSEM, Salmon, Kallisto, featureCounts) and implements DESeq2/edgeR statistical analysis with comprehensive visualization.

## Architecture

**Main components:**
- `Snakefile` - Main workflow entry point, loads config, defines global variables, includes rules
- `rules/differential_expression.smk` - Core analysis rules for table creation and DE computation
- `wrappers/` - Analysis modules containing scripts and conda environments:
  - `analysis_*_table/` - Scripts to aggregate raw quantification files into count tables
  - `DE_computation/` - Shared R utility functions (DESeq2_func.R, edgeR_func.R, load_data_func.R, DESeq2_edgeR_overlap_func.R)
  - `DE_computation_loading_data/` - Load count data and create count_data_original/txi objects
  - `DE_computation_PCA/` - Normalization and PCA/heatmap/MDS visualization; prepare comparison-specific normalized data
  - `DE_computation_DESeq2/` - Run DESeq2 for specific comparisons
  - `DE_computation_edgeR/` - Run edgeR for specific comparisons
  - `DE_computation_merge_results/` - Compare DESeq2 and edgeR results
  - `DE_report/` - HTML report generation

**Data flow:**
1. Load sample configuration from `config.json` via BioRoots utilities module
2. Aggregate per-sample quantification files → complete count tables (`.tsv` or `.RData`)
3. **Modular DE computation pipeline:**
   - `load_count_data` → creates `loading_data/count_data_original.RDS` and `txi.RDS`
   - `normalize_and_visualize` → creates `PCA/` with normalized objects and visualizations for all samples
   - `prepare_comparison_data` → copies or recomputes normalization per comparison (based on `normalize_data_per_comparison` config)
   - `deseq2_computation` / `edger_computation` → extract DE results per comparison
   - `deseq2_edger_overlap` → generates overlap visualizations
4. Generate consolidated HTML report (`final_report.html`)

## Configuration

Requires `config.json` based on `config.example.json`:
- **Required fields**: `organism`, `organism_gtf`, `samples` (with `sample_name`, `condition`)
- **Analysis selection**: `featureCount`, `RSEM`, `salmon_map`, `salmon_align`, `kallisto`, `is_mirna`
- **Comparison settings**: `conditions_to_compare` ("all" or comma-separated pairs like "treated:control")
- **Normalization mode**: `normalize_data_per_comparison` (FALSE = normalize all samples together; TRUE = normalize per comparison)
- **Filtering thresholds**: `remove_genes_with_sum_read_count_threshold`, `biotypes`, `pvalue_for_viz`, `fold_change_threshold`

The workflow also reads `workflow.config.json` for GUI integration parameters.

## Running the Pipeline

```bash
# Full run with default resources
snakemake --cores <N>

# Dry run to check workflow
snakemake --cores <N> --dryrun

# Run specific analysis type only
snakemake --cores <N> DE_featureCount_exon/final_report.html

# Use custom config file location
snakemake --cores <N> --configfile /path/to/config.json
```

## Common Development Tasks

**Add new quantification method:**
1. Create wrapper directory under `wrappers/` with `script.py` and `env.yaml`
2. Add rule in `rules/differential_expression.smk` to create count table
3. Update `Snakefile` analysis detection logic (lines 55-85)
4. Update `input_all()` function to include new analysis type

**Modify DE computation:**
- Edit utility functions in `wrappers/DE_computation/` (DESeq2_func.R, edgeR_func.R, load_data_func.R)
- Modify wrapper scripts in `DE_computation_DESeq2/`, `DE_computation_edgeR/`, etc. for pipeline changes
- Parameters passed via Python wrapper scripts (script.py)

**Debug a failed rule:**
```bash
snakemake --cores 1 <failed_output_file> --rerun-incomplete --verbose
```

## Key Dependencies

- **Snakemake** >= 5.18.0
- **BioRoots utilities** (external module from `BioIT-CEITEC/bioroots_utilities`)
- **R packages**: DESeq2, edgeR, limma, ggplot2, pheatmap, ggpubr, upsetr
- **Python**: pandas

Each wrapper has its own conda environment defined in `env.yaml`.

## Output Structure

```
DE_{analysis_type}/                    # e.g., DE_featureCount_exon, DE_RSEM
│
├── loading_data/                      # Created by load_count_data
│   ├── count_data_original.RDS
│   └── txi.RDS (NULL for featureCount)
│
├── PCA/                               # Created by normalize_and_visualize
│   ├── dds_normalized.RDS
│   ├── count_data_normalized.RDS
│   ├── edgeR_DGEList_normalized.RDS
│   ├── edgeR_fit_normalized.RDS
│   ├── DESeq2/report_data/           # Visualizations (PCA, heatmaps, pre/post norm)
│   └── edgeR/                        # MDS plots
│
└── {comparison}/                      # One per comparison (e.g., treated_vs_control)
    ├── normalized/                    # Created by prepare_comparison_data
    │   └── (normalized objects copied or recomputed)
    ├── DESeq2/
    │   ├── DESeq2_{comparison}.tsv   # Main DE results
    │   └── report_data/              # Volcano/MA plots
    ├── edgeR/
    │   ├── edgeR_{comparison}.tsv    # Main DE results
    │   └── detail_results/
    └── overlap/                       # DESeq2 vs edgeR overlap plots

final_report.html                      # Master report aggregating all analyses
logs/DE/                               # Execution logs
```
