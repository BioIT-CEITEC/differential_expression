# Differential Expression Analysis Pipeline

A Snakemake-based workflow for differential gene expression analysis supporting multiple quantification methods (featureCounts, RSEM, Salmon, Kallisto) with DESeq2 and edgeR statistical analysis.

## Overview

This pipeline performs differential expression (DE) analysis on RNA-seq count data:

1. **Aggregates** per-sample quantification files into a unified count table
2. **Normalizes** samples using both DESeq2 and edgeR methods
3. **Computes** differential expression for each comparison of conditions
4. **Generates** comprehensive visualizations (PCA, heatmaps, volcano plots, MA plots)
5. **Compares** DESeq2 and edgeR results to identify consensus DE genes
6. **Produces** an HTML report summarizing all results

## Supported Quantification Methods

- **featureCounts**: Gene/exon/transcript/UTR-level counts
- **RSEM**: Transcript and gene abundance estimates
- **Salmon**: Alignment-based and mapping-based quantification
- **Kallisto**: Pseudo-alignment quantification
- **miRBase**: miRNA counting (when `is_mirna = true`)

## Requirements

- **Snakemake** >= 5.18.0
- **Conda/Mamba** for environment management
- **config.json** file with sample information and parameters

## Quick Start

### 1. Prepare Configuration

Copy the example configuration and edit it for your experiment:

```bash
cp config.example.json config.json
```

Edit `config.json` to specify:
- Organism and reference genome/GTF
- Sample metadata (sample names, conditions, replicates)
- Which quantification methods were used
- Conditions to compare
- Filtering thresholds

### 2. Run the Pipeline

```bash
# Dry run to verify workflow
snakemake --cores 16 --dryrun

# Full execution
snakemake --cores 16

# With resource limits
snakemake --cores 16 --resources mem_mb=65000
```

## Configuration Options

### Key Parameters in `config.json`

| Parameter | Description | Default |
|-----------|-------------|---------|
| `organism` | Organism name | Required |
| `organism_gtf` | Path to GTF annotation file | Required |
| `samples` | Sample metadata table | Required |
| `conditions_to_compare` | Conditions to compare ("all" or "cond1:cond2") | "all" |
| `normalize_data_per_comparison` | Normalize all samples together (FALSE) or per comparison (TRUE) | false |
| `paired_replicates` | Samples are paired by replicate ID | false |
| `pvalue_for_viz` | P-value threshold for visualizations | 0.05 |
| `fold_change_threshold` | Log2 fold-change threshold | 2 |
| `remove_genes_with_sum_read_count_threshold` | Minimum total read count to keep gene | 0 |
| `biotypes` | Gene biotypes to include | "all" |

## Output Structure

```
DE_{method}/
├── loading_data/              # Raw count data objects
├── PCA/                       # Normalization visualizations (all samples)
│   ├── DESeq2/report_data/    # PCA, heatmaps, normalization plots
│   └── edgeR/                 # MDS plots
└── {comparison}/              # Per-comparison results
    ├── DESeq2/
    │   ├── DESeq2_{comparison}.tsv      # Main DE results
    │   └── report_data/                 # Volcano/MA plots
    ├── edgeR/
    │   ├── edgeR_{comparison}.tsv       # Main DE results
    │   └── detail_results/
    └── overlap/               # DESeq2 vs edgeR overlap

final_report.html              # Comprehensive HTML report
logs/DE/                       # Execution logs
```

## Output Files

### Main Results Tables

- **DESeq2_{comparison}.tsv**: DESeq2 differential expression results
  - Columns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, gene_name, gene_biotype, normalized counts, raw counts
  
- **edgeR_{comparison}.tsv**: edgeR differential expression results
  - Similar columns as DESeq2 output

### Visualizations

- **PCA plots**: Sample clustering based on variance-stabilized counts
- **Heatmaps**: Sample-to-sample similarity based on gene expression
- **Volcano plots**: Significance vs fold-change for DE genes
- **MA plots**: Mean expression vs fold-change
- **MDS plots**: Multidimensional scaling (edgeR)
- **Overlap plots**: Venn/UpSet diagrams showing DE gene overlap between methods

## Choosing Between DESeq2 and edgeR

| Method | Best For | Characteristics |
|--------|----------|-----------------|
| **DESeq2** | High specificity, many replicates | Conservative, robust to outliers, independent filtering |
| **edgeR** | Low replicates (<12), exploratory analysis | More sensitive, less conservative |

For most analyses, we recommend examining both methods and focusing on genes identified by both tools.

## Normalization Modes

### `normalize_data_per_comparison = false` (Default)
All samples are normalized together once. This is appropriate when:
- You want consistent normalization across all comparisons
- All samples are from the same experiment/batch

### `normalize_data_per_comparison = true`
Each comparison is normalized separately using only the samples in that comparison. This is appropriate when:
- Comparing very different condition groups
- You want to minimize batch effects between unrelated conditions

## Troubleshooting

### Common Issues

1. **"No rules needed"**: Check that your `config.json` specifies at least one quantification method
2. **Missing input files**: Verify that `qc_reports/{sample}/{method}/` directories contain expected files
3. **Memory errors**: Increase memory limits with `--resources mem_mb=XXXXX`
4. **Conda environment failures**: Try `snakemake --use-conda --rerun-incomplete`

### Getting Help

- Check `logs/DE/*.log` files for detailed error messages
- Run with `--verbose` for more detailed output
- Use `--dryrun` to verify the workflow before execution

## Citation

If you use this pipeline in your research, please cite:
- **DESeq2**: Love et al. (2014) [genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
- **edgeR**: Robinson et al. (2010) [bioinformatics.oxfordjournals.org/content/26/1/139](https://doi.org/10.1093/bioinformatics/btp616)

## License

This pipeline is provided as-is for research purposes.
