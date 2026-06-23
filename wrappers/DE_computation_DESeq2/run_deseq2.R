####################################################################################################
# Script to run DESeq2 differential expression analysis for a specific comparison
# Uses existing functions from DESeq2_func.R
####################################################################################################

library(data.table)

args <- commandArgs(trailingOnly = T)
dds_file <- args[1]                          # Pre-computed dds.RDS from normalization step
count_data_normalized_file <- args[2]        # count_data_normalized.RDS
experiment_design_file <- args[3]
comparison_name <- args[4]                   # e.g., "treated_vs_control"
pvalue_for_viz <- as.numeric(args[5])
fold_change_threshold <- as.numeric(args[6])
named_in_viz <- as.integer(args[7])
output_dir <- args[8]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/DESeq2_func.R"))

# Load data
dds <- readRDS(dds_file)
count_dt <- readRDS(count_data_normalized_file)

# Load experiment design from input file
experiment_design <- fread(experiment_design_file)

# Parse comparison
condsToCompare <- strsplit(comparison_name, "_vs_")[[1]]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Running DESeq2 result extraction for comparison:", comparison_name, "\n")

# Use existing function to get comparison-specific results
comp_res <- get_comparison_specific_DESeq2_table(
  dds, count_dt, experiment_design,
  condsToCompare,
  output_dir = output_dir,
  p_value_threshold = pvalue_for_viz,
  lfc_threshold = log2(fold_change_threshold)
)

# Rename DESeq2.tsv to DESeq2_<comparison>.tsv
deseq2_tsv <- file.path(output_dir, "DESeq2.tsv")
if (file.exists(deseq2_tsv)) {
  new_name <- paste0("DESeq2_", comparison_name, ".tsv")
  file.rename(deseq2_tsv, file.path(output_dir, new_name))
}

# Generate comparison-specific plots using existing function
create_comparison_specific_DESeq2_results(
  comp_res, dds, count_dt, condsToCompare,
  output_dir = output_dir,
  paired_samples = length(unique(experiment_design$patient)) != nrow(experiment_design),
  TOP = named_in_viz,
  p_value_threshold = pvalue_for_viz,
  lfc_threshold = log2(fold_change_threshold),
  condition_to_compare_vec = unique(experiment_design$condition)
)

cat("DESeq2 analysis complete for:", comparison_name, "\n")
