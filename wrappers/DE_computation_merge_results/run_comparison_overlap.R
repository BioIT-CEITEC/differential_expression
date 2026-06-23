####################################################################################################
# Script to compare DESeq2 and edgeR results and generate overlap visualization
# Uses existing functions from DESeq2_edgeR_edgeR_overlap_func.R
####################################################################################################

args <- commandArgs(trailingOnly = T)
deseq2_results_file <- args[1]             # DESeq2.tsv from DESeq2 rule
edger_results_file <- args[2]              # edgeR.tsv from edgeR rule
comparison_name <- args[3]                 # e.g., "treated_vs_control"
pvalue_for_viz <- as.numeric(args[4])
fold_change_threshold <- as.numeric(args[5])
output_dir <- args[6]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/DESeq2_edgeR_overlap_func.R"))

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Running DESeq2-edgeR overlap analysis for comparison:", comparison_name, "\n")

# Load results
deseq2_res <- fread(deseq2_results_file)
edger_res <- fread(edger_results_file)

# Use existing function to generate overlap analysis
comparison_specific_edgeR_DESeq2_overlap(
  output_dir,
  deseq2_res,
  edger_res,
  p_value_threshold = pvalue_for_viz,
  lfc_threshold = log2(fold_change_threshold)
)

cat("DESeq2-edgeR overlap analysis complete for:", comparison_name, "\n")
