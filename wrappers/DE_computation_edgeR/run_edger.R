####################################################################################################
# Script to run edgeR differential expression analysis for a specific comparison
# Uses existing functions from edgeR_func.R
####################################################################################################

library(data.table)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)
library(cowplot)
library(gplots)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(edgeR)
library(limma)
library(svglite)
library(plotly)
library(htmlwidgets)
library(rtracklayer)
library(ggupset)
library(UpSetR)

args <- commandArgs(trailingOnly = T)
edgeR_DGEList_file <- args[1]              # Pre-computed edgeR_DGEList.RDS from normalization step
fit_tagwise_dispersion_file <- args[2]     # edgeR fit object from normalization step
count_data_normalized_file <- args[3]      # count_data_normalized.RDS
comparison_name <- args[4]                 # e.g., "treated_vs_control"
pvalue_for_viz <- as.numeric(args[5])
fold_change_threshold <- as.numeric(args[6])
named_in_viz <- as.integer(args[7])
output_dir <- args[8]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/DESeq2_func.R"))
source(paste0(script.dir, "/../DE_computation/edgeR_func.R"))

# Load data
edgeR_DGEList <- readRDS(edgeR_DGEList_file)
fit_tagwise_dispersion_DGEGLM <- readRDS(fit_tagwise_dispersion_file)
count_dt <- readRDS(count_data_normalized_file)

# Parse comparison
condsToCompare <- strsplit(comparison_name, "_vs_")[[1]]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Running edgeR result extraction for comparison:", comparison_name, "\n")

# Create normalization-specific plots and data for this comparison
cat("Creating normalization-specific edgeR results...\n")
create_normalization_specific_edgeR_results(
  output_dir = output_dir,
  d = edgeR_DGEList,
  count_dt = count_dt[condition %in% condsToCompare],
  reduced_plot_design = TRUE
)

# Use existing function to get comparison-specific results
res_list <- get_comparison_specific_edgeR_table(
  fit_tagwise_dispersion_DGEGLM, edgeR_DGEList, count_dt,
  condsToCompare,
  output_dir = output_dir,
  p_value_threshold = pvalue_for_viz,
  lfc_threshold = log2(fold_change_threshold)
)
edgeR_comp_res <- res_list[[1]]
lrt_tgw <- res_list[[2]]

# Rename edgeR.tsv to edgeR_<comparison>.tsv
edger_tsv <- file.path(output_dir, "edgeR.tsv")
if (file.exists(edger_tsv)) {
  new_name <- paste0("edgeR_", comparison_name, ".tsv")
  file.rename(edger_tsv, file.path(output_dir, new_name))
}

# Generate comparison-specific plots using existing function
create_comparison_specific_edgeR_results(
  edgeR_comp_res, lrt_tgw, condsToCompare,
  output_dir = output_dir,
  paired_samples = length(unique(count_dt[,.(sample_name, patient)])) != nrow(unique(count_dt[,.(sample_name, patient)])),
  TOP = named_in_viz,
  p_value_threshold = pvalue_for_viz,
  lfc_threshold = log2(fold_change_threshold)
)

cat("edgeR analysis complete for:", comparison_name, "\n")
