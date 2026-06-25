####################################################################################################
# Script to normalize samples and generate PCA/heatmap/MDS plots
####################################################################################################

library(data.table)
library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(gplots)
library(pheatmap)
library(ggpubr)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)
library(cowplot)

args <- commandArgs(trailingOnly = T)
count_data_file <- args[1]
txi_file <- args[2]  # Can be empty for featureCount data
experiment_design_file <- args[3]
analysis_type <- args[4]
condition_to_compare_vec_file <- args[5]
output_dir <- args[6]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/DESeq2_func.R"))
source(paste0(script.dir, "/../DE_computation/edgeR_func.R"))

# Load count data
count_dt_original <- readRDS(count_data_file)
if(txi_file != "" && file.exists(txi_file)){
  txi <- readRDS(txi_file)
} else {
  txi <- NULL
}

# Load experiment design from input file (created by computation loading data rule)
experiment_design <- readRDS(experiment_design_file)
condition_to_compare_vec <- readRDS(condition_to_compare_vec_file)

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# DESeq2 normalization and visualization
# ============================================================================
cat("Running DESeq2 normalization and visualization...\n")

res_list <- DESeq2_computation(txi, count_dt_original, experiment_design, condition_to_compare_vec)
dds <- res_list[[1]]
count_dt_normalized <- res_list[[2]]

# Save normalized DESeq2 object
saveRDS(dds, file.path(output_dir, "dds_normalized.RDS"))
saveRDS(count_dt_normalized, file.path(output_dir, "count_data_normalized.RDS"))

# Generate all normalization-specific plots (PCA, heatmaps, etc.) using the standard function
cat("Creating DESeq2 normalization plots...\n")
create_normalization_specific_DESeq2_results(
  output_dir = output_dir,
  dds = dds,
  count_dt = count_dt_normalized,
  condition_to_compare_vec = condition_to_compare_vec
)

# ============================================================================
# edgeR normalization and visualization
# ============================================================================
cat("Running edgeR normalization and visualization...\n")

res_list <- edgeR_computation(txi, count_dt_original, experiment_design, analysis_type, condition_to_compare_vec)
edgeR_DGEList <- res_list[[1]]
fit_tagwise_dispersion_DGEGLM <- res_list[[2]]

# Save edgeR objects
saveRDS(edgeR_DGEList, file.path(output_dir, "edgeR_DGEList_normalized.RDS"))
saveRDS(fit_tagwise_dispersion_DGEGLM, file.path(output_dir, "edgeR_fit_normalized.RDS"))

# Generate all normalization-specific plots (MDS, heatmaps, etc.) using the standard function
cat("Creating edgeR normalization plots...\n")
create_normalization_specific_edgeR_results(
  output_dir = output_dir,
  d = edgeR_DGEList,
  count_dt = count_dt_normalized,
  reduced_plot_design = FALSE
)

cat("Normalization and visualization complete.\n")
cat("Output saved to:", output_dir, "\n")
