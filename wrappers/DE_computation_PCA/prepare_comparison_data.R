####################################################################################################
# Script to prepare comparison-specific normalized data
# If normalize_per_comparison = FALSE: copy from PCA folder
# If normalize_per_comparison = TRUE: recompute normalization with subsetted data
# Uses existing functions from ../DE_computation/DESeq2_func.R and edgeR_func.R
####################################################################################################

library(data.table)

args <- commandArgs(trailingOnly = T)
pca_dds_file <- args[1]
pca_count_file <- args[2]
pca_edger_file <- args[3]
pca_edger_fit_file <- args[4]
count_data_original_file <- args[5]
txi_file <- args[6]
experiment_design_file <- args[7]
normalize_per_comparison <- as.logical(toupper(args[8]))
comparison_conditions <- args[9]  # e.g., "treated_vs_control"
output_dir <- args[10]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/DESeq2_func.R"))
source(paste0(script.dir, "/../DE_computation/edgeR_func.R"))

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Preparing comparison-specific data for:", comparison_conditions, "\n")
cat("normalize_per_comparison:", normalize_per_comparison, "\n")

if (!normalize_per_comparison) {
  # Mode: FALSE - Copy normalized objects from PCA (all samples normalized together)
  cat("Copying normalized data from PCA folder (all samples together)...\n")

  file.copy(pca_dds_file, file.path(output_dir, "dds_normalized.RDS"), overwrite = TRUE)
  file.copy(pca_count_file, file.path(output_dir, "count_data_normalized.RDS"), overwrite = TRUE)
  file.copy(pca_edger_file, file.path(output_dir, "edgeR_DGEList_normalized.RDS"), overwrite = TRUE)
  file.copy(pca_edger_fit_file, file.path(output_dir, "edgeR_fit_normalized.RDS"), overwrite = TRUE)

  cat("Copied normalized data from PCA.\n")

} else {
  # Mode: TRUE - Recompute normalization with only samples from this comparison
  cat("Recomputing normalization with comparison-specific samples...\n")

  # Load original count data
  count_dt_original <- readRDS(count_data_original_file)
  if (txi_file != "" && file.exists(txi_file)) {
    txi <- readRDS(txi_file)
  } else {
    txi <- NULL
  }

  # Parse conditions from comparison name (e.g., "treated_vs_control" -> c("treated", "control"))
  condsToCompare <- strsplit(comparison_conditions, "_vs_")[[1]]

  # Subset data to only include samples from these conditions
  count_dt_comp <- count_dt_original[condition %in% condsToCompare]

  # Create comparison-specific experiment design
  experiment_design_comp <- unique(count_dt_comp[,.(sample_name, condition, patient)])
  setorder(experiment_design_comp, condition, patient)
  experiment_design_comp[,condition := factor(condition, levels = unique(condition))]
  experiment_design_comp[,patient := factor(patient, levels = unique(patient))]

  condition_to_compare_vec <- unique(experiment_design_comp$condition)

  cat("  Conditions:", paste(condition_to_compare_vec, collapse=", "), "\n")
  cat("  Samples in comparison:", nrow(experiment_design_comp), "\n")

  # Run DESeq2 normalization on subsetted data
  cat("  Running DESeq2 normalization...\n")
  res_list <- DESeq2_computation(txi, count_dt_comp, experiment_design_comp, condition_to_compare_vec)
  dds <- res_list[[1]]
  count_dt_comp <- res_list[[2]]

  # Run edgeR normalization on subsetted data
  cat("  Running edgeR normalization...\n")
  # Get analysis_type from the count data structure
  analysis_type <- if (!is.null(txi)) "RSEM" else "feature_count"
  res_list <- edgeR_computation(txi, count_dt_comp, experiment_design_comp, analysis_type, condition_to_compare_vec)
  edgeR_DGEList <- res_list[[1]]
  fit_tagwise_dispersion_DGEGLM <- res_list[[2]]

  # Save comparison-specific normalized objects
  saveRDS(dds, file.path(output_dir, "dds_normalized.RDS"))
  saveRDS(count_dt_comp, file.path(output_dir, "count_data_normalized.RDS"))
  saveRDS(edgeR_DGEList, file.path(output_dir, "edgeR_DGEList_normalized.RDS"))
  saveRDS(fit_tagwise_dispersion_DGEGLM, file.path(output_dir, "edgeR_fit_normalized.RDS"))

  cat("Comparison-specific normalization complete.\n")
}

cat("Output saved to:", output_dir, "\n")
