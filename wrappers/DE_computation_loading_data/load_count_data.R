####################################################################################################
# Script to load raw count data and create count_data_original and txi objects
####################################################################################################

library(data.table)

args <- commandArgs(trailingOnly = T)
experiment_design_file <- args[1]
counts_file <- args[2]
comparison_vec <- strsplit(args[3],split = ",")[[1]]
gtf_filename <- args[4]
analysis_type <- args[5]
paired_samples <- as.logical(toupper(args[6]))
use_custom_batch_effect_grouping <- as.logical(toupper(args[7]))
geneList <- args[8]
keepGene <- as.logical(toupper(args[9]))
chrmList <- args[10]
keepChrm <- as.logical(toupper(args[11]))
remove_genes_with_sum_read_count_threshold <- as.integer(args[12])
remove_genes_with_mean_read_count_threshold <- as.integer(args[13])
output_dir <- args[14]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/load_data_func.R"))

res_list <- read_and_prepare_design_data(comparison_vec,experiment_design_file,paired_samples,use_custom_batch_effect_grouping)
  experiment_design <- res_list[[1]]
  comparison_vec <- res_list[[2]]
  condition_to_compare_vec <- res_list[[3]]

# Read and prepare count data
res_list <- read_and_prepare_count_data(
  counts_file, experiment_design, gtf_filename, analysis_type,
  geneList, keepGene, chrmList, keepChrm,
  remove_genes_with_sum_read_count_threshold, remove_genes_with_mean_read_count_threshold
)
count_dt_original <- res_list[[1]]
txi <- res_list[[2]]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save objects
saveRDS(comparison_vec, file.path(output_dir, "comparison_vec.RDS"))
saveRDS(condition_to_compare_vec, file.path(output_dir, "condition_to_compare_vec.RDS"))
saveRDS(count_dt_original, file.path(output_dir, "count_data_original.RDS"))

# Always save txi - NULL for featureCount-based analyses, actual txi object for RSEM/kallisto/salmon
if(is.null(txi)){
  cat("Analysis type does not produce txi object (featureCount-based). Saving empty txi.\n")
  saveRDS(NULL, file.path(output_dir, "txi.RDS"))
} else {
  saveRDS(txi, file.path(output_dir, "txi.RDS"))
}

# Also save as RData for compatibility
save(count_dt_original, txi, file = file.path(output_dir, "count_data_objects.RData"))

cat("Successfully created count_data_original and txi objects\n")
cat("Output saved to:", output_dir, "\n")
