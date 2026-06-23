####################################################################################################
# Script to load raw count data and create count_data_original and txi objects
####################################################################################################

library(data.table)

args <- commandArgs(trailingOnly = T)
counts_file <- args[1]
experiment_design_file <- args[2]
gtf_filename <- args[3]
analysis_type <- args[4]
geneList <- args[5]
keepGene <- as.logical(toupper(args[6]))
chrmList <- args[7]
keepChrm <- as.logical(toupper(args[8]))
remove_genes_with_sum_read_count_threshold <- as.integer(args[9])
remove_genes_with_mean_read_count_threshold <- as.integer(args[10])
output_dir <- args[11]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/load_data_func.R"))

# Read experiment design
experiment_design <- fread(experiment_design_file)
if("patient" %in% colnames(experiment_design)){
  # Patient column already exists (from paired design)
} else {
  experiment_design[,patient := paste0("pat",seq_along(sample_name))]
}

# Make syntactically valid names
for (j in names(experiment_design)[sapply(experiment_design,class) == "character"]) {
  set(experiment_design, j = j, value = make.names(experiment_design[[j]]))
}
experiment_design[,sample_name:=make.names(sample_name)]
experiment_design[,condition := factor(condition,levels = unique(condition))]
experiment_design[,patient := factor(patient,levels = unique(patient))]

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
