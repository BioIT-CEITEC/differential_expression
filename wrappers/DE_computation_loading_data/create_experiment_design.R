####################################################################################################
# Script to create experiment design file from sample table
# Uses the existing read_and_prepare_design_data function from load_data_func.R
####################################################################################################

library(data.table)

args <- commandArgs(trailingOnly = T)
sample_tab_file <- args[1]
output_file <- args[2]
paired_samples <- as.logical(toupper(args[3]))
use_custom_batch_effect_grouping <- as.logical(toupper(args[4]))

# Read sample table (created by BioRoots load_sample())
experiment_design <- fread(sample_tab_file)

# Apply the same logic as read_and_prepare_design_data for patient assignment
if(paired_samples){
  if(use_custom_batch_effect_grouping && "batch_group" %in% colnames(experiment_design)){
    experiment_design[,patient := batch_group]
  } else {
    experiment_design[,patient := replicate]
  }
} else {
  experiment_design[,patient := paste0("pat",seq_along(sample_name))]
}

# Make syntactically valid names using make.names
for (j in names(experiment_design)[sapply(experiment_design,class) == "character"]) {
  set(experiment_design, j = j, value = make.names(experiment_design[[j]]))
}

experiment_design[,sample_name:=make.names(sample_name)]
experiment_design[,condition := factor(condition,levels = unique(condition))]
experiment_design[,patient := factor(patient,levels = unique(patient))]

# Save the experiment design
fwrite(experiment_design, output_file, sep = "\t")

cat("Experiment design created:", output_file, "\n")
