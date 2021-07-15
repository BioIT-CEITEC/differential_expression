library(data.table)
library(rjson)

run_all <- function(config_file,output_file){
  
  tab <- as.data.table(fromJSON(file = config_file))
  tab <- tab[,list(sample = full_name,name = full_name,condition,patient = donor)]
  
  if(any(condition) == "control"){
    
  }
  tab
  
  
}

# develop and test 2
config_file <- "/mnt/ssd/ssd_1/snakemake/CFBioinformatics/sequencing_results/projects/DE_test/mRNA_DE_featureCount/DE_test.DE.config.json"
output_file <- "/mnt/ssd/ssd_1/snakemake/Mrazlab/sequencing_results/projects"

# run as Rscript
# script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
# args <- commandArgs(trailingOnly = T)
# file_list <- tail(args,-1)
# output_file <- args[1]
# run_all(file_list)