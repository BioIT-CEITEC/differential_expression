library(tximport)
library(data.table)
library(readr)

run_all <- function(args){
  output_file <- args[1]
  tx2gene_file <- args[2]
  tx2gene <- read_delim(file.path(tx2gene_file),delim='\t')
  file_list <- tail(args,-2)
  sample_names <- gsub(".*salmon_.*/(.*).salmon.sf","\\1",file_list)
  names(file_list) <- sample_names
  txi <- tximport(file_list, type = "salmon", tx2gene = tx2gene) # Read gene counts
  # head(txi$counts) # See header

  save(txi,file = output_file)
}

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
