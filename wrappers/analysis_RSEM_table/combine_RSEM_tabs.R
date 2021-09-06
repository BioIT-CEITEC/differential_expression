library(tximport)
library(data.table)
library(rtracklayer)

run_all <- function(args){
  output_file <- args[1]
  file_list <- tail(args,-1)
  sample_names <- gsub(".*RSEM/(.*).genes.results","\\1",file_list)
  names(file_list) <- sample_names
  txi <- tximport(file_list, type = "rsem") # Read gene counts
  # head(txi$counts) # See header

  save(txi,file = output_file)
}

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
