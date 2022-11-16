library(tximport)
library(data.table)
library(rtracklayer)

run_all <- function(args){
  output_file <- args[1]
  gtf_file <- args[2]
  gtf <- as.data.table(rtracklayer::import(gtf_file, feature.type = "transcript"))
  tx2gene <- unique(gtf[,.(transcript_id , gene_id)])
  file_list <- tail(args,-2)
  sample_names <- gsub(".*kallisto/(.*).kallisto.tsv","\\1",file_list)
  names(file_list) <- sample_names
  txi <- tximport(file_list, type = "kallisto", tx2gene = tx2gene) # Read gene counts
  # head(txi$counts) # See header

  save(txi,file = output_file)
}

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
