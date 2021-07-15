library(data.table)

fread_vector_of_files <- function(file_list,regex = NULL,add_column = "sample"){
  rbindlist(lapply(file_list,function(x){
    res <- fread(x)
    if(is.null(regex)){
      res[,(add_column) := x]
    } else {
      res[,(add_column) := gsub(regex,"\\1",x)]
    }
    
  }))
}

run_all <- function(file_list,output_file){
  res_tab <- fread_vector_of_files(file_list,".*\\/(.*)\\.featureCounts.tsv$")
  res_tab[,sample := make.names(sample)]
  setnames(res_tab,tail(names(res_tab),2)[1],"count")
  res_tab <- dcast.data.table(res_tab,Geneid + Chr + Start + End + Strand + Length ~ sample,value.var = "count")
  write.table(res_tab,file = output_file,quote = F,row.names = F,col.names = T,sep = "\t")
}

# develop and test 2
# file_list <- list.files("/mnt/ssd/ssd_1/snakemake/Mrazlab/sequencing_results/primary_data/190306_MM_mRNA_idelalisib_run2/feature_counts/",pattern = ".featureCounts.tsv$",full.names = T)
# output_file <- "/mnt/ssd/ssd_1/snakemake/Mrazlab/sequencing_results/projects"

# run as Rscript
args <- commandArgs(trailingOnly = T)
file_list <- tail(args,-1)
output_file <- args[1]
run_all(file_list,output_file)