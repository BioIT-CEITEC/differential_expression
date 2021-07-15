library(tximport)
library(data.table)
library(rtracklayer)

run_all <- function(args){
  output_file <- args[1]
  assembly_gtf <- args[2]
  if(assembly_gtf == "ref_not_from_trans_assembly"){
    ref_from_trans_assembly <- F
  } else {
    ref_from_trans_assembly <- T
  }
  file_list <- tail(args,-2)
  sample_names <- gsub(".*rsem_counts/(.*).genes.results","\\1",file_list)
  names(file_list) <- sample_names
  txi.rsem <- tximport(file_list, type = "rsem") # Read gene counts
  # head(txi.rsem$counts) # See header
  txi<-txi.rsem
  
  if(ref_from_trans_assembly){
    feat_type <- "transcript"
    annotate_by<- c("gene_name","gene_id", "strand")  
    gtf_file <- "/mnt/ssd/ssd_3/references/termite/Peska_transcriptome/annot/Peska_transcriptome.gtf"
    gtf <- as.data.table(rtracklayer::import(gtf_file))[type == feat_type, c("seqnames","start","end",annotate_by), with=F]
    
    gtf[,first_annot := tstrsplit(gtf$gene_name,split = ",")[3]]
    gtf[,c("first_annot_ID","first_annot_info","first_annot_score") := tstrsplit(first_annot, split="\\|")]
    gtf[,first_annot_score := as.numeric(first_annot_score)]
    # gtf[,second_annot := tstrsplit(gtf$gene_name,split = ",")[4]]
    # gtf[,c("second_annot_ID","second_annot_info","second_annot_score") := tstrsplit(second_annot,split="\\|")]
    # gtf[,second_annot_score := as.numeric(second_annot_score)]
    suppressWarnings(gtf[,c("first_annot","second_annot") := NULL])
    
    gtf[grep("PF",first_annot_info), first_annot_DB := "Pfam"]
    gtf[!is.na(suppressWarnings(as.numeric(first_annot_info))), first_annot_DB := "UniProt"]
    gtf[is.na(first_annot_DB), first_annot_DB := "TransDecoder"]
    gtf[first_annot_DB == "Pfam", final_annot := paste(first_annot_info,first_annot_ID,sep = "|")]
    gtf[first_annot_DB == "UniProt", final_annot := first_annot_ID ]
    # gtf[first_annot_DB == "UniProt", final_annot := lapply(first_annot_ID, function(x) {
    #   name = system(paste0("wget -O - https://www.uniprot.org/uniprot/",x," 2> /dev/null | grep -oP \"title>.*</title\"|sed -e 's/title>//'|cut -f1 -d' '"))
    #   return(ifelse(name != "", paste0(first_annot_ID, name, "|"), first_annot_ID))
    # })]
    # gtf[first_annot_DB == "UniProt", final_annot := {
    #   name = system(paste0("wget -O - https://www.uniprot.org/uniprot/",first_annot_ID," 2> /dev/null | grep -oP \"title>.*</title\"|sed -e 's/title>//'|cut -f1 -d' '"))
    #   list(ifelse(name != "", paste0(first_annot_ID, name, "|"), first_annot_ID))
    #   }]
    # gtf[,final_annot := final_annot[1],by = seqnames]
    gtf[is.na(final_annot), final_annot := gsub("GENE.(.*?)~~(.*)","\\2",gene_id) ] #paste0("UNK_transcript_",seq_along(first_annot_DB))]
    fwrite(gtf[,.(gene_id, transcript_id=seqnames, protein_id=final_annot, protein_source=first_annot_DB)],
           file = gsub("complete.RSEM.RData","annot_to_seq_ids.tsv",output_file),sep = "\t")
    
    txi_mat_merge_by_gtf_annot <- function(mat,type,gtf){
      mat <- as.data.table(mat,keep.rownames = T)
      setnames(mat,"rn","gene_id")
      mat <- merge(gtf, mat, by = "gene_id")
      mat[,gene_id := NULL]
      setnames(mat,"final_annot","gene_id")
      if(type == "counts" | type == "abundance"){
        mat <- mat[,lapply(.SD,sum),by = gene_id]
      } else {
        mat <- mat[,lapply(.SD,max),by = gene_id]
      }
      rn <- mat$gene_id
      mat <- as.matrix(mat[,-1,with = F])
      row.names(mat) <- rn
      return(mat)
    }

    # gtf <- gtf[,.(gene_id,transcript_id)]
    gtf[,final_annot:=final_annot[1],by = seqnames]
    gtf = gtf[,.(final_annot, gene_id)]
    txi$counts <- txi_mat_merge_by_gtf_annot(txi$counts,"counts",gtf)
    txi$abundance <- txi_mat_merge_by_gtf_annot(txi$abundance,"abundance",gtf)
    txi$length <- txi_mat_merge_by_gtf_annot(txi$length,"length",gtf)

    # gtf[,ref_seq_id := gsub("GENE.(.*?)~~.*","\\1",gene_id)]
    # setnames(gtf,"gene_id","gtf_id")
    # setnames(gtf,"final_annot","final_annot_id")
    # setcolorder(gtf,c("final_annot_id","transcript_id","gtf_id"))
    # fwrite(gtf,file = gsub("complete.RSEM.RData","annot_to_seq_ids.tsv",output_file),sep = "\t")
  }
  
  save(txi,file = output_file)
}

# to test
# args <- "/mnt/nfs/shared/CFBioinformatics/all_projects/stage283_Termiti_peska.qseq_DE/mRNA_DE_RSEM/complete.RSEM.RData"
# args <- c(args,"/mnt/ssd/ssd_3/references/termite/Peska_transcriptome/annot/Peska_transcriptome.gtf")
# args <- c(args,grep("techrep",list.files("/mnt/nfs/shared/CFBioinformatics/all_projects/stage283_Termiti_peska.qseq_DE/input_files/rsem_counts",pattern = ".genes.results$",full.names = T),invert = T,value = T))


# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
