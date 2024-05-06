`%!in%` <- Negate(`%in%`)
`%!like%` <- Negate(`%like%`)

read_and_prepare_design_data <- function(comparison_vec,experiment_design_file,paired_samples,use_custom_batch_effect_grouping){
  experiment_design <- fread(experiment_design_file)
  if(paired_samples){
    if(use_custom_batch_effect_grouping){
      experiment_design[,patient := batch_group]
    } else {
      experiment_design[,patient := replicate]
    }
  } else {
    experiment_design[,patient := paste0("pat",seq_along(sample_name))]
  }
  
  # Make syntactically valid names using make.names
  for (j in names(experiment_design)[sapply(experiment_design,class) == "character"]) set(experiment_design, j = j, value = make.names(experiment_design[[j]]))
  comparison_vec <- make.names(comparison_vec)
  
  condition_to_compare_vec <- rev(unique(unlist(strsplit(comparison_vec,split= "_vs_"))))
  experiment_design[,condition_order := match(experiment_design$condition,condition_to_compare_vec)]
  setorder(experiment_design,condition_order,patient,na.last = T)
  experiment_design[,sample_name:=make.names(sample_name)]
  experiment_design[,condition_order := NULL]
  experiment_design[,condition := factor(condition,levels = unique(condition))]
  experiment_design[,patient := factor(patient,levels = unique(patient))]
  
  return(list(experiment_design,comparison_vec,condition_to_compare_vec))
}

filterGTF <- function(TSV, geneList = "all", keepGene = TRUE, chrmList = "all", keepChrm = TRUE){
    dt<-TSV
    dt_names<-copy(names(dt)) # c("Geneid","Chr","Feature_name","biotype")

    if(geneList == "all" & !keepGene){
      print("You are removing all genes from the table ! All genes will be kept !")
      filtered <- copy(dt)
    }else if(chrmList == "all" & !keepChrm){
      print("You are removing all chromosomes from the table ! All genes will be kept !")
      filtered <- copy(dt)
    }else{
      if(geneList != "all"){
        print("filtering genes")
        dg <- data.table(Geneid=unlist(strsplit(geneList,split=",",fixed = T)), geneList=TRUE)
        dt <- merge(dt, dg, by="Geneid", all=TRUE)
        dt[is.na(geneList) , geneList := FALSE] # change geneList value from NA to FALSE
      }else{
        dt[, geneList:=TRUE]
      }

      if(chrmList != "all"){
        print("filtering chromosomes")
        dc <- data.table(Chr=unlist(strsplit(chrmList,split=",",fixed = T)), chrmList=TRUE)
        dt <- merge(dt, dc, by="Chr", all=TRUE)
        dt[is.na(chrmList) , chrmList := FALSE] # change geneList value from NA to FALSE
      }else{
        dt[, chrmList:=TRUE]
      }

      setcolorder(dt, dt_names)

      filtered <- dt[geneList == keepGene & chrmList == keepChrm, ][, dt_names, with=FALSE]

      if(keepGene){textGene <- "Gene to keep: "}else{textGene <- "Gene to remove: "}
      if(keepChrm){textChrm <- "Chromosome to keep: "}else{textChrm <- "Chromosome to remove: "}

      print(paste0(textGene,geneList))
      print(paste0(textChrm,chrmList))

      print(paste0("Total genes before filtering: ",length(dt$Geneid)))
      print(paste0("Total genes after filtering: ",length(filtered$Geneid)))
    }

    return(filtered)
  }


read_and_prepare_count_data <- function(counts_file,experiment_design,gtf_filename,analysis_type,geneList,keepGene,chrmList,keepChrm,remove_genes_with_sum_read_count_threshold,remove_genes_with_mean_read_count_threshold){

  if(analysis_type %!like% "mirbase"){
    feat_type <- "gene"
    annotate_by<- c("gene_id","seqnames","gene_name", "gene_biotype")
    gtf_gene_tab <- as.data.table(rtracklayer::import(gtf_filename, feature.type = feat_type))[,annotate_by, with=F]
    setnames(gtf_gene_tab,c("Geneid","Chr","Gene_name","biotype"))
    gtf_gene_tab <- gtf_gene_tab[!is.na(biotype)]
    gtf_gene_tab[is.na(Gene_name) | Gene_name == "",Gene_name := Geneid]
    gtf_gene_tab[,duplicated := .N,by = "Gene_name"]
    gtf_gene_tab[duplicated > 1,Gene_name := paste0(Gene_name,"__",Geneid)]
    gtf_gene_tab[,duplicated := NULL]
    setnames(gtf_gene_tab,"Gene_name","Feature_name")
    setkey(gtf_gene_tab,"Geneid")

    gtf_gene_tab<-filterGTF(gtf_gene_tab,geneList,keepGene,chrmList,keepChrm)
  }

  if(analysis_type %like% "featureCount" | analysis_type %like% "mirbase"){
    txi <- NULL
    
    count_dt <- fread(counts_file)
    setnames(count_dt,make.names(colnames(count_dt)))
    setnames(count_dt,"mirna","Geneid", skip_absent = T)
    count_dt[,c("Chr","Start","End","Strand","Length") := NULL]
    count_dt <- melt(count_dt,measure.vars = experiment_design$sample_name,variable.name = "sample_name",value.name = "count")
    count_dt[,sample_name := factor(sample_name,levels = experiment_design$sample_name)]
    count_dt[,sum_count := sum(count),by = Geneid]
    count_dt[,mean_count := mean(count),by = Geneid]
    print(paste0("Total number of genes in the data: ",length(count_dt$Geneid)))
    filterM <- length(count_dt[mean_count > remove_genes_with_mean_read_count_threshold,]$Geneid)
    count_dt <- count_dt[mean_count > remove_genes_with_mean_read_count_threshold,]
    print(paste0("Total number of genes after RowMeans filtering: ",filterM," - mean cut-off ",remove_genes_with_mean_read_count_threshold))
    filterS <- length(count_dt[sum_count > remove_genes_with_sum_read_count_threshold,]$Geneid)
    count_dt <- count_dt[sum_count > remove_genes_with_sum_read_count_threshold,]
    print(paste0("Total number of genes after RowSums filtering: ",filterS," - sum cut-off ",remove_genes_with_sum_read_count_threshold))
    
  } else {
    load(counts_file)
    #load("DE_RSEM/complete_RSEM_table.RData")
    colnames(txi$counts) <- make.names(colnames(txi$counts))
    colnames(txi$length) <- make.names(colnames(txi$length))
    colnames(txi$abundance) <- make.names(colnames(txi$abundance))
    
    # Remove "gene:" or "transcript": added by RSEM - doesn't go well when merged with Ensembl or other annotation-based information
    rownames(txi$abundance) <- gsub("^gene:", "", rownames(txi$abundance))
    rownames(txi$counts) <- gsub("^gene:", "", rownames(txi$counts))
    rownames(txi$length) <- gsub("^gene:", "", rownames(txi$length))
    rownames(txi$abundance) <- gsub("^transcript:","", rownames(txi$abundance))
    rownames(txi$counts) <- gsub("^transcript:", "", rownames(txi$counts))
    rownames(txi$length) <- gsub("^transcript:", "", rownames(txi$length))
    
    # Keep only non-zero expressed genes that are in ensembl DB
    #   If gene/isoform has only counts 0.x it will be kept in the output but will have 0 counts in final DE tables
    keep<-rowSums(txi$counts)>0 & rownames(txi$counts) %in% gtf_gene_tab$Geneid
    print(paste0("Total number of genes in the data: ",length(keep)))

    txi$abundance<-txi$abundance[keep,]
    txi$counts<-txi$counts[keep,]
    txi$length<-txi$length[keep,]

    # filter for low counts per mean of counts
    filterM<-rowMeans(txi$counts)>remove_genes_with_mean_read_count_threshold & rownames(txi$counts) %in% gtf_gene_tab$Geneid
    print(paste0("Total number of genes after RowMeans filtering: ",length(filterM)," - mean cut-off ",remove_genes_with_mean_read_count_threshold))

    txi$abundance<-txi$abundance[filterM,]
    txi$counts<-txi$counts[filterM,]
    txi$length<-txi$length[filterM,]

    # filter for low counts per sum of counts
    filterS<-rowSums(txi$counts)>remove_genes_with_sum_read_count_threshold & rownames(txi$counts) %in% gtf_gene_tab$Geneid
    print(paste0("Total number of genes after RowSums filtering: ",length(filterS)," - sum cut-off ",remove_genes_with_sum_read_count_threshold))

    txi$abundance<-txi$abundance[filterS,]
    txi$counts<-txi$counts[filterS,]
    txi$length<-txi$length[filterS,]

    #set correct col order
    txi$counts<-txi$counts[,match(experiment_design$sample_name, colnames(txi$counts))]
    txi$abundance<-txi$abundance[,match(experiment_design$sample_name, colnames(txi$abundance))]
    txi$length<-txi$length[,match(experiment_design$sample_name, colnames(txi$length))]
    
    count_dt<-as.data.table(txi$counts,keep.rownames = T)
    setnames(count_dt,"rn","Geneid")
    count_dt <- melt(count_dt,measure.vars = experiment_design$sample_name,variable.name = "sample_name",value.name = "count")
    count_dt[,sample_name := factor(sample_name,levels = experiment_design$sample_name)]
    count_dt[,sum_count := sum(count),by = Geneid]
    count_dt[,mean_count := mean(count),by = Geneid]
    
    #rename to Gene names from Ensembl
    # rownames(txi$abundance) <- gtf_gene_tab[rownames(txi$abundance)]$Feature_name
    # rownames(txi$counts) <- gtf_gene_tab[rownames(txi$counts)]$Feature_name
    # rownames(txi$length) <- gtf_gene_tab[rownames(txi$length)]$Feature_name
    #
    txi$length[txi$length == 0] <- 1 # If gene has length 0 replace it with 1 to avoid error later, might be source of bias and/or error; https://support.bioconductor.org/p/84304/
    
    # cts <- txi$counts
    # anyNA(cts) # Should be false
    # normMat <- txi$length
    # normMat <- as.data.frame(normMat/exp(rowMeans(log(normMat))))
  }

  if(analysis_type %like% "mirbase"){
    count_dt[,`:=`(Feature_name=Geneid, biotype="miRNA")]
  }else{
    count_dt <- merge(gtf_gene_tab,count_dt,by = "Geneid")
  }
  setnames(count_dt,"Geneid","Ensembl_Id")
  count_dt <- merge(experiment_design[,.(sample_name, condition, patient)],count_dt,by = "sample_name")
  setcolorder(count_dt,c("Ensembl_Id","Feature_name","biotype","sample_name","condition","patient","count","sum_count","mean_count"))
  setkey(count_dt,Ensembl_Id,condition,patient,sample_name)
  
  return(list(count_dt,txi))
}
