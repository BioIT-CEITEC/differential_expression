####################################################################################################
####################################################################################################
#
# Script to calculate differential gene expression using DESeq2 and edgeR package
# Designed for mRNA differential gene expression analysis based on Ensembl results using feature_count, htseq-count, STAR counts
#
####################################################################################################
############################## ALWAYS CHECK DESIGNS AND CONTRASTS!!!! ##############################
####################################################################################################
#
# INPUT_COUNTS - table with raw genes counts -> first column = gene id, first row = patient id
#
# ALWAYS CHECK DESIGN FOR BOTH EDGER AND DESEQ2 AND SELELECT CORRECT COEF/CONTRAST BASED ON YOUR BIOL. QUESTION
# CHECK THE INTERCEPT/NON-INTERCEPT AS WELL!
#
# parsed_ensembl - table with gene_id, gene_name and gene_biotype based in Ensembl
# biotypes - table with name and biotype_group where name equals gene biotype
####################################################################################################
# TODO: Add tSNE clustering
# TODO: Add shared DE genes visualization (inspiration at VBCF BioComp)
####################################################################################################

#devtools::install_github("r-lib/later")



run_all <- function(args){
  library(data.table)

  config_file <- args[1]
  counts_file <- args[2]
  OUTPUT_DIR <- args[3]
  sqlite_path <- args[4]
  biotypes_file <- args[5]
  condsToCompare <- rev(strsplit(args[6],"_vs_")[[1]])
  counts_type <- args[7]
  organism <- args[8]
  use_tag_to_pair_samples <- as.logical(toupper(args[9]))
  ref_from_trans_assembly <- as.logical(toupper(args[10]))

  config_tab <- as.data.table(rjson::fromJSON(file = config_file))
  config_tab[,condition := sapply(1:length(config_tab$samples),function(x) config_tab$samples[[x]]$condition)]
  config_tab[,replicate := sapply(1:length(config_tab$samples),function(x) config_tab$samples[[x]]$replicate)]
  config_tab[,full_name := sapply(1:length(config_tab$samples),function(x) config_tab$samples[[x]]$sample_name)]


  #remove unused samples
  condition_design <- config_tab$conditions_to_compare[1]
  if(condition_design != "all"){
    #condition_design <- strsplit(condition_design,",")[[1]]
    #condition_design <- condition_design[which(grepl(condsToCompare[1],condition_design) & grepl(condsToCompare[2],condition_design))]
    #condition_design <- strsplit(condition_design,":")[[1]]
    setkey(config_tab,condition)
    #config_tab <- config_tab[condition_design,]
    config_tab <- config_tab[condsToCompare,]
  } else {
    setorder(config_tab,condition)
  }
  #print(config_tab)

  #set control samples as first
  if(any(config_tab$condition == "control")){
    config_tab <- rbind(config_tab[condition == "control"],config_tab[condition != "control"])
  }

  #create sample description tab
  if(use_tag_to_pair_samples == T){
    samples_desc <- config_tab[,list(sample = full_name, condition, replicate, patient = replicate)]
  } else {
    if(length(unique(config_tab$donor)) > 1){
      samples_desc <- config_tab[,list(sample = full_name, condition, replicate, patient = donor)]
    } else {
      samples_desc <- config_tab[,list(sample = full_name, condition, replicate, patient = info)]
    }
  }
  samples_desc <- unique(samples_desc)

  # Make syntactically valid names using make.names
  for (j in names(samples_desc)[sapply(samples_desc,class) == "character"]) set(samples_desc, j = j, value = make.names(samples_desc[[j]]))
  condsToCompare <- make.names(condsToCompare)

  # select biotypes
  biotypes <- read.table(biotypes_file, header=T, sep="\t", stringsAsFactors = F)
  featuresToAnalyze <- config_tab$biotypes[1]

  if(featuresToAnalyze=="all"){
    print("Analyzing all the biotypes.")
    featuresToAnalyze2<-c(biotypes[, 1])
  }else{
    featuresToAnalyze <- strsplit(featuresToAnalyze,",")[[1]]
    print(paste0("Analyzing following biotypes", paste(biotypes[biotypes$biotype_group %in% featuresToAnalyze, 1],sep = ", "), "."))
    featuresToAnalyze2<-c(biotypes[biotypes$biotype_group %in% featuresToAnalyze, 1])
  }

  #set other parameters
  P_THRESHOLD <- as.numeric(config_tab$pvalue_for_viz[1])
  FOLD_CHANGE <- as.numeric(config_tab$fold_change_threshold[1])
  LFC_THRESHOLD <- log(FOLD_CHANGE, 2)
  TOP <- as.integer(config_tab$named_in_viz[1]) # How many top DE genes should be plotted

  #is intercept
  INTERCEPT <- samples_desc$condition[1] == "control" # c('TRUE', 'FALSE') - do we want to calculate with intercept - common starting point for all comparisons - are the comparisons relative to something? (the intercept); changes shrunking of logFC; https://support.bioconductor.org/p/69341/; turning off intercept for paired samples doesn't make any sense

  #is paired
  PAIRED <- all(sort(unique(samples_desc[condition == condsToCompare[1]]$patient)) == sort(unique(samples_desc[condition == condsToCompare[2]]$patient))) # c('TRUE', 'FALSE') - do we have paired data? This could be data in different time points from the same patient/samples
  PAIRED <- PAIRED && length(unique(samples_desc[condition == condsToCompare[1]]$patient)) > 1

  # Correct intercept in case PAIRED==T but INTERCEPT==F - it doesn't make too much sense
  if(PAIRED==TRUE){
    print("Making sure INTERCEPT==T because PAIRED==T. If you have specific reason not to set INTERCEPT==T if PAIRED==T please comment this section.")
    INTERCEPT<-TRUE
  } else {

    samples_desc[,patient := paste(patient,seq_along(patient),sep = "_")]
  }
  # ####################################################################################################
  # # Function to check installed packages and installing them if they are missing for CRAN and Bioconductor packages
  # # WARNING: This style of bioc installation works for R < 3.5 Bioconductor installation
  # packages<-c("RColorBrewer", "gplots", "rgl", "pheatmap", "ggrepel", "ggplot2", "dplyr", "ggpubr") # CRAN packages
  # packages_bioc<-c("ensembldb", "DESeq2", "edgeR", "vsn", "limma") # Bioconductor packages
  #
  # if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  #   print(paste("Installing CRAN packages ", setdiff(packages, rownames(installed.packages())), " because they are missing!", sep=""))
  #   install.packages(setdiff(packages, rownames(installed.packages())))
  # }
  #
  # if (length(setdiff(packages_bioc, rownames(installed.packages()))) > 0) {
  #   print(paste("Installing Bioconductor packages ", setdiff(packages_bioc, rownames(installed.packages())), " because they are missing!", sep=""))
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite(setdiff(packages_bioc, rownames(installed.packages())), suppressUpdates=T)
  # }
  #
  # ####################################################################################################

  library("RColorBrewer")

  hmcol<-colorRampPalette(brewer.pal(9, "GnBu"))(100)

  pdf<-NULL
  if (!(is.null(pdf))){
    pdf(pdf, paper="a4")
  }

  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  #start analysis create output dir and read

  dir.create(OUTPUT_DIR, recursive = T,showWarnings = F)
  conds<-factor(c(samples_desc$condition), levels=c(unique(samples_desc$condition)))
  coldata <- as.data.frame(samples_desc[,list(condition = factor(condition,levels = unique(condition)),patient = as.factor(patient))])
  rownames(coldata) <- samples_desc$sample

  if(counts_type == "feature_count"){
    mrcounts <- as.data.frame(fread(counts_file))
    rownames(mrcounts) <- mrcounts$Geneid
  } else {
    load(counts_file)
    colnames(txi$counts) <- make.names(colnames(txi$counts))
    colnames(txi$length) <- make.names(colnames(txi$length))
    colnames(txi$abundance) <- make.names(colnames(txi$abundance))
  }

  ####################################################################################################
  ### Select only desired gene biotypes; replaced parsed Ensembl from Biomart
  # Prepare Ensembl annotation
  library("R.utils")
  if(ref_from_trans_assembly != T){
    if(organism == "homo_sapiens") {
      library("ensembldb")

      EnsDb <- ensembldb::EnsDb(sqlite_path)
      #ensembldb::listColumns(EnsDb.Hsapiens.v87) # See all available columns
      #parsedEnsembl <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(GeneBiotypeFilter(featuresToAnalyze2))) # Get only genes with specificed biotype(s)
      parsedEnsembl <- ensembldb::genes(EnsDb)
      parsedEnsembl <- as.data.frame(parsedEnsembl)
      parsedEnsembl <- parsedEnsembl[,c("gene_id", "gene_name", "gene_biotype")]
    } else {
      parsedEnsembl <- as.data.frame(fread(sqlite_path, sep="\t", header = T))
    }
    parsedEnsembl[is.na(parsedEnsembl$gene_name) | parsedEnsembl$gene_name %in% c("", " "), "gene_name"]<-parsedEnsembl[is.na(parsedEnsembl$gene_name) | parsedEnsembl$gene_name %in% c("", " "), "gene_id"] # Replace missing gene names by gene ids
    rownames(parsedEnsembl)<-parsedEnsembl$gene_id
  } else {
    if(counts_type == "feature_count"){
      parsedEnsembl <- data.frame(gene_id = rownames(mrcounts),gene_name = rownames(mrcounts),gene_biotype = "protein_coding",row.names = rownames(mrcounts))
    }else{
      parsedEnsembl <- data.frame(gene_id = rownames(txi$counts),gene_name = rownames(txi$counts),gene_biotype = "protein_coding",row.names = rownames(txi$counts))
    }
  }


  if(counts_type == "feature_count"){
    # Reorder columns according to the description otherwise DESeq2 calls an error, if we read from a sample sheet
    mrcounts<-mrcounts[,match(rownames(coldata), colnames(mrcounts))]

    # Get only selected genes for a biotype
    mrcounts<-mrcounts[rownames(mrcounts) %in% parsedEnsembl[parsedEnsembl$gene_biotype %in% featuresToAnalyze2,"gene_id"],] # Select only those in parsed Ensembl selection
    ####################################################################################################
    # Remove not expressed genes in any sample
    mrcounts<-mrcounts[rowSums(mrcounts)!=0,]
  } else {

    # Get only selected genes for a biotype
    # TODO - find out how to filter txi by one command
    txi$abundance<-txi$abundance[rownames(txi$abundance) %in% parsedEnsembl[parsedEnsembl$gene_biotype %in% featuresToAnalyze2,
                                                                            "gene_id"],] # Select only those in parsed Ensembl selection
    txi$counts<-txi$counts[rownames(txi$counts) %in% parsedEnsembl[parsedEnsembl$gene_biotype %in% featuresToAnalyze2,
                                                                   "gene_id"],] # Select only those in parsed Ensembl selection
    txi$length<-txi$length[rownames(txi$length) %in% parsedEnsembl[parsedEnsembl$gene_biotype %in% featuresToAnalyze2,
                                                                   "gene_id"],] # Select only those in parsed Ensembl selection

    # Keep only non-zero expressed genes
    #   If gene/isoform has only counts 0.x it will be kept in the output but will have 0 counts in final DE tables
    keep<-rowSums(txi$counts)>0
    txi$abundance<-txi$abundance[keep,]
    txi$counts<-txi$counts[keep,]
    txi$length<-txi$length[keep,]

    txi$counts<-txi$counts[,match(rownames(coldata), colnames(txi$counts))]
    txi$abundance<-txi$abundance[,match(rownames(coldata), colnames(txi$abundance))]
    txi$length<-txi$length[,match(rownames(coldata), colnames(txi$length))]

    # Remove "gene:" or "transcript": added by RSEM - doesn't go well when merged with Ensembl or other annotation-based information
    rownames(txi$abundance)<-gsub("^gene:", "", rownames(txi$abundance))
    rownames(txi$counts)<-gsub("^gene:", "", rownames(txi$counts))
    rownames(txi$length)<-gsub("^gene:", "", rownames(txi$length))
    rownames(txi$abundance)<-gsub("^transcript:", "", rownames(txi$abundance))
    rownames(txi$counts)<-gsub("^transcript:", "", rownames(txi$counts))
    rownames(txi$length)<-gsub("^transcript:", "", rownames(txi$length))

    cts <- txi$counts
    txi$length[txi$length == 0] <- 1 # If gene has length 0 replace it with 1 to avoid error later, might be source of bias and/or error; https://support.bioconductor.org/p/84304/
    anyNA(cts) # Should be false
    normMat <- txi$length
    normMat <- as.data.frame(normMat/exp(rowMeans(log(normMat))))
    mrcounts<-txi$counts
  }
  #print(mrcounts)

  ####################################################################################################
  setwd(OUTPUT_DIR)

  ####################################################################################################
  library(ggplot2)
  if(length(unique(conds)) >= 3){
    num.conds <- length(unique(conds))
  }else{
    num.conds <- 3
  }

  cond_colours<-brewer.pal(num.conds, "Set1")[as.factor(conds)]
  names(cond_colours)<-conds
  # cond_shapes<-c("\u25A0","\u25B2","\u25C6","\u25CF","\u25BC","\u25B6","\u25C0","\u25A3","\u25C8","\u25C9","\u25E9","\u25EA")[]
  # # cond_shapes<-c("\u25CF","\u25B2","\u25FE","\u25C6","\u25BC","\u25BA","\u25C4","\u25C9","\u25C8","\u25A3","\u25E9","\u25EA")[]
  # cond_shapes<-cond_shapes[unique(coldata$patient)]
  # names(cond_shapes)<-unique(coldata$patient)

  bp <- ggplot(data.table(sample = colnames(mrcounts), value = colSums(mrcounts), condition=conds), aes(sample,value, fill=condition)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = cond_colours) +
    theme_bw() + theme(legend.position="bottom") +
    xlab("") +
    ylab("") +
    ggtitle("Total Counts") +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_text(data=data.table(sample = colnames(mrcounts), value = colSums(mrcounts), condition=conds),aes(x=sample,y=value,label=value),vjust=-0.2)

  ggsave(filename = "counts_barplot.pdf", plot = bp, width = 7, height = 7, dpi = 600, device = "pdf")

  sink("DESeq2_design_control.txt")
  print(coldata)
  sink()

  ####################################################################################################
  ####################################################################################################
  # DESeq2 part
  # Make the count object, normalise, dispersion and testing
  library("DESeq2")

  if(counts_type == "feature_count"){
    # Design "design = ~0+condition" would be when we compare without a "common" intercept; when we have dependent before and after treatment patients intercept is on place; when
    #   we compare independent groups "~0+" should be on place https://support.bioconductor.org/p/69374/
    if(INTERCEPT==T){
      print("Using intercept in DESeq2")
      if(PAIRED==T){
        print("Using paired design in DESeq2")
        dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~patient+condition)
      }else{
        print("Using simple design in DESeq2")
        dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~condition)
      }
    }else{
      print("Not using intercept in DESeq2")
      if(PAIRED==T){
        print("Using paired design in DESeq2")
        dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~0+patient+condition)
      }else{
        print("Using simple design in DESeq2")
        dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~0+condition)
      }
    }
  } else {
    #RSEM
    if(INTERCEPT==T){
      print("Using intercept in DESeq2")
      if(PAIRED==T){
        print("Using paired design in DESeq2")
        dds<-DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~patient+condition)
      }else{
        print("Using simple design in DESeq2")
        dds<-DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~condition)
      }
    }else{
      print("Not using intercept in DESeq2")
      if(PAIRED==T){
        print("Using paired design in DESeq2")
        dds<-DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~0+patient+condition)
      }else{
        print("Using simple design in DESeq2")
        dds<-DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~0+condition)
      }
    }
  }


  # Make sure the levels are correct
  if(PAIRED==T){
    print("Correcting paired design factors")
    dds$patient<-factor(dds$patient, levels<-levels(coldata$patient))
  }else{
    print("Correcting simple design factors")
    dds$condition<-factor(dds$condition, levels=levels(coldata$condition))# Make
  }

  # Remove very low count genes
  #keep <- rowSums(counts(dds)) >= 10
  #dds <- dds[keep,]

  # The same thing which follows be calculated by >DESeq(dds) instead of three separate commands
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds)
  dds<-nbinomWaldTest(dds)
  cds<-dds

  ####################################################################################################

  pdf(file="DESeq2_disperison_plot.pdf", width=7, height=7)
  plotDispEsts(cds, main="Dispersion Plot")
  dev.off()

  rawcounts<-counts(cds, normalized=FALSE) # Save raw counts
  normcounts<-counts(cds, normalized=TRUE) # Save normalized counts
  log2counts<-log2(normcounts+1) # Save log2 of normalized counts

  vsd<-DESeq2::varianceStabilizingTransformation(cds) # Save counts transformed with variance Stabilizing Transformation
  vstcounts<-assay(vsd)

  if(PAIRED==TRUE){ # If paired - remove batch effect
    vsd_batch <- vsd
    assay(vsd_batch) <- limma::removeBatchEffect(assay(vsd_batch), vsd_batch$patient) # Designed for log transform so we approximate it with vsd (less negative values than with rlog); we might get negative values which we 'round' to
    vstcounts_batch <- assay(vsd_batch)
    log2counts_batch<-limma::removeBatchEffect(log2(normcounts+1), cds$patient) # Save log2 of normalized counts
    rawcounts_batch<-limma::removeBatchEffect(counts(cds, normalized=FALSE), cds$patient) # Save raw counts

  }


  # Normalization check

  rawcountsum <- colSums(rawcounts)
  rawcountsum <- as.data.frame(rawcountsum)
  rawcountsum$Sample <- row.names(rawcountsum)
  rawcountsum <- as.data.table(rawcountsum)
  rawcountsum <- melt(rawcountsum, id.vars = c("Sample"))
  rawcountsum <- rawcountsum[, .(Sample, condition=conds, value)]

  normcountsum <- colSums(normcounts)
  normcountsum <- as.data.frame(normcountsum)
  normcountsum$Sample <- row.names(normcountsum)
  normcountsum <- as.data.table(normcountsum)
  normcountsum <- melt(normcountsum, id.vars = c("Sample"))
  normcountsum <- normcountsum[, .(Sample, condition=conds, value)]

  rawcountsum$type <- "Raw counts"
  normcountsum$type <- "Normalised counts"
  allcountsum <- rbind(rawcountsum,normcountsum)
  fwrite(allcountsum,"pre_post_norm_counts.tsv", sep="\t")

  library(cowplot)
  # rcs <- ggplot(rawcountsum, aes(Sample, value, fill=condition))+
  #   geom_bar(stat = "identity", width = 0.8)+
  #   scale_fill_manual(values = unique(cond_colours))+
  #   theme_bw() +
  #   ylab("") +
  #   xlab("") +
  #   ggtitle("Pre Normalised Counts") +
  #   theme(axis.text.x = element_text(angle = 90))+
  #   theme(plot.title = element_text(face="bold"))
  #
  # ncs <- ggplot(normcountsum, aes(Sample, value, fill=condition))+
  #   geom_bar(stat = "identity", width = 0.8)+
  #   scale_fill_manual(values = unique(cond_colours))+
  #   theme_bw() +
  #   ylab("") +
  #   xlab("") +
  #   ggtitle("Post Normalised Counts") +
  #   theme(axis.text.x = element_text(angle = 90))+
  #   theme(plot.title = element_text(face="bold"))

  # rcs_ncs <- plot_grid(rcs,ncs,nrow=2)

  # now add the title
  if(condition_design == "all"){
    count.title <- "All samples"
  }else{
    count.title <- paste0(condsToCompare[2]," vs ",condsToCompare[1])
  }

  rcs_ncs <- ggplot(allcountsum, aes(Sample, value, fill=condition))+
    geom_bar(stat = "identity", width = 0.8)+
    scale_fill_manual(values = unique(cond_colours))+
    theme_bw() +
    ylab("") +
    xlab("") +
    facet_grid(factor(type, levels=c("Raw counts","Normalised counts")) ~ .) +
    ggtitle(paste0("Pre & Post Normalised Counts\n",count.title)) +
    theme(axis.text.x = element_text(angle = 90))+
    theme(plot.title = element_text(face="bold"))


  # title <- ggdraw() +
  #   draw_label(count.title,
  #     fontface = 'bold',
  #     x = 0,
  #     hjust = 0
  #   ) +
  #   theme(
  #     # add margin on the left of the drawing canvas,
  #     # so title is aligned with left edge of first plot
  #     plot.margin = margin(0, 0, 0, 7)
  #   )
  #
  # # rel_heights values control vertical title margins
  # title_rcs_ncs <- plot_grid(title, rcs_ncs, ncol = 1,rel_heights = c(0.1, 1))

  ggsave(filename = "pre_post_norm_counts.png", rcs_ncs, units = "in", dpi = 200, width = 7, height = 7, device = "png")
  ggsave(filename = "pre_post_norm_counts.svg", rcs_ncs, width = 7, height = 7, device = "svg")
  ggsave(filename = "pre_post_norm_counts.pdf", rcs_ncs, width = 7, height = 7, device = "pdf")

  # Heatmaps
  library("gplots")

  if(condition_design == "all"){
    hm.log.title <- "Sample to Sample Correlation (Log2)"
    hm.vst.title <- "Sample to Sample Correlation (VST)"
    hm.raw.title <- "Sample to Sample Correlation (Raw Counts)"
  }else{
    hm.log.title <- paste0("Sample to Sample Correlation (Log2)\n",condsToCompare[2]," vs ",condsToCompare[1])
    hm.vst.title <- paste0("Sample to Sample Correlation (VST)\n",condsToCompare[2]," vs ",condsToCompare[1])
    hm.raw.title <- paste0("Sample to Sample Correlation (Raw Counts)\n",condsToCompare[2]," vs ",condsToCompare[1])
  }


  pdf(file="heatmaps_samples.pdf")
  heatmap.2(cor(log2counts), trace="none", col=hmcol, main=hm.log.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
  heatmap.2(cor(vstcounts), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
  heatmap.2(cor(rawcounts), trace="none", col=hmcol, main=hm.raw.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
  dev.off()

  png(file="heatmaps_samples_log.png", width = 7, height = 7, unit = "in", res=200)
  heatmap.2(cor(log2counts), trace="none", col=hmcol, main=hm.log.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
  dev.off()
  png(file="heatmaps_samples_vst.png", width = 7, height = 7, unit = "in", res=200)
  heatmap.2(cor(vstcounts), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
  dev.off()

  svg(file="heatmaps_samples_log.svg", width = 7, height = 7)
  heatmap.2(cor(log2counts), trace="none", col=hmcol, main=hm.log.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
  dev.off()
  svg(file="heatmaps_samples_vst.svg", width = 7, height = 7)
  heatmap.2(cor(vstcounts), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
  dev.off()

  if(PAIRED==TRUE){ # If paired - remove batch effect
    print("Plotting sample heatmaps with batch correction as well.")

    if(condition_design == "all"){
      hm.log.title <- "Sample to Sample Correlation (Log2)\nwith a batch effect removed"
      hm.vst.title <- "Sample to Sample Correlation (VST)\nwith a batch effect removed"
      hm.raw.title <- "Sample to Sample Correlation (Raw Counts)\nwith a batch effect removed"
    }else{
      hm.log.title <- paste0("Sample to Sample Correlation (Log2)\n",condsToCompare[2]," vs ",condsToCompare[1],"\nwith a batch effect removed")
      hm.vst.title <- paste0("Sample to Sample Correlation (VST)\n",condsToCompare[2]," vs ",condsToCompare[1],"\nwith a batch effect removed")
      hm.raw.title <- paste0("Sample to Sample Correlation (Raw Counts)\n",condsToCompare[2]," vs ",condsToCompare[1],"\nwith a batch effect removed")
    }

    pdf(file="heatmaps_samples_batch.pdf")
    heatmap.2(cor(log2counts_batch), trace="none", col=hmcol, main=hm.log.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
    heatmap.2(cor(vstcounts_batch), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
    heatmap.2(cor(rawcounts_batch), trace="none", col=hmcol, main=hm.raw.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
    dev.off()

    png(file="heatmaps_samples_log_batch.png", width = 7, height = 7, unit = "in", res=200)
    heatmap.2(cor(log2counts_batch), trace="none", col=hmcol, main=hm.log.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
    dev.off()
    png(file="heatmaps_samples_vst_batch.png", width = 7, height = 7, unit = "in", res=200)
    heatmap.2(cor(vstcounts_batch), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
    dev.off()

    svg(file="heatmaps_samples_log_batch.svg", width = 7, height = 7)
    heatmap.2(cor(log2counts_batch), trace="none", col=hmcol, main=hm.log.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
    dev.off()
    svg(file="heatmaps_samples_vst_batch.svg", width = 7, height = 7)
    heatmap.2(cor(vstcounts_batch), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=cond_colours, margins=c(9.5,9.5))
    dev.off()

  }

  vstcounts<-vstcounts[apply(vstcounts, 1, max) != 0,]
  # pca<-princomp(vstcounts)
  #
  # pdf(file="sample_to_sample_PCA.pdf")
  # # First two components
  # par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  # plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
  # text(pca$loadings, as.vector(colnames(mrcounts)), pos=3)
  # legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  # # Three components
  # par(mfrow=c(1,3), oma=c(2,0,0,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  # plot(pca$loadings[,c(1,2)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC2", xlab="PC1")
  # text(pca$loadings[,c(1,2)], as.vector(colnames(mrcounts)), pos=3)
  # plot(pca$loadings[,c(1,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC1")
  # text(pca$loadings[,c(1,3)], as.vector(colnames(mrcounts)), pos=3)
  # plot(pca$loadings[,c(2,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC2")
  # text(pca$loadings[,c(2,3)], as.vector(colnames(mrcounts)), pos=3)
  # legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  # dev.off()


  library(ggrepel)

  if(condition_design == "all"){
    pca.title <- "PCA (DESeq2 VST)"
  }else{
    pca.title <- paste0("PCA (DESeq2 VST) ",condsToCompare[2]," vs ",condsToCompare[1])
  }

  # calculate the variance for each gene
  rv <- rowVars(vstcounts)
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  # perform a PCA on the data in vstcounts for the selected genes
  pcaData <- prcomp(t((vstcounts)[select,]))
  # the contribution to the total variance for each component
  percentVar <- pcaData$sdev^2 / sum( pcaData$sdev^2 )
  # assembly the data for the plot
  pca <- data.frame(PC1=pcaData$x[,1], PC2=pcaData$x[,2], PC3=pcaData$x[,3])
  pca$sample <- rownames(pca)
  pca <- merge(samples_desc, pca, by="sample")

  pca1 <- ggplot(pca, aes(PC1, PC2, color=condition))
  if(length(unique(pca$patient)) == length(pca$patient)){
    pca1 <- pca1 + geom_point(size=4, alpha=0.8)
  }else{
    pca1 <- pca1 + geom_point(size=4, alpha=0.8, aes(shape = patient))
  }
  pca1 <- pca1 + scale_color_manual(values = unique(cond_colours), name="") +
    theme_bw() +
    #scale_shape_manual(values=cond_shapes) +
    scale_shape_manual(values=1:length(unique(pca$patient))) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    theme(plot.title = element_text(face="bold")) +
    theme(legend.position="bottom") +
    ggtitle(pca.title)

  pca1S <- pca1 + ggrepel::geom_text_repel(aes(PC1, PC2, label = sample), color="black", max.overlaps = length(pca$sample))
  pca1P <- pca1 + ggrepel::geom_text_repel(aes(PC1, PC2, label = replicate), color="black", max.overlaps = length(pca$sample))

  ggsave(filename = "sample_to_sample_PCA.png", pca1S, units = "in", dpi=200, width = 7, height = 7, device="png")
  ggsave(filename = "sample_to_sample_PCA2.png", pca1P, units = "in", dpi=200, width = 7, height = 7, device="png")
  ggsave(filename = "sample_to_sample_PCA.svg", pca1S, width = 7, height = 7, device="svg")
  ggsave(filename = "sample_to_sample_PCA2.svg", pca1P, width = 7, height = 7, device="svg")
  ggsave(filename = "sample_to_sample_PCA.pdf", pca1S, width = 7, height = 7, device=cairo_pdf)
  ggsave(filename = "sample_to_sample_PCA2.pdf", pca1P, width = 7, height = 7, device=cairo_pdf)

  pca1v3 <- ggplot(pca, aes(PC1, PC3, color=condition))
  if(length(unique(pca$patient)) == length(pca$patient)){
    pca1v3 <- pca1v3 + geom_point(size=4, alpha=0.8)
  }else{
    pca1v3 <- pca1v3 + geom_point(size=4, alpha=0.8, aes(shape = patient))
  }
  pca1v3 <- pca1v3 + scale_color_manual(values = unique(cond_colours), name="") +
    theme_bw() +
    scale_shape_manual(values=1:length(unique(pca$patient))) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC3: ",round(percentVar[3] * 100),"% variance")) +
    theme(plot.title = element_text(face="bold")) +
    ggrepel::geom_text_repel(aes(PC1, PC3, label = replicate), color="black", max.overlaps = length(pca$sample)) +
    theme(legend.position = "none") +
    ggtitle(pca.title)

  pca2v3 <- ggplot(pca, aes(PC2, PC3, color=condition))
  if(length(unique(pca$patient)) == length(pca$patient)){
    pca2v3 <- pca2v3 + geom_point(size=4, alpha=0.8)
  }else{
    pca2v3 <- pca2v3 + geom_point(size=4, alpha=0.8, aes(shape = patient))
  }
  pca2v3 <- pca2v3 + scale_color_manual(values = unique(cond_colours), name="") +
    theme_bw() +
    scale_shape_manual(values=1:length(unique(pca$patient))) +
    xlab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    ylab(paste0("PC3: ",round(percentVar[3] * 100),"% variance")) +
    theme(plot.title = element_text(face="bold")) +
    ggrepel::geom_text_repel(aes(PC2, PC3, label = replicate), color="black", max.overlaps = length(pca$sample)) +
    theme(legend.position = "none") +
    ggtitle(pca.title)

  pca1v2 <- pca1P + theme(legend.position = "none")

  pcax3 <- plot_grid(pca1v2, pca1v3, pca2v3, nrow = 1)
  pcax3.legend <- get_legend( pca1 )
  pcax3wlegend <- plot_grid(pcax3, pcax3.legend, ncol = 1, rel_heights = c(1, .1))
  ggsave(filename = "sample_to_sample_PCAx3.pdf", pcax3wlegend, width = 7, height = 7, device=cairo_pdf)

  ## 3D PCA plot
  library(plotly)
  library(htmlwidgets)

  pca3D <- plot_ly(pca, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, colors=unique(cond_colours))
  pca3D <- pca3D %>% layout(scene = list(xaxis = list(title = paste0("PC1: ",round(percentVar[1] * 100,2),"% variance")),
                                         yaxis = list(title = paste0("PC2: ",round(percentVar[2] * 100,2),"% variance")),
                                         zaxis = list(title = paste0("PC3: ",round(percentVar[3] * 100,2),"% variance"))))
  pca3D <- pca3D %>% add_trace(text = pca$sample, type="scatter3d", mode = "markers", hoverinfo = 'text')
  pca3D <- pca3D %>% layout( title = list(text=pca.title, size = 10))

  htmlwidgets::saveWidget(widget=pca3D ,"sample_to_sample_PCA_3D.html")
  fwrite(pca, "sample_to_sample_PCA.tsv", sep="\t")

  # Get PCA with batch effect from DESeq2 results https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples
  if(PAIRED == T){
    library("ggplot2")
    library("ggrepel")

    if(length(condition_design) == 1){
      pca_batch.title <- "PCA (DESeq2 VST) with a batch effect removed"
    }else{
      pca_batch.title <- paste0("PCA (DESeq2 VST) ",condsToCompare[2]," vs ",condsToCompare[1]," with a batch effect removed")
    }

    # calculate the variance for each gene
    rv_batch <- rowVars(vstcounts_batch)
    # select the ntop genes by variance
    select_batch <- order(rv_batch, decreasing=TRUE)[seq_len(min(500, length(rv_batch)))]
    # perform a PCA on the data in vstcounts_batch for the selected genes
    pcaData_batch <- prcomp(t((vstcounts_batch)[select_batch,]))
    # the contribution to the total variance for each component
    percentVar_batch <- pcaData_batch$sdev^2 / sum( pcaData_batch$sdev^2 )
    # assembly the data for the plot
    pca_batch <- data.frame(PC1=pcaData_batch$x[,1], PC2=pcaData_batch$x[,2], PC3=pcaData_batch$x[,3])
    pca_batch$sample <- rownames(pca_batch)
    pca_batch <- merge(samples_desc, pca_batch, by="sample")

    pca1_batch <- ggplot(pca_batch, aes(PC1, PC2, color=condition))
    if(length(unique(pca_batch$patient)) == length(pca_batch$patient)){
      pca1_batch <- pca1_batch + geom_point(size=4, alpha=0.8)
    }else{
      pca1_batch <- pca1_batch + geom_point(size=4, alpha=0.8, aes(shape = patient))
    }
    pca1_batch <- pca1_batch + scale_color_manual(values = unique(cond_colours), name="") +
      theme_bw() +
      scale_shape_manual(values=1:length(unique(pca$patient))) +
      xlab(paste0("PC1: ",round(percentVar_batch[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar_batch[2] * 100),"% variance")) +
      theme(plot.title = element_text(face="bold")) +
      theme(legend.position="bottom") +
      ggtitle(pca_batch.title)

    pca1_batchS <- pca1_batch + ggrepel::geom_text_repel(aes(PC1, PC2, label = sample), color="black", max.overlaps = length(pca_batch$patient))
    pca1_batchP <- pca1_batch + ggrepel::geom_text_repel(aes(PC1, PC2, label = replicate), color="black", max.overlaps = length(pca_batch$patient))

    ggsave(filename = "sample_to_sample_PCA_batch.png", pca1_batchS, units = "in", dpi=200, width = 7, height = 7, device="png")
    ggsave(filename = "sample_to_sample_PCA_batch2.png", pca1_batchP, units = "in", dpi=200, width = 7, height = 7, device="png")

    ggsave(filename = "sample_to_sample_PCA_batch.svg", pca1_batchS, width = 7, height = 7, device="svg")
    ggsave(filename = "sample_to_sample_PCA_batch2.svg", pca1_batchP, width = 7, height = 7, device="svg")

    ggsave(filename = "sample_to_sample_PCA_batch.pdf", pca1_batchS, width = 7, height = 7, device=cairo_pdf)
    ggsave(filename = "sample_to_sample_PCA_batch2.pdf", pca1_batchP, width = 7, height = 7, device=cairo_pdf)

    ## 3D PCA plot
    pca3D_batch <- plot_ly(pca_batch, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, colors=unique(cond_colours))
    pca3D_batch <- pca3D_batch %>% layout(scene = list(xaxis = list(title = paste0("PC1: ",round(percentVar_batch[1] * 100,2),"% variance")),
                                         yaxis = list(title = paste0("PC2: ",round(percentVar_batch[2] * 100,2),"% variance")),
                                         zaxis = list(title = paste0("PC3: ",round(percentVar_batch[3] * 100,2),"% variance"))))
    pca3D_batch <- pca3D_batch %>% add_trace(text = pca_batch$sample, type="scatter3d", mode = "markers", hoverinfo = 'text')
    pca3D_batch <- pca3D_batch %>% layout( title = list(text=pca_batch.title, size = 10))

    htmlwidgets::saveWidget(widget=pca3D_batch ,"sample_to_sample_PCA_3D_batch.html")
    fwrite(pca_batch, "sample_to_sample_PCA_batch.tsv", sep="\t")


    #pcaData2_batch <- plotPCA(vsd_batch, intgroup=c("condition", "patient"), returnData=TRUE)
    #pcaData2_batch$condition <- as.vector(pcaData2_batch$condition)
    #percentVar <- round(100 * attr(pcaData2_batch, "percentVar"))
    #pca2_batch<-ggplot(pcaData2_batch, aes(PC1, PC2, color=condition)) +
    #  geom_point(size=4, aes(shape = patient)) +
    #  scale_color_manual(values = unique(cond_colours), name="") +
    #  theme_bw() +
    #  ggrepel::geom_text_repel(aes(PC1, PC2, label = rownames(pcaData2_batch)), color="black") +
    #  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    #  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #  theme(plot.title = element_text(face="bold")) +
    #  theme(legend.position="bottom") +
    #  ggtitle(pca_batch.title)

    #ggsave(filename = "sample_to_sample_PCA_batch.png", pca1_batch, units = "in", dpi=200, width = 7, height = 7, device="png")
    #ggsave(filename = "sample_to_sample_PCA_batch.pdf", pca2_batch, width = 7, height = 7, device="pdf")

    write.table(x = assay(vsd_batch), file = "norm_counts_batch.tsv", sep = "\t", col.names = NA)
  }

  #pdf(file="contributions_PCA.pdf", width = 7, height = 7)
  #plot(pca, type = "l", main="Principal Component Contributions")
  #dev.off()

  ####################################################################################################
  # DESeq2 results
  #res<-results(dds, contrast=c("condition", levels(conds)[2], levels(conds)[1])) # Extract contrasts of interest
  res<-results(dds, contrast=c("condition", condsToCompare[2], condsToCompare[1])) # Extract contrasts of interest

  # NEW in DESeq2 >=1.6.0
  # Produce shrunk lfc; used to be default, now in separate command to avoid overshinking in some specialized cases https://support.bioconductor.org/p/95695/; https://rdrr.io/bioc/DESeq2/man/lfcShrink.html; https://github.com/Bioconductor-mirror/DESeq2/blob/master/NEWS; http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
  # this removes column lfcSE
  if(INTERCEPT==T){
    print("Using intercept in lfc shrink")
    if(length(grep(paste0("_vs_", condsToCompare[1], "\\b"), resultsNames(dds)))==(length(resultsNames(dds))-1)){ # Scan how many times we see condition 1 (intercept) - MOST LIKELY WILL CAUSE A LOT OF TROUBLES!!!
      res<- lfcShrink(dds, coef=grep(paste0("condition_", condsToCompare[2], "_vs"), resultsNames(dds)), res=res, type="normal")
      # resLFC <- lfcShrink(dds, coef=2, res=res) # Gives different results than contrast=c("condition", condsToCompare[2], condsToCompare[1]
    }else{
      res<- lfcShrink(dds, contrast=c("condition", condsToCompare[2], condsToCompare[1]), res=res, type="normal") # this removes column lfcSE
    }
  }

  res<-res[order(res$padj, res$pvalue),]
  sig<-rownames(res[(abs(res$log2FoldChange) >= LFC_THRESHOLD) & (res$padj < P_THRESHOLD) & !is.na(res$padj),]) # Extract significant results

  res2<-as.data.frame(res)
  res2<-merge(res2, parsedEnsembl, by="row.names", all.x = T) # Add Ensembl genotype info
  rownames(res2)<-res2$Row.names
  res2$Row.names<-NULL

  res2<-merge(res2, counts(dds, normalized=TRUE), by="row.names") # Add norm. counts
  rownames(res2)<-res2$Row.names
  res2$Row.names<-NULL

  res2<-merge(res2, counts(dds, normalized=FALSE), by="row.names") # Add raw counts
  rownames(res2)<-res2$Row.names
  res2$Row.names<-NULL
  res2$gene_id<-NULL
  res2<-res2[with(res2, order(padj, pvalue, abs(log2FoldChange), -baseMean)), ]

  colnames(res2)<-gsub(pattern = ".x", replacement = "_normCounts", fixed = T, colnames(res2))
  colnames(res2)<-gsub(pattern = ".y", replacement = "_rawCounts", fixed = T, colnames(res2))

  # Quick check of DE genes
  tmpMatrix<-matrix(ncol=1, nrow=5)
  rownames(tmpMatrix)<-c("total genes", paste("LFC >= ", round(LFC_THRESHOLD, 2), " (up)", sep=""), paste("LFC <= ", -(round(LFC_THRESHOLD, 2)), " (down)", sep=""), "not de", "low counts")
  tmpMatrix[1,1]<-nrow(res2)
  tmpMatrix[2,1]<-nrow(res2[res2$log2FoldChange >= (LFC_THRESHOLD) & (res2$padj < P_THRESHOLD) & !is.na(res2$padj),])
  tmpMatrix[3,1]<-nrow(res2[(res2$log2FoldChange <= (-LFC_THRESHOLD)) & (res2$padj < P_THRESHOLD) & !is.na(res2$padj),])
  #tmpMatrix[4,1]<-nrow(res2[(res2$log2FoldChange > (-LFC_THRESHOLD) & (res2$log2FoldChange < (LFC_THRESHOLD))) & !is.na(res2$padj),])
  tmpMatrix[4,1]<-nrow(res2[((res2$padj >= P_THRESHOLD) | ((res2$log2FoldChange > (-LFC_THRESHOLD)) & (res2$log2FoldChange < (LFC_THRESHOLD)))) & !is.na(res2$padj),])
  tmpMatrix[5,1]<-sum(is.na(res2$padj))
  tmpMatrix[,1]<-paste(": ", tmpMatrix[,1], ", ",round(tmpMatrix[,1]/((tmpMatrix[1,1]/100)), 1), "%", sep="")

  sink("DESeq2_de_genes_check.txt")
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " Without LogFC Cut-off", sep=""))
  summary(res, alpha=P_THRESHOLD)
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " and LogFC >= ", round(LFC_THRESHOLD, 2), sep=""))
  print(noquote(tmpMatrix))
  print("Used tests")
  write.table(mcols(res)$description, col.names=F, row.names=F)
  sink()

  write.table(x = res2, file = "DESeq2.tsv", sep = "\t", col.names = NA)

  # Results without independent filtering and cooks cutoff - raw results
  #resNoFil<-results(dds, contrast=c("condition", levels(conds)[2], levels(conds)[1]), independentFiltering=F, cooksCutoff=F) # Turn off independent filtering and cooks cutoff for outliers
  resNoFil<-results(dds, contrast=c("condition", condsToCompare[2], condsToCompare[1]), independentFiltering=F, cooksCutoff=F) # Turn off independent filtering and cooks cutoff for outliers

  # NEW in DESeq2 >=1.6.0
  # Produce shrunk lfc; used to be default, now in separate command to avoid overshinking in some specialized cases https://support.bioconductor.org/p/95695/; https://rdrr.io/bioc/DESeq2/man/lfcShrink.html; https://github.com/Bioconductor-mirror/DESeq2/blob/master/NEWS; http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
  # this removes column lfcSE
  if(INTERCEPT==T){
    print("Using intercept in lfc shrink")
    if(length(grep(paste0("_vs_", condsToCompare[1], "\\b"), resultsNames(dds)))==(length(resultsNames(dds))-1)){ # Scan how many times we see condition 1 (intercept) - MOST LIKELY WILL CAUSE A LOT OF TROUBLES!!!
      resNoFil<- lfcShrink(dds, coef=grep(paste0("condition_", condsToCompare[2], "_vs"), resultsNames(dds)), res=resNoFil, type="normal")
      # resNoFilLFC <- lfcShrink(dds, coef=2, res=resNoFil) # Gives different results than contrast=c("condition", condsToCompare[2], condsToCompare[1]
    }else{
      resNoFil<- lfcShrink(dds, contrast=c("condition", condsToCompare[2], condsToCompare[1]), res=resNoFil, type="normal") # this removes column lfcSE
    }
  }

  resNoFil<-resNoFil[order(resNoFil$padj, resNoFil$pvalue),]
  sigNoFil<-rownames(resNoFil[(abs(resNoFil$log2FoldChange) >= LFC_THRESHOLD) & (resNoFil$padj < P_THRESHOLD) & !is.na(resNoFil$padj),])

  resNoFil2<-as.data.frame(resNoFil)
  resNoFil2<-merge(resNoFil2, parsedEnsembl, by="row.names", all.x = T)
  rownames(resNoFil2)<-resNoFil2$Row.names
  resNoFil2$Row.names<-NULL

  resNoFil2<-merge(resNoFil2, counts(dds, normalized=TRUE), by="row.names")
  rownames(resNoFil2)<-resNoFil2$Row.names
  resNoFil2$Row.names<-NULL

  resNoFil2<-merge(resNoFil2, counts(dds, normalized=FALSE), by="row.names")
  rownames(resNoFil2)<-resNoFil2$Row.names
  resNoFil2$Row.names<-NULL
  resNoFil2$gene_id<-NULL
  resNoFil2<-resNoFil2[with(resNoFil2, order(padj, pvalue, abs(log2FoldChange), -baseMean)), ]

  colnames(resNoFil2)<-gsub(pattern = ".x", replacement = "_normCounts", fixed = T, colnames(resNoFil2))
  colnames(resNoFil2)<-gsub(pattern = ".y", replacement = "_rawCounts", fixed = T, colnames(resNoFil2))

  # Quick check of DE genes
  tmpMatrix<-matrix(ncol=1, nrow=5)
  rownames(tmpMatrix)<-c("total genes", paste("LFC >= ", round(LFC_THRESHOLD, 2), " (up)", sep=""), paste("LFC <= ", -(round(LFC_THRESHOLD, 2)), " (down)", sep=""), "not de", "low counts")
  tmpMatrix[1,1]<-nrow(resNoFil2)
  tmpMatrix[2,1]<-nrow(resNoFil2[resNoFil2$log2FoldChange >= (LFC_THRESHOLD) & (resNoFil2$padj < P_THRESHOLD) & !is.na(resNoFil2$padj),])
  tmpMatrix[3,1]<-nrow(resNoFil2[(resNoFil2$log2FoldChange <= (-LFC_THRESHOLD)) & (resNoFil2$padj < P_THRESHOLD) & !is.na(resNoFil2$padj),])
  #tmpMatrix[4,1]<-nrow(resNoFil2[(resNoFil2$log2FoldChange > (-LFC_THRESHOLD) & (resNoFil2$log2FoldChange < (LFC_THRESHOLD))) & !is.na(resNoFil2$padj),])
  tmpMatrix[4,1]<-nrow(resNoFil2[((resNoFil2$padj >= P_THRESHOLD) | ((resNoFil2$log2FoldChange > (-LFC_THRESHOLD)) & (resNoFil2$log2FoldChange < (LFC_THRESHOLD)))) & !is.na(resNoFil2$padj),])
  tmpMatrix[5,1]<-0 # Turned off by turning off cooks cutoff and ind. filtering
  tmpMatrix[,1]<-paste(": ", tmpMatrix[,1], ", ",round(tmpMatrix[,1]/((tmpMatrix[1,1]/100)), 1), "%", sep="")

  sink("DESeq2_de_genes_check_noIndFilt.txt")
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " Without LogFC Cut-off", sep=""))
  summary(resNoFil, alpha=P_THRESHOLD)
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " and LogFC >= ", round(LFC_THRESHOLD, 2), sep=""))
  print(noquote(tmpMatrix))
  print("Used tests")
  write.table(mcols(resNoFil)$description, col.names=F, row.names=F)
  sink()

  write.table(x = resNoFil2, file ="DESeq2_noIndFilt.tsv", sep = "\t", col.names = NA)

  # sink("used_test.txt") # Save used tests and groups
  #   print(mcols(res)$description)
  # sink()

  ####################################################################################################
  TOP_BCKP<-TOP

  if(length(sig)<TOP){ # Avoid error by ploting more TOP genes than significantly DE genes
    TOP<-length(sig)
  }

  if(TOP==0){ # Set range for naming the samples - help to avoid naming of samples even if there is no DE gene; THIS ERROR IS OK WHEN WE HAVE 0 DE GENES Error in text.default(res[RANGE, ]$log2FoldChange, -log(res[RANGE, ]$padj,  : zero-length 'labels' specified
    RANGE<-0
    print("ZERO DE GENES FOUND, expecting error \"Error in text.default(res[RANGE, ]$log2FoldChange, -log(res[RANGE, ]$padj,  : zero-length \'labels\' specified\"!")
  }else{
    if(TOP==1){ # Helps to avoid naming if ONE DE is called
      RANGE<-1
    }else{
      RANGE<-c(1:TOP) # Helps to avoid naming if MORE THAN ONE DE is called
    }
  }

  res<-res[order(res$padj, res$pvalue),] # Make sure res is ordered by adj.p-value

  # Add custom gene names to volcano plot
  #selectedGenes<-which(rownames(res) %in% parsedEnsembl[parsedEnsembl$gene_name %in% c("IL2RA", "IL2RB"), "gene_id"])
  #RANGE<-c(RANGE, selectedGenes)
  if(TOP>0){

    pdf(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1],".pdf", sep=""))
    par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
    plot(res$log2FoldChange, -log(res$padj, 10), main=paste("Volcanoplot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""), cex=0.4, pch=19)
    #   text(res[1:length(sig),]$log2FoldChange, -log(res[1:length(sig),]$padj,10),
    #        labels=parsedEnsembl[parsedEnsembl$gene_id %in% rownames(res[1:length(sig),]), "gene_name"], cex=0.3, pos=3)
    text(res[RANGE,]$log2FoldChange, -log(res[RANGE,]$padj,10),
         labels=parsedEnsembl[rownames(res)[RANGE], "gene_name"], cex=0.3, pos=3)
    legend("bottomleft", condsToCompare[1], cex=0.5)
    legend("bottomright", condsToCompare[2], cex=0.5)
    abline(h=-log(P_THRESHOLD, 10), col="red", lty=2)
    abline(v=c(-LFC_THRESHOLD, LFC_THRESHOLD), col="darkblue", lty=2)
    dev.off()

  }

  pdf(file=paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],".pdf", sep=""))
  par(xpd=TRUE) # Allow labels to be out of the frame
  DESeq2::plotMA(res, alpha=P_THRESHOLD, main=paste("MA plot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""),
                 ylim=c(-max(abs(res$log2FoldChange), na.rm=T), max(abs(res$log2FoldChange), na.rm=T)))
  abline(h=c(-LFC_THRESHOLD, LFC_THRESHOLD), col="darkblue", lty=2)
  text(res[1:TOP,]$baseMean, res[1:TOP,]$log2FoldChange, labels=parsedEnsembl[rownames(res)[1:TOP], "gene_name"], cex=0.3, pos=3)
  legend("bottomright", condsToCompare[1], cex=0.5)
  legend("topright", condsToCompare[2], cex=0.5)
  dev.off()

  # Plot Volcano plot with ggplot2 http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html
  library("dplyr")
  library("ggplot2")
  library("ggrepel")
  resForPlot<-as.data.frame(res)
  resForPlot$Gene<-substr(parsedEnsembl[rownames(res), "gene_name"],1,20)
  results <- mutate(resForPlot, sig=ifelse(resForPlot$padj<P_THRESHOLD, paste0("padj<", P_THRESHOLD), "Not Sig"))
  volkan <- ggplot(results, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(col=sig), size=0.5) +
    scale_color_manual(values=c("black", "red"))

  results<-results[order(results$padj, results$pvalue),] # make sure it is correctly ordered

  volkan <- volkan + geom_text_repel(data=dplyr::filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)+geom_vline(xintercept = 0)+
    geom_vline(xintercept = c(LFC_THRESHOLD, -LFC_THRESHOLD), linetype = "longdash", colour="blue")+
    ggtitle(paste("Volcanoplot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""))+
    theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + theme(plot.title = element_text(face="bold"))

  pdf(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggplot2.pdf", sep=""))
  print(volkan)
  dev.off()

  png(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggplot2.png", sep=""), units = "in",width = 7, height = 7, res = 200)
  print(volkan)
  dev.off()

  svg(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggplot2.svg", sep=""),width = 7, height = 7)
  print(volkan)
  dev.off()

  #ggsave(paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggplot2.svg", sep=""), p, width = 7, height = 7, device = svg)

  # Great MA plot http://www.sthda.com/english/rpkgs/ggpubr/reference/ggmaplot.html
  library("ggpubr")

  #ggmaplot(resForPlot, main =  expression("Group 1" %->% "Group 2"),
  ma <- ggmaplot(resForPlot, main =  paste0("MA plot ", condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes"),
                fdr = P_THRESHOLD, fc = FOLD_CHANGE, size = 0.4,
                palette = c("#B31B21", "#1465AC", "darkgray"),
                genenames = as.vector(resForPlot$Gene),
                legend = "top", top = TOP,
                font.label = c("bold", 11),
                font.legend = "bold",
                font.main = "bold",
                ggtheme = ggplot2::theme_minimal())+
    theme(plot.title = element_text(hjust = 0.5))

  pdf(file=paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggpubr.pdf", sep=""))
  print(ma)
  dev.off()

  png(file=paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggpubr.png", sep=""), units = "in",width = 7, height = 7, res = 200)
  print(ma)
  dev.off()

  svg(file=paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggpubr.svg", sep=""),width = 7, height = 7)
  print(ma)
  dev.off()

  #ggsave(paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggpubr.svg", sep=""), ma, width = 7, height = 7, device=svg)

  TOP<-TOP_BCKP

  if(length(sigNoFil)<TOP){  # Avoid error by ploting more TOP genes than significantly DE genes
    TOP<-length(sigNoFil)
  }

  if(TOP==0){ # Set range for naming the samples - help to avoid naming of samples even if there is no DE gene
    RANGE<-0
    print("ZERO DE GENES FOUND, expecting error \"Error in text.default(res[RANGE, ]$log2FoldChange, -log(res[RANGE, ]$padj,  : zero-length \'labels\' specified\"!")
  }else{
    if(TOP==1){ # Helps to avoid naming if ONE DE is called
      RANGE<-1
    }else{
      RANGE<-c(1:TOP) # Helps to avoid naming if MORE THAN ONE DE is called
    }
  }

  resNoFil<-resNoFil[order(resNoFil$padj, resNoFil$pvalue),] # Make sure resNoFil is ordered by adj.p-value

  if(TOP>0){
    pdf(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1], "_noIndFilt.pdf", sep=""))
    par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
    plot(resNoFil$log2FoldChange,-log(resNoFil$padj,10),main=paste("Volcanoplot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""), cex=0.4, pch=19)
    #   text(resNoFil[1:length(sigNoFil),]$log2FoldChange, -log(resNoFil[1:length(sigNoFil),]$padj,10),
    #        labels=parsedEnsembl[parsedEnsembl$gene_id %in% rownames(resNoFil[1:length(sigNoFil),]), "gene_name"], cex=0.3, pos=3)
    text(resNoFil[RANGE,]$log2FoldChange, -log(resNoFil[RANGE,]$padj,10),
         labels=parsedEnsembl[rownames(resNoFil)[RANGE], "gene_name"], cex=0.3, pos=3)
    legend("bottomleft", condsToCompare[1], cex=0.5)
    legend("bottomright", condsToCompare[2], cex=0.5)
    abline(h=-log(P_THRESHOLD, 10), col="red", lty=2)
    abline(v=c(-LFC_THRESHOLD, LFC_THRESHOLD), col="darkblue", lty=2)
    dev.off()

    pdf(file=paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_noIndFilt.pdf", sep=""))
    par(xpd=TRUE) # Allow labels to be out of the frame
    DESeq2::plotMA(res, alpha=P_THRESHOLD, main=paste("MA plot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""),
                   ylim=c(-max(abs(res$log2FoldChange), na.rm=T), max(abs(res$log2FoldChange), na.rm=T)))
    abline(h=c(-LFC_THRESHOLD, LFC_THRESHOLD), col="darkblue", lty=2)
    text(resNoFil[1:TOP,]$baseMean, resNoFil[1:TOP,]$log2FoldChange, labels=parsedEnsembl[rownames(res)[1:TOP], "gene_name"], cex=0.3, pos=3)
    legend("bottomright", condsToCompare[1], cex=0.5)
    legend("topright", condsToCompare[2], cex=0.5)
    dev.off()
  }
  # Plot Volcano plot with ggplot2 http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html
  # library(dplyr)
  # library(ggplot2)
  # library(ggrepel)
  resForPlot<-as.data.frame(resNoFil)
  resForPlot$Gene<-substr(parsedEnsembl[rownames(resNoFil), "gene_name"],1,20)
  results <- mutate(resForPlot, sig=ifelse(resForPlot$padj<P_THRESHOLD, paste0("padj<", P_THRESHOLD), "Not Sig"))
  volkan_nofilt <- ggplot(results, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(col=sig), size=0.5) +
    scale_color_manual(values=c("black", "red"))
  #p
  results<-results[order(results$padj, results$pvalue),] # make sure it is correctly ordered
  #p+geom_text(data=filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)

  pdf(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1], "_noIndFilt_ggplot2.pdf", sep=""))
  volkan_nofilt <- volkan_nofilt + geom_text_repel(data=dplyr::filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)+geom_vline(xintercept = 0)+
    geom_vline(xintercept = c(LFC_THRESHOLD, -LFC_THRESHOLD), linetype = "longdash", colour="blue")+
    ggtitle(paste("Volcanoplot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""))+
    theme(plot.title = element_text(hjust = 0.5))
  print(volkan_nofilt)
  dev.off()

  # Trying MA plot with ggplot2 http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html
  # p = ggplot(results, aes(baseMean, log2FoldChange)) +
  #   geom_point(aes(col=sig), size=0.5) +
  #   scale_color_manual(values=c("black", "red"))
  # p
  # results<-results[order(results$padj, results$pvalue),] # make sure it is correctly ordered
  # p+geom_text(data=filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)
  # p+geom_text_repel(data=filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)+geom_vline(xintercept = 0)+
  #   geom_vline(xintercept = c(LFC_THRESHOLD, -LFC_THRESHOLD), linetype = "longdash", colour="blue")+
  #   ggtitle(paste("MA plot ", condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""))+
  #   theme(plot.title = element_text(hjust = 0.5))

  # Great MA plot http://www.sthda.com/english/rpkgs/ggpubr/reference/ggmaplot.html
  # library("ggpubr")

  pdf(file=paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_noIndFilt_ggpubr.pdf", sep=""))
  #ggmaplot(resForPlot, main =  expression("Group 1" %->% "Group 2"),
  ma_nofilt <- ggmaplot(resForPlot, main =  paste0("MA plot ", condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes"),
                fdr = P_THRESHOLD, fc = FOLD_CHANGE, size = 0.4,
                palette = c("#B31B21", "#1465AC", "darkgray"),
                genenames = as.vector(resForPlot$Gene),
                legend = "top", top = TOP,
                font.label = c("bold", 11),
                font.legend = "bold",
                font.main = "bold",
                ggtheme = ggplot2::theme_minimal())+
    theme(plot.title = element_text(hjust = 0.5))
  print(ma_nofilt)
  dev.off()

  TOP<-TOP_BCKP

  # Heatmaps of selected genes
  library("pheatmap")

  res_tmp<-res[order(res$padj, res$pvalue, -abs(res$log2FoldChange), -(res$baseMean)),]
  select<-rownames(res_tmp[(abs(res_tmp$log2FoldChange) >= LFC_THRESHOLD) & (res_tmp$padj < P_THRESHOLD) & !is.na(res_tmp$padj),])
  selectAll<-select

  if(length(select)>TOP){
    select<-select[1:TOP]
  }

  # Add custom genes to the heatmap
  #select <- c(select, rownames(resForPlot)[resForPlot$Gene %in% c("IL2RA", "IL2RB")])

  # save(dds,select,file = "test_file.Rdata")
  nt<-DESeq2::normTransform(dds) # defaults to log2(x+1)
  # nt<-DESeq2::rlog(dds, blind=FALSE)
  # nt<-DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

  if(PAIRED==TRUE){
    print("Correcting the expressions for batch effect for the heatmaps")
    assay(nt) <- limma::removeBatchEffect(assay(nt), nt$patient)
  }

  log2.norm.counts<-assay(nt)[select,,drop = F]
  # rownames(log2.norm.counts)<-parsedEnsembl[rownames(log2.norm.counts), "gene_name"]
  rownames(log2.norm.counts)<-substr(parsedEnsembl[rownames(log2.norm.counts), "gene_name"],1,20)
  log2.norm.counts<-log2.norm.counts[order(rowMeans(log2.norm.counts)),]

  if(PAIRED==T){
    print("Using paired design")
    df<-as.data.frame(colData(dds)[,c("condition","patient")])
  }else{
    print("Using simple design")
    df<-as.data.frame(colData(dds)[,c("condition")])
  }

  TOP_BCKP<-TOP

  if(TOP>length(select)){
    TOP<-length(select)
  }

  rownames(df)<-colnames(log2.norm.counts) # Rename rows to fit to column names

  if(ncol(df)==1){ # Change naming of the df for pheatmap - nicer visualization
    colnames(df)<-"condition"
  }

  if(TOP >  1){
    #pdf(file="heatmap_selected_orderBaseMeanCluster.pdf", onefile=FALSE, height= 2 + (TOP/2))
    pdf(file="heatmap_selected_orderBaseMeanCluster.pdf", onefile=FALSE, height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[2], " vs ", condsToCompare[1], sep=""))
    dev.off()

    png(filename = "heatmap_selected_orderBaseMeanCluster.png", res = 200, units = "in", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[2], " vs ", condsToCompare[1], sep=""))
    dev.off()

    svg(filename = "heatmap_selected_orderBaseMeanCluster.svg", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[2], " vs ", condsToCompare[1], sep=""))
    dev.off()

    pdf(file="heatmap_selected_orderBaseMean.pdf", onefile=FALSE, height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[2], " vs ", condsToCompare[1], sep=""))
    dev.off()

    png(filename = "heatmap_selected_orderBaseMean.png", res = 200, units = "in", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[2], " vs ", condsToCompare[1], sep=""))
    dev.off()

    svg(filename = "heatmap_selected_orderBaseMean.svg", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[2], " vs ", condsToCompare[1], sep=""))
    dev.off()

  }

  TOP<-TOP_BCKP

  graphics.off()

  #system("for i in heatmap_selected_*; do pdftk $i cat 2-end output tmp.pdf; mv tmp.pdf $i; done") # Cut first empty page from heatmap plot
  tryCatch({
    if(TOP > 1){
      # Plot counts for all significant genes
      pdf(file="all_sig_genes_normCounts.pdf")
      for(i in 1:length(selectAll)){
        plotCounts(dds, gene=selectAll[i], intgroup="condition", main=parsedEnsembl[selectAll[i], "gene_name"])
        mtext(paste0("adj. p-value < ", P_THRESHOLD, " logFC >= ", round(LFC_THRESHOLD,3)))
        # axis(1, at=seq_along(levels(coldata$condition)), levels(coldata$condition), las=2) # Ugly but works; I am not able to turn off axis() setting in plotCounts function
      }
      dev.off()
    }
  }, error=function(e){})
  # Write all normalized counts
  tmp.table<-as.data.frame(normcounts)
  tmp.table$gene_id<-rownames(tmp.table)
  tmp.table<-merge(tmp.table, parsedEnsembl, by="gene_id", all.x=T)
  write.table(tmp.table, file="norm_counts.tsv", sep="\t", row.names=T, col.names=NA)
  write.table(rownames(res), file="background.txt", row.names=F, quote = F) # Write background file

  # Make gene background/universe
  background<-as.data.frame(res2[, "gene_name"])
  colnames(background)<-"gene_name"
  background$gene_id<-rownames(res2)
  write.table(background, file="background.txt", row.names=F, quote = F, sep="\t") # Write background file
  ####################################################################################################
  ####################################################################################################
  # edgeR part
  library("edgeR")

  if(counts_type == "feature_count"){
    d<-DGEList(counts=mrcounts, group=coldata$condition) # edgeR DGE object
    d<-calcNormFactors(d) # Calculate normalization factors

    # Filtering of non-informative genes - very low count-per-million; recommended by edgeR vignette
    # http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf -
    # 4.1.3 Filter low expression tags
    keep<-rowSums(cpm(d)>1) >= 3 # Change to 1/3 of the samples? (length(colnames(d))/3)
    d<-d[keep, , keep.lib.sizes=FALSE]
    d<-calcNormFactors(d) # Re-calculate normalization factors after filtering
  } else {
    o <- log(calcNormFactors(cts/normMat, na.rm = T)) + log(colSums(cts/normMat, na.rm = T)) # Must not contain NA values
    d <- DGEList(cts, group=coldata$condition)
    d$offset <- t(t(log(normMat)) + o) # d is now ready for estimate dispersion functions see edgeR User's Guide
    #d<-calcNormFactors(d) # Re-calculate normalization factors after filtering - should be included in offset! Also in tximport manual they say "...is now ready for estimate dispersion functions see edgeR User's Guide..."

    # Filtering of non-informative genes - very low count-per-million; recommended by edgeR vignette
    # http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf -
    # 4.1.3 Filter low expression tags
    keep<-rowSums(cpm(d)>1) >= 3 # Change to 1/3 of the samples? (length(colnames(d))/3)
    d<-d[keep, , keep.lib.sizes=FALSE]
    # Recalculate cts and normMat to exclude removed genes
    cts<-cts[rownames(cts) %in% rownames(d$counts), ]
    normMat<-txi$length[rownames(txi$length) %in% rownames(d$counts),]
    normMat <- as.data.frame(normMat/exp(rowMeans(log(normMat))))
    anyNA(normMat)
    o <- log(calcNormFactors(cts/normMat, na.rm = T)) + log(colSums(cts/normMat, na.rm = T)) # Must not contain NA values
    d <- DGEList(cts, group=coldata$condition)
    d$offset <- t(t(log(normMat)) + o)
  }

  # Design "design = ~0+condition" would be when we compare without a "common" intercept; when we have dependent before and after treatment patients intercept is on place; when
  #   we compare independent groups "~0+" should be on place https://support.bioconductor.org/p/69374/
  #design<-model.matrix(~coldata$patient+coldata$condition)
  #design<-model.matrix(+coldata$condition)
  #colnames(design) <- gsub("coldata$condition", "", colnames(design), fixed = T)

  if(INTERCEPT==T){
    print("Using intercept in edgeR")
    if(PAIRED==T){
      print("Using paired design in edgeR")
      design<-model.matrix(~patient+condition, data=coldata)
    }else{
      print("Using simple design in edgeR")
      design<-model.matrix(~condition, data=coldata)
    }
  }else{
    print("Not using intercept in edgeR")
    if(PAIRED==T){
      print("Using paired design in edgeRs")
      design<-model.matrix(~0+patient+condition, data=coldata)
    }else{
      print("Using simple design in edgeR")
      design<-model.matrix(~0+condition, data=coldata)
    }
  }

  # The same thing which follows can be done using >d<-estimateDisp(d, design) which should be also used if you plan to use QL methods bellow
  # https://support.bioconductor.org/p/79149/
  # You can use  estimateGLMRobustDisp() is you expect outliers in the data not associated with any sample or particular gene
  d<-estimateGLMCommonDisp(d, design) # Calculate GLM for common dispersion
  d<-estimateGLMTrendedDisp(d, design) # Calculate GLM for trended dispersion
  d<-estimateGLMTagwiseDisp(d, design) # Calculate GLM for tagwise dispersion

  # https://support.bioconductor.org/p/76790/
  fit_tgw<-glmFit(d, design, dispersion=d$tagwise.dispersion) # Fit tagwise dispersion;  fit_tgw<-glmQLFit(d, design, dispersion=d$tagwise.dispersion) can be used if the number of replicates is low; QL (glmQLFit and glmQLFTest) is more strict in the assumptions ~ increases adj.pvalues

  ### Replaced by following generation of contrast
  # if(INTERCEPT==T){
  # 	lrt_tgw<-glmLRT(fit_tgw, coef=grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit_tgw$design))) # coef=2 should be the same as contrast=c(0, -1, 0) if we want to compare intercept and second condition https://www.biostars.org/p/102036/
  # #	lrt_tgw<-glmLRT(fit_tgw, contrast=c(0, -1, 1)) # If we want to compare other two (and more) conditions without the intercept
  # }else{ # Intercepted design; if we want to compare 3 vs 1 or 2 vs 1 we just set coef=3 or coef=2, respectively; fit_tgw<-glmQLFTest(fit_tgw, coef=ncol(design)) can be used if the number of replicates is low; we can also make a bit more precise calculation for genes with certain logFC - we set logFC under which genes ARE NOT interesting for us (not the logFC we expect to see!) - >lrt_tgw<-glmTreat(fit_tgw, coef=ncol(design), lfc=0.5)
  # # In case we have more comparisons and/or non-intercepted design
  # 	my.contrasts <- makeContrasts(
  #   	postvspre = paste0(colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[1], "\\b"), colnames(fit_tgw$design))], "-", colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit_tgw$design))]), levels=fit_tgw$design) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
  # 	colnames(my.contrasts)<-"contrast"
  # 	lrt_tgw<-glmLRT(fit_tgw, contrast=my.contrasts[, "contrast"]) # If we have 3 conditions and want to compare 3 vs 1 we set contrast=c(0, -1, 1), if we want to compare 3 vs 1 or 2 vs 1 we just set coef=3 or coef=2, respectively; some more examples of contrast https://www.biostars.org/p/110861/
  # }
  ### Replacement of coef and contrast, should do the same
  if(!length(grep(paste0("condition", condsToCompare[1], "\\b"), colnames(fit_tgw$design)))){# If I cannot find condsToCompare[1] (usually intercept) set contrast as 1 for condsToCompare2
    cond_label <<- paste0(colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit_tgw$design))])
    my.contrasts <- makeContrasts(postvspre = cond_label, levels=fit_tgw$design) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
    # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
  }else{ # Else make proper contrast
    cond_label <<- paste0(colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit_tgw$design))], "-", colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[1], "\\b"), colnames(fit_tgw$design))])
    # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
    my.contrasts <- makeContrasts(postvspre = cond_label, levels=fit_tgw$design) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
  }
  colnames(my.contrasts)<-"contrast"
  lrt_tgw<-glmLRT(fit_tgw, contrast=my.contrasts[, "contrast"]) # If we have 3 conditions and want to compare 3 vs 1 we set contrast=c(0, -1, 1), if we want to compare 3 vs 1 or 2 vs 1 we just set coef=3 or coef=2, respectively; some more examples of contrast https://www.biostars.org/p/110861/

  resultsTbl.tgw<-topTags(lrt_tgw, n=nrow(lrt_tgw$table), adjust.method = "BH", sort.by = "p.value")$table # Extract all genes

  ####################################################################################################
  # Store row number that match between results and raw counts to include raw counts later in the
  #   output
  wh.rows.tgw<-match(rownames(resultsTbl.tgw), rownames(d$counts))

  # Combine results with extracted DE genes, common dispersion, UpDown values, normalized counts and
  #   raw counts

  decidetesttab <- decideTestsDGE(lrt_tgw,adjust.method="BH", p.value=P_THRESHOLD, lfc=LFC_THRESHOLD)
  decidetesttab <- decidetesttab@.Data[wh.rows.tgw]
  combResults.tgw<-cbind(resultsTbl.tgw,
                         "tgw.Disp"=d$tagwise.dispersion[wh.rows.tgw],
                         "UpDown"=decidetesttab,
                         cpm(d$counts[wh.rows.tgw,]),
                         d$counts[wh.rows.tgw,])

  # For all dispersions modify the results with renaming the columns for norm and raw counts
  reformateFinalTable<-function(inputTable, ...){
    colnames(inputTable)[8:((nrow(coldata))+7)]<-paste(colnames(inputTable)[8:((nrow(coldata))+7)], 'normCounts', sep="_")
    colnames(inputTable)[((nrow(coldata))+8):(((nrow(coldata))+7)+(nrow(coldata)))]<-paste(colnames(inputTable)[((nrow(coldata))+8):(((nrow(coldata))+7)+(nrow(coldata)))],'rawCounts',sep="_")
    inputTable<-inputTable[order(inputTable$FDR, inputTable$PValue, -abs(inputTable$logFC), -inputTable$logCPM),]
    return(inputTable)
  }

  combResults.tgw<-reformateFinalTable(combResults.tgw)

  # Merge with gene names
  combResults.tgw<-merge(combResults.tgw, parsedEnsembl, by="row.names", all.x = T) # Add Ensembl genotype info
  rownames(combResults.tgw)<-combResults.tgw$Row.names
  combResults.tgw$Row.names<-NULL
  combResults.tgw$gene_id<-NULL

  # Reformate the final table a bit
  columnsTmp<-c("logCPM", "logFC", "PValue", "FDR", "LR", "tgw.Disp", "UpDown", "gene_name", "gene_biotype")
  columnsTmp2<-colnames(combResults.tgw)[!(colnames(combResults.tgw) %in% columnsTmp)]
  combResults.tgw<-combResults.tgw[c(columnsTmp, columnsTmp2)]
  combResults.tgw<-combResults.tgw[order(combResults.tgw$FDR, combResults.tgw$PValue, -abs(combResults.tgw$logFC), -combResults.tgw$logCPM),]

  colnames(combResults.tgw)[colnames(combResults.tgw)=="logFC"]<-"log2FoldChange"
  colnames(combResults.tgw)[colnames(combResults.tgw)=="PValue"]<-"pvalue"
  colnames(combResults.tgw)[colnames(combResults.tgw)=="FDR"]<-"padj"
  combResults.tgw$UpDown<-NULL

  write.table(x = combResults.tgw, file = "edgeR.tsv", sep="\t", row.names=T, col.names=NA)

  ####################################################################################################

  sink("edgeR_design_control.txt") #save sample names to control with lib.sizes
  print(subset(cbind(coldata, d$samples), select=-c(group)))
  sink()

  points<-c(0, 1, 2, 5, 6, 15, 16, 17, 18)

  if(length(unique(d$samples$group)) >= 3){
    num.conds <- length(unique(conds))
  }else{
    num.conds <- 3
  }

  cond_colours<-brewer.pal(num.conds, "Paired")[d$samples$group]
  names(cond_colours)<-d$samples$group

  dfMDS<-data.frame(x=plotMDS(d,method="bcv")$x, y=plotMDS(d,method="bcv")$y)
  dfMDS$sample <- rownames(dfMDS)
  dfMDS$condition <- conds

  mds <- ggplot(dfMDS, aes(x, y, color=condition)) +
    geom_point(size = 3) +
    scale_color_manual(values = unique(cond_colours), name="") +
    theme_bw() +
    xlab("BCV distance 1") +
    ylab("BCV distance 2") +
    theme(legend.position="bottom") +
    ggtitle("MDSPlot_BCV_distance") +
    ggrepel::geom_text_repel(aes(x, y, label = sample), color="black") +
    theme(plot.title = element_text(face="bold"))

  ggsave("MDS_plot.pdf", mds, units = "in", width = 7, height = 7, dpi = 200)

  # See the effect of batch - should be included in the formula
  A<-aveLogCPM(d)
  d2<-d[A>1,]
  d2<-calcNormFactors(d2)
  logCPM<-cpm(d2, log=TRUE, prior.count=5)
  logCPMc<-removeBatchEffect(logCPM, coldata$patient)

  dflogCPM<-data.frame(x=plotMDS(logCPM)$x, y=plotMDS(logCPM)$y)
  dflogCPM$sample <- rownames(dflogCPM)
  dflogCPM$condition <- conds

  dflogCPMc<-data.frame(x=plotMDS(logCPMc)$x, y=plotMDS(logCPMc)$y)
  dflogCPMc$sample <- rownames(dflogCPMc)
  dflogCPMc$condition <- conds

  mds.logcpm <- ggplot(dflogCPM, aes(x, y, color=condition)) +
    geom_point(size = 3) +
    scale_color_manual(values = unique(cond_colours), name="") +
    theme_bw() +
    xlab("Leading logFC dim 1") +
    ylab("Leading logFC dim 2") +
    theme(aspect.ratio=1) +
    theme(legend.position="bottom") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    ggtitle("MDS (logCPM)\nwithout sample pairing") +
    ggrepel::geom_text_repel(aes(x, y, label = sample), color="black") +
    theme(plot.title = element_text(face="bold"))

  mds.logcpmc <- ggplot(dflogCPMc, aes(x, y, color=condition)) +
    geom_point(size = 3) +
    scale_color_manual(values = unique(cond_colours), name="") +
    theme_bw() +
    xlab("Leading logFC dim 1") +
    ylab("Leading logFC dim 2") +
    theme(aspect.ratio=1) +
    theme(legend.position="bottom") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    ggtitle("MDS (logCPM)\nwith sample pairing") +
    ggrepel::geom_text_repel(aes(x, y, label = sample), color="black") +
    theme(plot.title = element_text(face="bold"))

  pmds <- plot_grid(mds.logcpm, mds.logcpmc, ncol = 2,align = "hv")

  ggsave("MDS_plot_batchEffect.pdf", pmds, units = "in", width = 7, height = 7, dpi=200)

  ### Plot expression profiles
  pdens <- ggplot(reshape2::melt(logCPM), aes(value, color=Var2)) +
    geom_density() +
    #theme_bw() +
    scale_color_brewer(palette = "Set3") +
    ggtitle("Expression profiles") +
    xlab("Log10 of normalized expression per gene (DESeq2)") +
    ylab("Density") +
    theme(plot.title = element_text(face="bold")) +
    theme(legend.position="bottom") +
    labs(color = "")

  ggsave("normalized_gene_expression_check.pdf", pdens, units = "in", width = 7, height = 7, dpi=200)

  # Plot BCV plot
  pdf(file="BCV_plot.pdf")
  plotBCV(d, col.tagwise = "black")
  dev.off()

  # Plot tagwise mean variation
  pdf(file="edgeR_mean_variation_tgwDisp.pdf")
  mv<-plotMeanVar(d, show.raw.vars=TRUE, show.tagwise.vars=TRUE, NBline=TRUE, main="Mean Variation Tagwise")
  dev.off()

  sink("edgeR_de_genes_check.txt")
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " Without LogFC Cut-off", sep=""))
  dt<-decideTestsDGE(lrt_tgw, adjust.method="BH", p.value=P_THRESHOLD)
  summary(dt)
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " and LogFC >= ", round(LFC_THRESHOLD, 2), sep=""))
  dt<-decideTestsDGE(lrt_tgw, adjust.method="BH", p.value=P_THRESHOLD, lfc=LFC_THRESHOLD)
  summary(dt)
  sink()


  # Store names of DE genes with FDR filtration
  de.genes.tgw<-rownames(resultsTbl.tgw)[resultsTbl.tgw$FDR<=P_THRESHOLD]

  TOP_BCKP<-TOP

  if(TOP>summary(dt)[1,]+summary(dt)[3,]){
    TOP<-summary(dt)[1,]+summary(dt)[3,]
  }

  if(TOP > 1){
    # For labels in the MA plot
    de.tags<-rownames(topTags(lrt_tgw, n=TOP)$table)
    gene.labels<-combResults.tgw[rownames(combResults.tgw) %in% de.tags,]

    pdf(file=paste0("edgeR_MAplot_", condsToCompare[2], "_vs_", condsToCompare[1], ".pdf"))
    plotSmear(lrt_tgw, de.tags=de.genes.tgw, main="FC Plot With Tagwise Dispersion")
    abline(h=c(-LFC_THRESHOLD, LFC_THRESHOLD), col="dodgerblue")
    text(x=gene.labels$logCPM, y=gene.labels$log2FoldChange, labels=gene.labels$gene_name, cex=0.7, pos=1)
    dev.off()
  }


  TOP<-TOP_BCKP

  graphics.off()

  ### EdgeR volcano
  ## need genes in rows, P-values, adj. P-values (FDR), logFC and gene labels
  res.tgwForPlot<-combResults.tgw[c(1:5,7)]

  edge.results <- mutate(res.tgwForPlot, sig=ifelse(res.tgwForPlot$padj<P_THRESHOLD, paste0("padj<", P_THRESHOLD), "Not Sig"))
  p <- ggplot(edge.results, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(col=sig), size=0.5) +
    scale_color_manual(values=c("black", "red"))

  edge.results<-edge.results[order(edge.results$padj, edge.results$pvalue),] # make sure it is correctly ordered

  p <- p + geom_text_repel(data=dplyr::filter(edge.results[1:TOP,], padj<P_THRESHOLD), aes(label=gene_name), size=3)+geom_vline(xintercept = 0)+
    geom_vline(xintercept = c(LFC_THRESHOLD, -LFC_THRESHOLD), linetype = "longdash", colour="blue")+
    ggtitle(paste("Volcanoplot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes (edgeR)", sep=""))+
    theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + theme(plot.title = element_text(face="bold"))

  pdf(file=paste("edgeR_volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggplot2.pdf", sep=""))
  print(p)
  dev.off()

  ###
  res.tgwForPlot$baseMean <- exp(res.tgwForPlot$logCPM)
  edge.ma <- ggmaplot(res.tgwForPlot, main =  paste0("MA plot ", condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes (edgeR)"),
                 fdr = P_THRESHOLD, fc = FOLD_CHANGE, size = 0.4,
                 palette = c("#B31B21", "#1465AC", "darkgray"),
                 genenames = as.vector(res.tgwForPlot$gene_name),
                 legend = "top", top = TOP,
                 font.label = c("bold", 11),
                 font.legend = "bold",
                 font.main = "bold",
                 ggtheme = ggplot2::theme_minimal())+
    theme(plot.title = element_text(hjust = 0.5))

  pdf(file=paste("edgeR_MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggpubr.pdf", sep=""))
  print(edge.ma)
  dev.off()


  # # Filtering of all results based on set filtering
  # filterFinalTable<-function(inputTable, ...){
  #   inputTable.fil<-inputTable[inputTable$UpDown!=0,] # Keep only DE genes
  #   inputTable.fil<-inputTable.fil[inputTable.fil$FDR<=P_THRESHOLD,] # Filter by FDR value
  #   inputTable.fil<-inputTable.fil[abs(inputTable.fil$logFC)>=LFC_THRESHOLD,] # Filter by logFC value
  #   inputTable.fil<-inputTable.fil[order(-inputTable.fil$UpDown,inputTable.fil$FDR,
  #                                        -inputTable.fil$logFC),] # ,-inputTable.fil$logCPM),]
  #   if (nrow(inputTable.fil)==0){inputTable.fil[1,1]<-"No results"}
  #   return(inputTable.fil)
  # }
  #
  # combResults.tgw.fil<-filterFinalTable(combResults.tgw)
  #
  # write.table(x = combResults.tgw.fil, file = "edgeR_onlyDE.tsv", sep="\t", row.names=T, col.names=NA)

  ####################################################################################################
  ####################################################################################################
  # Gene ontology - TODO
  # You can use goseq
  # Visualisation - TODO
  # You can use glimma
  ####################################################################################################
  ####################################################################################################
  # Overlap between DESeq2 and edgeR and Venn diagrams - TODO
  # http://www.ats.ucla.edu/stat/r/faq/venn.htm

  library("limma")

  selectDeseq2<-rownames(res2[(abs(res2$log2FoldChange) >= LFC_THRESHOLD) & (res2$padj < P_THRESHOLD) & !is.na(res2$padj),]) # with LFC cut-off
  selectEdger<-rownames(combResults.tgw[(abs(combResults.tgw$log2FoldChange) >= LFC_THRESHOLD) & (combResults.tgw$padj < P_THRESHOLD) & !is.na(combResults.tgw$padj),]) # with LFC cut-off
  #selectDeseq2<-rownames(res2[(res2$padj < P_THRESHOLD) & !is.na(res2$padj),]) # without LFC cut-off
  #selectEdger<-rownames(combResults.tgw[(combResults.tgw$padj < P_THRESHOLD) & !is.na(combResults.tgw$padj),]) # without LFC cut-off

  vennTable<-matrix(ncol=2, nrow=length(unique(c(selectDeseq2, selectEdger))))
  rownames(vennTable)<-unique(c(selectDeseq2, selectEdger))
  colnames(vennTable)<-c("DESeq2", "edgeR")

  vennTable[,1]<-rownames(vennTable) %in% selectDeseq2
  vennTable[,2]<-rownames(vennTable) %in% selectEdger

  vennTable2<-vennCounts(vennTable)

  pdf(file="overlap_DESeq2_edgeR_venn.pdf")
  par(oma=c(0,0,1,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  vennDiagram(vennTable2, counts.col="black", main=paste("Overlap Between DESeq2 and edgeR results \nwith LogFC >= ",
                                                         round(LFC_THRESHOLD, 2), " and adj.pval ", P_THRESHOLD, sep=""))
  dev.off()

  # Add gene names
  vennTable<-as.data.frame(vennTable)
  vennTable<-merge(vennTable, parsedEnsembl, by="row.names", all.x=T) # Add Ensembl genotype info
  rownames(vennTable)<-vennTable$Row.names
  vennTable$Row.names<-NULL
  vennTable$gene_id<-NULL
  vennTable<-vennTable[order(rowSums(vennTable[,c(1:2)]), vennTable[,1], decreasing = T),]

  write.table(vennTable, file = "overlap_DESeq2_edgeR_de.tsv", sep="\t", col.names=NA)

  graphics.off()


}

# # Save session info, history and image
# sink("session_info.txt")
# sessionInfo()
# sink()
#
# save.image(file = "image.RData")
# savehistory(file = "history.Rhistory")

# develop and test 2
#  STAGE <- "stage408"
#  ANALYSIS <- "Vacek-Soucek_RNA-Seq.swapped_samples_-_paired"
#  mRNA_DE <- "mRNA_DE_RSEM"
#
#  #complete.feature_count.tsv
#  #COMPARISON <- c("PCO_H_vs_PCO_T","PT_H_vs_PT_T","PCO_H_vs_PT_H","PCO_T_vs_PT_T","PCO_T_vs_COM_H","PCO_H_vs_COM_H","PT_T_vs_COM_H","PT_H_vs_COM_H")
#  COMPARISON <- c("Trop_2_BirA_DOX_plus_vs_Trop_2_BirA_DOX_minus")
#
#  for(comp in 1:length(COMPARISON)){
#  args <- character(9)
#  args[1] <- paste0("/mnt/ssd/ssd_1/snakemake/",STAGE,"_",ANALYSIS,"/",mRNA_DE,"/",ANALYSIS,".differential_expresion.config.json")
#  args[2] <- paste0("/mnt/ssd/ssd_1/snakemake/",STAGE,"_",ANALYSIS,"/",mRNA_DE,"/complete.RSEM.RData")
#  args[3] <- paste0("/mnt/ssd/ssd_1/snakemake/",STAGE,"_",ANALYSIS,"/",mRNA_DE,"/",COMPARISON[comp],"/all")
#  args[4] <- "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/GRCh38-p10.sqlite.gz"
#  args[5] <- "/mnt/ssd/ssd_3/references/general/default/annot/biotypes_list_mod.txt"
#  args[6] <- COMPARISON[comp]
#  args[7] <- "RSEM"
#  args[8] <- "homsap"
#  args[9] <- "True"
#  args[10] <- "False"
#
# # run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
# print(COMPARISON[comp])
#  }
#
