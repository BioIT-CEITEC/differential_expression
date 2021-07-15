####################################################################################################
####################################################################################################
#
# Script to calculate differential gene expression using DESeq2 and edgeR package
# Designed for mRNA differential gene expression analysis based on Ensembl results using featureCounts, htseq-count, STAR counts
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
# TODO: Add nice volcano and MA plots also for edgeR results if possible
# TODO: Add tSNE clustering
# TODO: Add shared DE genes visualization (inspiration at VBCF BioComp)
####################################################################################################

#devtools::install_github("r-lib/later")

library(data.table)

run_all <- function(args){

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

  orig_colnames_order <- config_tab$full_name

  #remove unused samples
  condition_design <- config_tab$conditions_to_compare[1]
  if(condition_design != "all"){
    condition_design <- strsplit(condition_design,",")[[1]]
    condition_design <- condition_design[which(grepl(condsToCompare[1],condition_design) & grepl(condsToCompare[2],condition_design))]
    condition_design <- strsplit(condition_design,":")[[1]]
    setkey(config_tab,condition)
    config_tab <- config_tab[condition_design,]
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
    samples_desc <- config_tab[,list(sample = full_name,name = full_name,condition,patient = tag)]
  } else {
    if(length(unique(config_tab$donor)) > 1){
      samples_desc <- config_tab[,list(sample = full_name,name = full_name,condition,patient = donor)]
    } else {
      samples_desc <- config_tab[,list(sample = full_name,name = full_name,condition,patient = info)]
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

  if(counts_type == "featureCounts"){
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

  if(ref_from_trans_assembly != T){
    if(organism == "homsap") {
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
    if(counts_type == "featureCounts"){
      parsedEnsembl <- data.frame(gene_id = rownames(mrcounts),gene_name = rownames(mrcounts),gene_biotype = "protein_coding",row.names = rownames(mrcounts))
    }else{
      parsedEnsembl <- data.frame(gene_id = rownames(txi$counts),gene_name = rownames(txi$counts),gene_biotype = "protein_coding",row.names = rownames(txi$counts))
    }
  }
  

  if(counts_type == "featureCounts"){
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
  cond_colours<-brewer.pal(length(unique(conds)), "Paired")[as.factor(conds)]
  names(cond_colours)<-conds

  pdf(file="counts_barplot.pdf", width=10, height=8)
  par(mar = c(7,5,4,2) + 0.1)
  bp<-barplot(apply(mrcounts, 2, sum), las=2, col=cond_colours, ylim=c(0,(max(apply(mrcounts,2,sum)))*1.1))
  text(bp, apply(mrcounts,2,sum), labels=apply(mrcounts, 2, sum), cex=1, pos=3) # https://stackoverflow.com/questions/27466035/adding-values-to-barplot-of-table-in-r
  legend("topleft", levels((conds)), cex=0.6, fill=cond_colours[levels(conds)])
  #  abline(h = 20000000, col="red")
  dev.off()

  sink("DESeq2_design_control.txt")
  print(coldata)
  sink()

  ####################################################################################################
  ####################################################################################################
  # DESeq2 part
  # Make the count object, normalise, dispersion and testing
  library("DESeq2")

  if(counts_type == "featureCounts"){
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

  pdf(file="DESeq2_disperison_plot.pdf")
  plotDispEsts(cds, main="Dispersion Plot")
  dev.off()

  rawcounts<-counts(cds, normalized=FALSE) # Save raw counts
  normcounts<-counts(cds, normalized=TRUE) # Save normalized counts
  log2counts<-log2(normcounts+1) # Save log2 of normalized counts

  vsd<-DESeq2::varianceStabilizingTransformation(cds) # Save counts tranformed with variance Stabilizing Transformation
  vstcounts<-assay(vsd)

  # Normalization check
  pdf(file="pre_post_norm_counts.pdf")
  par(mfrow=c(3,1))
  #  par(mar = c(10,5,4,2) + 0.1)
  barplot(colSums(rawcounts), col=cond_colours, las=2,cex.names=1.3,main="Pre Normalised Counts", ylim=c(0,(max(apply(rawcounts,2,sum)))*1.1))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend("center", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=0.6, horiz=TRUE)
  barplot(colSums(normcounts), col=cond_colours, las=2, cex.names=1.3, main="Post Normalised Counts", ylim=c(0,(max(apply(normcounts,2,sum)))*1.1))
  dev.off()

  # Heatmaps
  library("gplots")

  pdf(file="heatmaps_samples.pdf")
  heatmap.2(cor(vstcounts), trace="none", col=hmcol, main="Sample to Sample Correlation (VST)", RowSideColors=cond_colours, margins=c(9.5,9.5))
  heatmap.2(cor(log2counts), trace="none", col=hmcol, main="Sample to Sample Correlation (Log2)", RowSideColors=cond_colours, margins=c(9.5,9.5))
  heatmap.2(cor(rawcounts), trace="none", col=hmcol, main="Sample to Sample Correlation (Raw Counts)", RowSideColors=cond_colours, margins=c(9.5,9.5))
  dev.off()

  if(PAIRED==TRUE){ # If paired - remove batch effect
    print("Plotting sample heatmaps with batch correction as well.")

    log2counts_batch<-limma::removeBatchEffect(log2(normcounts+1), cds$patient) # Save log2 of normalized counts
    rawcounts_batch<-limma::removeBatchEffect(counts(cds, normalized=FALSE), cds$patient) # Save raw counts

    vsd_batch <- vsd
    assay(vsd_batch) <- limma::removeBatchEffect(assay(vsd_batch), vsd_batch$patient) # Designed for log transform so we approximate it with vsd (less negative values than with rlog); we might get negative values which we 'round' to
    vstcounts_batch <- assay(vsd_batch)

    pdf(file="heatmaps_samples_batch.pdf")
    heatmap.2(cor(vstcounts_batch), trace="none", col=hmcol, main="Sample to Sample Correlation (VST)", RowSideColors=cond_colours, margins=c(9.5,9.5))
    heatmap.2(cor(log2counts_batch), trace="none", col=hmcol, main="Sample to Sample Correlation (Log2)", RowSideColors=cond_colours, margins=c(9.5,9.5))
    heatmap.2(cor(rawcounts_batch), trace="none", col=hmcol, main="Sample to Sample Correlation (Raw Counts)", RowSideColors=cond_colours, margins=c(9.5,9.5))
    dev.off()
  }

  vstcounts<-vstcounts[apply(vstcounts, 1, max) != 0,]
  pca<-princomp(vstcounts)

  pdf(file="sample_to_sample_PCA.pdf")
  # First two components
  par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
  text(pca$loadings, as.vector(colnames(mrcounts)), pos=3)
  legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  # Three components
  par(mfrow=c(1,3), oma=c(2,0,0,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  plot(pca$loadings[,c(1,2)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC2", xlab="PC1")
  text(pca$loadings[,c(1,2)], as.vector(colnames(mrcounts)), pos=3)
  plot(pca$loadings[,c(1,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC1")
  text(pca$loadings[,c(1,3)], as.vector(colnames(mrcounts)), pos=3)
  plot(pca$loadings[,c(2,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC2")
  text(pca$loadings[,c(2,3)], as.vector(colnames(mrcounts)), pos=3)
  legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  dev.off()

  #if(INTERCEPT == T){
  #pdf(file="sample_to_sample_PCA_batchEffect.pdf")
  #  # First two components
  #  par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  #  plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
  #  text(pca$loadings, as.vector(colnames(mrcounts)), pos=3)
  #  legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  #  # Three components
  #  par(mfrow=c(1,3), oma=c(2,0,0,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  #  plot(pca$loadings[,c(1,2)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC2", xlab="PC1")
  #  text(pca$loadings[,c(1,2)], as.vector(colnames(mrcounts)), pos=3)
  #  plot(pca$loadings[,c(1,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC1")
  #  text(pca$loadings[,c(1,3)], as.vector(colnames(mrcounts)), pos=3)
  #  plot(pca$loadings[,c(2,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC2")
  #  text(pca$loadings[,c(2,3)], as.vector(colnames(mrcounts)), pos=3)
  #  legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  #dev.off()
  #}

  # Get PCA with batch effect from DESeq2 results https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples
  if(PAIRED == T){
    library("ggplot2")
    library("ggrepel")

    pdf(file="sample_to_sample_PCA_batchEffect.pdf")
    vsd_batch <- DESeq2::varianceStabilizingTransformation(dds)
    pcaData <- plotPCA(vsd_batch, intgroup=c("condition", "patient"), returnData=TRUE)
    pcaData$condition <- as.vector(pcaData$condition)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    a<-ggplot(pcaData, aes(PC1, PC2, color=condition, shape=patient)) +
      geom_point(size=3) +
      ggrepel::geom_text_repel(aes(PC1, PC2, label = rownames(pcaData))) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      ggtitle("PCA (DESeq2 VST) without a batch effect removed.")

    assay(vsd_batch) <- limma::removeBatchEffect(assay(vsd_batch), vsd_batch$patient) # Designed for log transform so we approximate it with vsd (less negative values than with rlog); we might get negative values which we 'round' to 0
    pcaData <- plotPCA(vsd_batch, intgroup=c("condition", "patient"), returnData=TRUE)
    pcaData$condition <- as.vector(pcaData$condition)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    b<-ggplot(pcaData, aes(PC1, PC2, color=condition, shape=patient)) +
      geom_point(size=3) +
      ggrepel::geom_text_repel(aes(PC1, PC2, label = rownames(pcaData))) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      ggtitle("PCA (DESeq2 VST) with a batch effect removed.")
    print(a)
    print(b)
    dev.off()

    write.table(x = assay(vsd_batch), file = "norm_counts_batch_corrected.tsv", sep = "\t", col.names = NA)
  }

  pdf(file="contributions_PCA.pdf")
  plot(pca, type = "l", main="Principal Component Contributions")
  dev.off()

  # pca2<-prcomp(t(vstcounts),scale=TRUE,center=TRUE)
  #
  # pdf(file="sample_to_sample_PCA2.pdf")
  #   par(mfrow=c(1,3), oma=c(2,0,0,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  #   plot(pca2$x[,1],pca2$x[,2], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC2",xlab="PC1")
  #   text(pca2$x, as.vector(colnames(mrcounts)), pos=3, cex=1)
  #   plot(pca2$x[,1],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC1")
  #   text(pca2$x[,1],pca2$x[,3], as.vector(colnames(mrcounts)), pos=3, cex=1)
  #   plot(pca2$x[,2],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC2")
  #   text(pca2$x[,2],pca2$x[,3], as.vector(colnames(mrcounts)), pos=3, cex=1)
  #   legend(-75,-38, levels(conds), fill=cond_colours[levels(conds)], cex=1)
  # dev.off()
  #
  # pdf(file="contributions_PCA.pdf")
  #   plot(pca2, type = "l", main="Principal Component Contributions")
  # dev.off()

#  library("plotly")
#
#  tab <- as.data.table(pca$loadings[,1:3])
#  tab[,cond := names(cond_colours)]
#  tab[,cond_colour := cond_colours]
#  tab[,sample := samples_desc$name]
#
#  #c('#BF382A', '#0C4B8E') maybe better colors
#  p <- plot_ly(tab, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3, color = ~cond, colors = unique(tab$cond_colour),text = ~paste('Sample:', sample)) %>% add_markers()
#  htmlwidgets::saveWidget(plotly::as_widget(p),"samples_to_sample_3D_PCA.html")



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
      res<- lfcShrink(dds, coef=grep(paste0("condition_", condsToCompare[2], "_vs"), resultsNames(dds)), res=res)
      # resLFC <- lfcShrink(dds, coef=2, res=res) # Gives different results than contrast=c("condition", condsToCompare[2], condsToCompare[1]
    }else{
      res<- lfcShrink(dds, contrast=c("condition", condsToCompare[2], condsToCompare[1]), res=res) # this removes column lfcSE
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
      resNoFil<- lfcShrink(dds, coef=grep(paste0("condition_", condsToCompare[2], "_vs"), resultsNames(dds)), res=resNoFil)
      # resNoFilLFC <- lfcShrink(dds, coef=2, res=resNoFil) # Gives different results than contrast=c("condition", condsToCompare[2], condsToCompare[1]
    }else{
      resNoFil<- lfcShrink(dds, contrast=c("condition", condsToCompare[2], condsToCompare[1]), res=resNoFil) # this removes column lfcSE
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
  resForPlot$Gene<-parsedEnsembl[rownames(res), "gene_name"]
  results = mutate(resForPlot, sig=ifelse(resForPlot$padj<P_THRESHOLD, paste0("padj<", P_THRESHOLD), "Not Sig"))
  p = ggplot(results, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(col=sig), size=0.5) +
    scale_color_manual(values=c("black", "red"))

  results<-results[order(results$padj, results$pvalue),] # make sure it is correctly ordered

  p = p + geom_text_repel(data=filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)+geom_vline(xintercept = 0)+
    geom_vline(xintercept = c(LFC_THRESHOLD, -LFC_THRESHOLD), linetype = "longdash", colour="blue")+
    ggtitle(paste("Volcanoplot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""))+
    theme(plot.title = element_text(hjust = 0.5))

  pdf(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggplot2.pdf", sep=""))
  print(p)
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
  library("ggpubr")

  # Add custom gene names to MA plot - there is a issue with coloring of the extra genes which are being considered as DE (and colored) even if they are not
  # selectedGenes<-which(resForPlot$Gene %in% c("IL2RA", "IL2RB"))
  # resForPlot[selectedGenes, "padj"] <- 0
  #
  # TOP<-TOP+length(selectedGenes)
  #
  # pdf(file=paste("MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],"_ggpubr.pdf", sep=""))
  # #ggmaplot(resForPlot, main =  expression("Group 1" %->% "Group 2"),
  #   ggmaplot(resForPlot, main =  paste0("MA plot ", condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes"),
  #            fdr = P_THRESHOLD, fc = 0, size = 0.4,
  #            palette = c("#B31B21", "#1465AC", "darkgray"),
  #            genenames = as.vector(resForPlot$Gene),
  #            legend = "top", top = TOP,
  #            font.label = c("bold", 11),
  #            font.legend = "bold",
  #            font.main = "bold",
  #            ggtheme = ggplot2::theme_minimal())+
  #     theme(plot.title = element_text(hjust = 0.5))
  # dev.off()

  pdf(file=paste("MAplot_", condsToCompare[2], "_vs_", condsToCompare[1],"_ggpubr.pdf", sep=""))
  #ggmaplot(resForPlot, main =  expression("Group 1" %->% "Group 2"),
  p <- ggmaplot(resForPlot, main =  paste0("MA plot ", condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes"),
           fdr = P_THRESHOLD, fc = FOLD_CHANGE, size = 0.4,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           genenames = as.vector(resForPlot$Gene),
           legend = "top", top = TOP,
           font.label = c("bold", 11),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())+
    theme(plot.title = element_text(hjust = 0.5))

  print(p)
  dev.off()

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
  resForPlot$Gene<-parsedEnsembl[rownames(resNoFil), "gene_name"]
  results = mutate(resForPlot, sig=ifelse(resForPlot$padj<P_THRESHOLD, paste0("padj<", P_THRESHOLD), "Not Sig"))
  p = ggplot(results, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(col=sig), size=0.5) +
    scale_color_manual(values=c("black", "red"))
  #p
  results<-results[order(results$padj, results$pvalue),] # make sure it is correctly ordered
  #p+geom_text(data=filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)

  pdf(file=paste("volcanoplot_", condsToCompare[2], "_vs_", condsToCompare[1], "_noIndFilt_ggplot2.pdf", sep=""))
  p <- p + geom_text_repel(data=filter(results[1:TOP,], padj<P_THRESHOLD), aes(label=Gene), size=3)+geom_vline(xintercept = 0)+
    geom_vline(xintercept = c(LFC_THRESHOLD, -LFC_THRESHOLD), linetype = "longdash", colour="blue")+
    ggtitle(paste("Volcanoplot ",condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes", sep=""))+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
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
  p <- ggmaplot(resForPlot, main =  paste0("MA plot ", condsToCompare[2], " vs ", condsToCompare[1], " top ", TOP, " genes"),
           fdr = P_THRESHOLD, fc = FOLD_CHANGE, size = 0.4,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           genenames = as.vector(resForPlot$Gene),
           legend = "top", top = TOP,
           font.label = c("bold", 11),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
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
  rownames(log2.norm.counts)<-parsedEnsembl[rownames(log2.norm.counts), "gene_name"]
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
    pdf(file="heatmap_selected_orderBaseMeanCluster.pdf", onefile=FALSE, height= 2 + (TOP/2))
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)", sep=""))
    dev.off()

    pdf(file="heatmap_selected_orderBaseMean.pdf", onefile=FALSE, height= 2 + (TOP/3))
    pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE,
             cluster_cols=FALSE, annotation_col=df, main = paste("Top ", TOP, " significantly DE genes (log2norm)", sep=""))
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

  if(counts_type == "featureCounts"){
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
  combResults.tgw<-cbind(resultsTbl.tgw,"tgw.Disp"=d$tagwise.dispersion[wh.rows.tgw],
                         "UpDown"=decideTestsDGE(lrt_tgw,adjust.method="BH", p.value=P_THRESHOLD,
                                                 lfc=LFC_THRESHOLD)[wh.rows.tgw], cpm(d$counts[wh.rows.tgw,]),d$counts[wh.rows.tgw,])

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

  cond_colours<-brewer.pal(length(unique(d$samples$group)), "Paired")[d$samples$group]
  names(cond_colours)<-d$samples$group

  pdf(file="MDS_plot.pdf") # Plot biological coefficient of variation. The closer samples are together the more they are similar - should be True for replicates
  par(mfrow=c(1,1),xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  plotMDS(d, main="MDSPlot_BCV_distance", xlab="BCV distance 1", ylab="BCV distance 2", labels = rownames(d$samples),
          method="bcv", col=cond_colours, pch=points[as.numeric(coldata$patient)], font=2)
  legend("topleft", legend=levels(d$samples$group), fill=cond_colours[levels(d$samples$group)], ncol=1, cex = 0.75)
  # legend("topright", legend=rownames(d$samples), pch=points[as.numeric(coldata$patient)], ncol=1, cex = 0.75)
  dev.off()

  # See the effect of batch - should be included in the formula
  A<-aveLogCPM(d)
  d2<-d[A>1,]
  d2<-calcNormFactors(d2)
  logCPM<-cpm(d2, log=TRUE, prior.count=5)
  logCPMc<-removeBatchEffect(logCPM, coldata$patient)

  pdf(file="MDS_plot_batchEffect.pdf")
  par(mfrow=c(1,2), oma=c(1,0,0,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  plotMDS(logCPM, col=cond_colours, main="MDS without sample pairing (logCPM)", font=2)
  #	legend("topleft", levels(d$samples$group), fill=cond_colours[levels(d$samples$group)], cex=0.6, horiz = FALSE)
  #	plot(1, type="n", axes=F, xlab="", ylab="")
  #	legend("center", levels(d$samples$group), fill=colors[as.numeric(unique(d$samples$group))], cex=0.6, horiz = FALSE) # This places legend in the middle but set par(mfrow=c(1,3))
  plotMDS(logCPMc, col=cond_colours, main="MDS with sample pairing (logCPM)", font=2)
  legend(-1,-0.45, levels(d$samples$group), fill=cond_colours[levels(d$samples$group)], cex=0.6)
  dev.off()

  ### Plot expression profiles
  color<-rainbow(n = ncol(logCPM))
  density_plot <- density(logCPM[, 1])
  # Get min and max for the plot
  minForPlotX<-min(density_plot$x)
  maxForPlotX<-max(density_plot$x)
  maxForPlotY<-max(density_plot$y)

  # TODO Replace with recode() from dplyr
  for (s in 2:ncol(logCPM)){
    logcounts <- logCPM[,s]+1
    density_plot <- density(logcounts)
    if(min(density_plot$x)<minForPlotX){
      minForPlotX<-min(density_plot$x)
    }
    if(max(density_plot$x)>maxForPlotX){
      maxForPlotX<-max(density_plot$x)
    }
    if(max(density_plot$y)>maxForPlotY){
      maxForPlotY<-max(density_plot$y)
    }
  }

  density_plot <- density(logCPM[, 1])

  pdf(paste0("normalized_gene_expression_check.pdf"), width=10)
  plot(density_plot, main="Expression profiles", xlim=c(minForPlotX*1.3, maxForPlotX*1.3),
       ylim=c(0, maxForPlotY*1.3), xlab=paste0("Log10 of normalized expression per gene (DESeq2)"),
       ylab="Density", col=color[1])

  for (s in 2:ncol(logCPM)){
    logcounts <- logCPM[,s]
    density_plot <- density(logcounts)
    lines(density_plot, col=color[s])
  }

  legend("topright", legend=colnames(logCPM), col = color, cex = .5, fill=color)
  dev.off()

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
# args <- character(9)
# args[1] <- "/mnt/ssd/ssd_1/snakemake/stage283_Termiti_peska.qseq_DE/mRNA_DE_RSEM/Termiti_peska.qseq_DE.differential_expresion.config.json"
# args[2] <- "/mnt/ssd/ssd_1/snakemake/stage283_Termiti_peska.qseq_DE/mRNA_DE_RSEM/complete.RSEM.RData"
# args[3] <- "/mnt/ssd/ssd_1/snakemake/stage283_Termiti_peska.qseq_DE/mRNA_DE_RSEM/Young_k_vs_Worker/all/"
# args[4] <- "/mnt/ssd/ssd_3/references/termite/Peska_transcriptome/annot/Peska_transcriptome.sqlite.gz"
# args[5] <- "/mnt/ssd/ssd_3/references/general/default/annot/biotypes_list_mod.txt"
# args[6] <- "Young_k_vs_Worker"
# args[7] <- "RSEM"
# args[8] <- "termite"
# args[9] <- "False"
# args[10] <- "True"

# # run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)
