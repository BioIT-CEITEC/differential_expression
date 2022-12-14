
hmcol <<- colorRampPalette(brewer.pal(9, "GnBu"))(100)

count_matrix_from_dt <- function(count_dt, value_var = "count", condition_to_compare_vec = condition_to_compare_vec){
  res <- dcast.data.table(count_dt,Ensembl_Id ~ condition + sample_name,value.var = value_var, sep="::")
  names(res) <- gsub(paste(unlist(paste0(condition_to_compare_vec,"::")), collapse = "|"), "", names(res))
  #names(res) <- names(res)
  mat <- as.matrix(res[,-1,with = F],rownames.force = T)
  rownames(mat) <- res$Ensembl_Id
  return(mat)
}


get_title_from_design <- function(experiment_design,prefix = "",suffix = "",non_dual_text = "",connection = " "){
  if(length(unique(experiment_design$condition)) > 2) {
    main_text <- non_dual_text
  } else {
    main_text <- paste0(unique(experiment_design$condition)[2]," vs ",unique(experiment_design$condition)[1])
  }
  if(main_text != ""){
    if(prefix != ""){
      prefix <- paste0(prefix,connection)
    }
    if(suffix != ""){
      suffix <- paste0(connection,suffix)
    }
    return(paste0(prefix,main_text,suffix))
  } else {
    if(suffix != ""){
      suffix <- paste0(connection,suffix)
    }
    return(paste0(prefix,suffix))
  }
}


prepare_colors_and_shapes <- function(experiment_design){

  if(length(unique(levels(experiment_design$condition))) >= 3){
    num.conds <- length(unique(levels(experiment_design$condition)))
  }else{
    num.conds <- 3
  }
  experiment_design[,cond_colours := c(brewer.pal(num.conds, "Dark2"),brewer.pal(num.conds, "Set2"),brewer.pal(num.conds, "Pastel2"))[experiment_design$condition]]
  # cond_shapes<-c("\u25A0","\u25B2","\u25C6","\u25CF","\u25BC","\u25B6","\u25C0","\u25A3","\u25C8","\u25C9","\u25E9","\u25EA")[]
  cond_shapes<-rep(c(16,15,18,3,4,1,0,5,17,2,6,8,11,12),4)
  cond_shapes_plotly <- rep(c("circle","square","diamond","cross","x","circle-open","square-open","diamond-open"),4)
  experiment_design[,pat_shapes := cond_shapes[experiment_design$patient]]
  experiment_design[,pat_shapes_plotly := cond_shapes_plotly[experiment_design$patient]]


  return(experiment_design)
}


DESeq2_computation <- function(txi = NULL,count_dt = NULL,experiment_design,remove_genes_with_mean_read_count_threshold,condition_to_compare_vec){

  paired_samples <- length(unique(experiment_design$patient)) != nrow(experiment_design)

  if(!is.null(txi)){
    analysis_type <- "RSEM"
  } else {
    analysis_type <- "feature_count"
  }

  if(analysis_type == "feature_count"){
    #FeatureCount
    if(paired_samples==T){
      print("Using paired design in DESeq2")
      dds<-DESeqDataSetFromMatrix(countData = count_matrix_from_dt(count_dt,"count",condition_to_compare_vec = condition_to_compare_vec), colData = experiment_design, design = ~patient+condition)
    }else{
      print("Using simple design in DESeq2")
      dds<-DESeqDataSetFromMatrix(countData = count_matrix_from_dt(count_dt,"count",condition_to_compare_vec = condition_to_compare_vec), colData = experiment_design, design = ~condition)
    }
  } else {
    #RSEM
    if(paired_samples==T){
      print("Using paired design in DESeq2")
      dds<-DESeqDataSetFromTximport(txi = txi, colData = experiment_design, design = ~patient+condition)
    }else{
      print("Using simple design in DESeq2")
      dds<-DESeqDataSetFromTximport(txi = txi, colData = experiment_design, design = ~condition)
    }
  }

  # Remove very low count genes
  keep <- rowSums(counts(dds)) >= remove_genes_with_mean_read_count_threshold
  dds <- dds[keep,]

  # The same thing which follows be calculated by >DESeq(dds) instead of three separate commands
  dds<-estimateSizeFactors(dds)
  dds<-estimateDispersions(dds)
  dds<-nbinomWaldTest(dds)

  count_dt[,rawcounts := as.vector(t(counts(dds, normalized=FALSE)))]
  count_dt[,normcounts := as.vector(t(counts(dds, normalized=TRUE)))]
  count_dt[,log2counts := log2(normcounts+1)]
  count_dt[,vstcounts := as.vector(t(assay(DESeq2::varianceStabilizingTransformation(dds))))]

  if(paired_samples==TRUE){ # If paired - remove batch effect
    count_dt[,vstcounts_batch := as.vector(t(limma::removeBatchEffect(count_matrix_from_dt(count_dt,"vstcounts",condition_to_compare_vec = condition_to_compare_vec), dds$patient)))]
    count_dt[,rawcounts_batch := as.vector(t(limma::removeBatchEffect(count_matrix_from_dt(count_dt,"rawcounts",condition_to_compare_vec = condition_to_compare_vec), dds$patient)))]
    count_dt[,log2counts_batch := as.vector(t(limma::removeBatchEffect(count_matrix_from_dt(count_dt,"log2counts",condition_to_compare_vec = condition_to_compare_vec), dds$patient)))]
  }

  #count_dt <- merge(experiment_design[,.(sample_name,condition,patient)],count_dt,by = "sample_name")

  #setkey(count_dt,Feature_name,condition,patient,sample_name)

  return(list(dds,count_dt))
}

create_normalization_specific_DESeq2_results <- function(output_dir,dds,count_dt,condition_to_compare_vec){

  #set output dir but remember where to return
  orig_dir <- getwd()
  dir.create(paste0(output_dir,"/report_data"),showWarnings = F,recursive = T)
  setwd(output_dir)

  # get experiment_design_dt from count_dt
  experiment_design_dt <- unique(count_dt[,.(sample_name,condition,patient)])
  setorder(experiment_design_dt,condition,patient)
  experiment_design_dt <- prepare_colors_and_shapes(experiment_design_dt)

  #set paired samples if more then one sample belongs to one patient(batch)
  paired_samples <- length(unique(experiment_design_dt$patient)) != nrow(experiment_design_dt)


  if(paired_samples){
    fwrite(experiment_design_dt[,.(`Sample name` = sample_name,`Condition` = condition,`Batch/patient pairing` = patient)],"DESeq2_experiment_design.tsv",sep = "\t")
  } else {
    fwrite(experiment_design_dt[,.(`Sample name` = sample_name,`Condition` = condition)],"DESeq2_experiment_design.tsv",sep = "\t")
  }



  ####################################################################################################

  pdf(file="DESeq2_disperison_plot.pdf", width=7, height=7)
  plotDispEsts(dds, main="Dispersion Plot")
  dev.off()


  ####################################################################################################

  # Normalization check
  rcs<-count_dt[,.(value = sum(rawcounts), type = "Raw counts"),by = .(sample_name,condition)]
  ncs<-count_dt[,.(value = sum(normcounts), type = "Normalised counts"),by = .(sample_name,condition)]
  bcs<-rbind(rcs,ncs)
  fwrite(bcs,"pre_post_norm_counts.tsv", sep="\t")

  # prepare the title
  count.title <- get_title_from_design(experiment_design_dt,"","","All samples")

  rcs_ncs <- ggplot(bcs, aes(sample_name, value, fill=condition))+
    geom_bar(stat = "identity", width = 0.8)+
    scale_fill_manual(values = unique(experiment_design_dt$cond_colours))+
    theme_bw() +
    ylab("") +
    xlab("") +
    facet_grid(factor(type, levels=c("Raw counts","Normalised counts")) ~ .) +
    ggtitle(paste0("Pre & Post Normalised Counts\n",count.title)) +
    theme(axis.text.x = element_text(angle = 90))+
    theme(plot.title = element_text(face="bold"))

  ggsave(filename = "report_data/pre_post_norm_counts.png", rcs_ncs, units = "in", dpi = 200, width = 7, height = 7, device = "png")
  ggsave(filename = "report_data/pre_post_norm_counts.svg", rcs_ncs, width = 7, height = 7, device = "svg")
  ggsave(filename = "pre_post_norm_counts.pdf", rcs_ncs, width = 7, height = 7, device = "pdf")


  #?? organise by

  ####################################################################################################
  # Heatmaps

  hm.log.title <- get_title_from_design(experiment_design_dt,"Sample to Sample Correlation (Log2)",connection = "\n")
  hm.vst.title <- get_title_from_design(experiment_design_dt,"Sample to Sample Correlation (VST)",connection = "\n")
  hm.raw.title <- get_title_from_design(experiment_design_dt,"Sample to Sample Correlation (Raw Counts)",connection = "\n")


  pdf(file="heatmaps_samples.pdf")
  heatmap.2(cor(count_matrix_from_dt(count_dt,"log2counts",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.log.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
  heatmap.2(cor(count_matrix_from_dt(count_dt,"vstcounts",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
  heatmap.2(cor(count_matrix_from_dt(count_dt,"rawcounts",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.raw.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
  dev.off()

  png(file="report_data/heatmaps_samples_log2.png", width = 7, height = 7, unit = "in", res=200)
  heatmap.2(cor(count_matrix_from_dt(count_dt,"log2counts",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.log.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
  dev.off()
  png(file="report_data/heatmaps_samples_vst.png", width = 7, height = 7, unit = "in", res=200)
  heatmap.2(cor(count_matrix_from_dt(count_dt,"vstcounts",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
  dev.off()

  svg(file="report_data/heatmaps_samples_log2.svg", width = 7, height = 7)
  heatmap.2(cor(count_matrix_from_dt(count_dt,"log2counts",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.log.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
  dev.off()
  svg(file="report_data/heatmaps_samples_vst.svg", width = 7, height = 7)
  heatmap.2(cor(count_matrix_from_dt(count_dt,"vstcounts",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
  dev.off()

  ####################################################################################################

  if(paired_samples==TRUE){ # If paired - remove batch effect
    print("Plotting sample heatmaps with batch correction as well.")

    hm.log.title <- get_title_from_design(experiment_design_dt,"Sample to Sample Correlation (Log2)","with a batch effect removed",connection = "\n")
    hm.vst.title <- get_title_from_design(experiment_design_dt,"Sample to Sample Correlation (VST)","with a batch effect removed",connection = "\n")
    hm.raw.title <- get_title_from_design(experiment_design_dt,"Sample to Sample Correlation (Raw Counts)","with a batch effect removed",connection = "\n")

    pdf(file="heatmaps_samples_batch.pdf")
    heatmap.2(cor(count_matrix_from_dt(count_dt,"log2counts_batch",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.log.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
    heatmap.2(cor(count_matrix_from_dt(count_dt,"vstcounts_batch",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
    heatmap.2(cor(count_matrix_from_dt(count_dt,"rawcounts_batch",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.raw.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
    dev.off()

    png(file="report_data/heatmaps_samples_log_batch.png", width = 7, height = 7, unit = "in", res=200)
    heatmap.2(cor(count_matrix_from_dt(count_dt,"log2counts_batch",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.log.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
    dev.off()
    png(file="report_data/heatmaps_samples_vst_batch.png", width = 7, height = 7, unit = "in", res=200)
    heatmap.2(cor(count_matrix_from_dt(count_dt,"vstcounts_batch",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
    dev.off()

    svg(file="report_data/heatmaps_samples_log_batch.svg", width = 7, height = 7)
    heatmap.2(cor(count_matrix_from_dt(count_dt,"log2counts_batch",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.log.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
    dev.off()
    svg(file="report_data/heatmaps_samples_vst_batch.svg", width = 7, height = 7)
    heatmap.2(cor(count_matrix_from_dt(count_dt,"vstcounts_batch",condition_to_compare_vec = condition_to_compare_vec)), trace="none", col=hmcol, main=hm.vst.title, RowSideColors=experiment_design_dt$cond_colours, margins=c(9.5,9.5))
    dev.off()

  }

  ####################################################################################################

  # count_dt[,max_vstcounts := max(vstcounts),by = Feature_name]
  # pca<-princomp(count_matrix_from_dt(count_dt[max_vstcounts > 0],"vstcounts",condition_to_compare_vec = condition_to_compare_vec))
  # count_dt[,max_vstcounts := NULL]

  # pdf(file="sample_to_sample_PCA.pdf")
  # # First two components
  # par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  # plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
  # text(pca$loadings, as.vector(colnames(count_dt)), pos=3)
  # legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  # # Three components
  # par(mfrow=c(1,3), oma=c(2,0,0,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  # plot(pca$loadings[,c(1,2)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC2", xlab="PC1")
  # text(pca$loadings[,c(1,2)], as.vector(colnames(count_dt)), pos=3)
  # plot(pca$loadings[,c(1,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC1")
  # text(pca$loadings[,c(1,3)], as.vector(colnames(count_dt)), pos=3)
  # plot(pca$loadings[,c(2,3)], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)", ylab="PC3", xlab="PC2")
  # text(pca$loadings[,c(2,3)], as.vector(colnames(count_dt)), pos=3)
  # legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
  # dev.off()

  ####################################################################################################




  get_pca_plot <- function(prcomp_data,experiment_design_dt,pca_title,comp_x_id = 1,comp_y_id = 2){

    pcaData <- as.data.table(prcomp_data$x,keep.rownames = T)
    setnames(pcaData,"rn","sample_name")
    pcaData <- merge(pcaData, experiment_design_dt, by="sample_name")

    PC1 <- paste0("PC",comp_x_id)
    PC2 <- paste0("PC",comp_y_id)

    pca1 <- ggplot(pcaData, aes_string(PC1, PC2, color="condition"))
    if(paired_samples){
      pca1 <- pca1 + geom_point(size=3, aes(shape = patient))
    }else{
      pca1 <- pca1 + geom_point(size=3)
    }
    pca1 <- pca1 + scale_color_manual(values = unique(experiment_design_dt$cond_colours)) +
      theme_bw() +
      scale_shape_manual(values=experiment_design_dt$pat_shapes) +
      xlab(paste0(PC1," - ",round(summary(prcomp_data)$importance[2,comp_x_id] * 100,1),"% variance explained")) +
      ylab(paste0(PC2," - ",round(summary(prcomp_data)$importance[2,comp_y_id] * 100,1),"% variance explained")) +
      theme(plot.title = element_text(face="bold")) +
      theme(legend.position="bottom") +
      ggtitle(pca_title)


    # pca1P = pca1 + ggrepel::geom_text_repel(aes(PC1, PC2, label = patient), color="black", max.overlaps = length(pcaData$sample_name))
    return(pca1 + ggrepel::geom_text_repel(aes_string(PC1, PC2, label = "sample_name"), color="black", max.overlaps = length(pcaData$sample_name)))
  }

  get_prcomp <- function(count_dt,experiment_design_dt,paired=FALSE){
    if(!paired){
      count_dt[,max_vstcounts := max(vstcounts),by = Feature_name]
      vstcounts <- count_matrix_from_dt(count_dt[max_vstcounts > 0],"vstcounts",condition_to_compare_vec = condition_to_compare_vec)
      count_dt[,max_vstcounts := NULL]
    }else{
      count_dt[,max_vstcounts_batch := max(vstcounts_batch),by = Feature_name]
      vstcounts <- count_matrix_from_dt(count_dt[max_vstcounts_batch > 0],"vstcounts_batch",condition_to_compare_vec = condition_to_compare_vec)
      count_dt[,max_vstcounts_batch := NULL]
    }

    # calculate the variance for each gene
    rv <- rowVars(vstcounts)
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
    # perform a PCA on the data in vstcounts for the selected genes
    prcomp_data <- prcomp(t((vstcounts)[select,]))

    return(prcomp_data)
  }

  prcomp_data <- get_prcomp(count_dt,experiment_design_dt,paired=FALSE)
  pca_title <- get_title_from_design(experiment_design_dt,"PCA (DESeq2 VST)")

  pca.dt <- as.data.table(prcomp_data$x, keep.rownames = T)
  setnames(pca.dt,"rn","sample_name")
  fwrite(pca.dt,"sample_to_sample_PCA.tsv", sep="\t")

  pca1S <- get_pca_plot(prcomp_data,experiment_design_dt,pca_title)

  ggsave(filename = "report_data/sample_to_sample_PCA.png", pca1S, units = "in", dpi=200, width = 7, height = 7, device="png")
  # ggsave(filename = "sample_to_sample_PCA2.png", pca1P, units = "in", dpi=200, width = 7, height = 7, device="png")
  ggsave(filename = "report_data/sample_to_sample_PCA.svg", pca1S, width = 7, height = 7, device="svg")
  # ggsave(filename = "sample_to_sample_PCA2.svg", pca1P, width = 7, height = 7, device="svg")
  ggsave(filename = "sample_to_sample_PCA.pdf", pca1S, width = 7, height = 7, device="pdf")

  ####################################################################################################
  # 3D PCA plot with plotly
  ####################################################################################################

  pcaData <- as.data.table(prcomp_data$x,keep.rownames = T)
  setnames(pcaData,"rn","sample_name")
  pcaData <- merge(pcaData, experiment_design_dt, by="sample_name")
  pcaData[,condition := as.character(condition)]
  setorder(pcaData,condition)

  if(paired_samples){
    fig <- plot_ly(pcaData, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, colors = unique(pcaData$cond_colours),hoverinfo = "text",
                   marker = list(symbol = ~pat_shapes_plotly),text = ~paste('Sample name:',sample_name,'<br>Condition:',condition, '<br>Batch pair',patient))
  } else {
    fig <- plot_ly(pcaData, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, colors = unique(pcaData$cond_colours),hoverinfo = "text",
                   text = ~paste('Sample name:',sample_name,'<br>Condition:',condition))
  }
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(title = pca_title,
                        scene = list(xaxis = list(title = paste0("PC1 - ",round(summary(prcomp_data)$importance[2,1] * 100,1),"% variance explained")),
                                     yaxis = list(title = paste0("PC2 - ",round(summary(prcomp_data)$importance[2,2] * 100,1),"% variance explained")),
                                     zaxis = list(title = paste0("PC3 - ",round(summary(prcomp_data)$importance[2,3] * 100,1),"% variance explained"))))
  saveWidget(as_widget(fig), "sample_to_sample_PCA_3D.html")

  ####################################################################################################

  # Get PCA with batch effect from DESeq2 results https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples
  if(paired_samples){

    prcomp_data <- get_prcomp(count_dt,experiment_design_dt,paired=TRUE)

    pca.dt <- as.data.table(prcomp_data$x, keep.rownames = T)
    setnames(pca.dt,"rn","sample_name")
    fwrite(pca.dt,"sample_to_sample_PCA_batch.tsv", sep="\t")

    pca_title <- get_title_from_design(experiment_design_dt,"PCA (DESeq2 VST)","with a batch effect removed")

    pca1_batchS <- get_pca_plot(prcomp_data,experiment_design_dt,pca_title)

    ggsave(filename = "report_data/sample_to_sample_PCA_batch.png", pca1_batchS, units = "in", dpi=200, width = 7, height = 7, device="png")
    # ggsave(filename = "sample_to_sample_PCA_batch2.png", pca1_batchP, units = "in", dpi=200, width = 7, height = 7, device="png")

    ggsave(filename = "report_data/sample_to_sample_PCA_batch.svg", pca1_batchS, width = 7, height = 7, device="svg")
    # ggsave(filename = "sample_to_sample_PCA_batch2.svg", pca1_batchP, width = 7, height = 7, device="svg")

    ggsave(filename = "sample_to_sample_PCA_batch.pdf", pca1_batchS, width = 7, height = 7, device="pdf")


    ####################################################################################################
    # 3D PCA plot with plotly
    ####################################################################################################

    pcaData <- as.data.table(prcomp_data$x,keep.rownames = T)
    setnames(pcaData,"rn","sample_name")
    pcaData <- merge(pcaData, experiment_design_dt, by="sample_name")
    pcaData[,condition := as.character(condition)]
    setorder(pcaData,condition)


    fig <- plot_ly(pcaData, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, colors = unique(pcaData$cond_colours),hoverinfo = "text",
                   marker = list(symbol = ~pat_shapes_plotly),text = ~paste('Sample name:',sample_name,'<br>Condition:',condition, '<br>Batch pair',patient))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(title = pca_title,
                          scene = list(xaxis = list(title = paste0("PC1 - ",round(summary(prcomp_data)$importance[2,1] * 100,1),"% variance explained")),
                                       yaxis = list(title = paste0("PC2 - ",round(summary(prcomp_data)$importance[2,2] * 100,1),"% variance explained")),
                                       zaxis = list(title = paste0("PC3 - ",round(summary(prcomp_data)$importance[2,3] * 100,1),"% variance explained"))))
    saveWidget(as_widget(fig), "sample_to_sample_PCA_3D_batch.html")

  }

  ####################################################################################################

  # if(paired_samples){
  #   fwrite(as.data.table(count_matrix_from_dt(count_dt,"vstcounts_batch")), file = "norm_counts_batch.tsv", sep = "\t", row.names = T)
  # }

  setwd(orig_dir)

  #pdf(file="contributions_PCA.pdf", width = 7, height = 7)
  #plot(pca, type = "l", main="Principal Component Contributions")
  #dev.off()
}

get_comparison_specific_DESeq2_table <- function(dds,count_dt,experiment_design,condsToCompare,output_dir,p_value_threshold,lfc_threshold){
  #set output dir but remember where to return
  orig_dir <- getwd()
  dir.create(paste0(output_dir,"/detail_results"),showWarnings = F,recursive = T)
  setwd(output_dir)

  #create copy of dds not to change orig object
  dds <- copy(dds)

  deseq_obj_comp_res <- results(dds, contrast=c("condition", condsToCompare[1], condsToCompare[2]), independentFiltering=T)
  deseq_obj_comp_res <- lfcShrink(dds, contrast=c("condition", condsToCompare[1], condsToCompare[2]), res=deseq_obj_comp_res, type="normal", lfcThreshold=0)
  comp_res <- as.data.table(deseq_obj_comp_res,keep.rownames=T)
  deseq_obj_comp_res_no_filt <- results(dds, contrast=c("condition", condsToCompare[1], condsToCompare[2]), independentFiltering=F, cooksCutoff=F)
  deseq_obj_comp_res_no_filt <- lfcShrink(dds, contrast=c("condition", condsToCompare[1], condsToCompare[2]), type="normal", lfcThreshold=0,res = deseq_obj_comp_res_no_filt)
  comp_res <- merge(comp_res,as.data.table(deseq_obj_comp_res_no_filt,keep.rownames=T)[,.(rn,no_filter_log2FoldChange = log2FoldChange,no_filter_pvalue = pvalue,no_filter_padj = padj)],by = "rn")
  comp_res[,abs_log2FoldChange := abs(log2FoldChange)]
  comp_res[,abs_no_filter_log2FoldChange := abs(no_filter_log2FoldChange)]
  setnames(comp_res,"rn","Ensembl_Id")
  comp_res <- merge(unique(count_dt[,.(Ensembl_Id,Feature_name,biotype)]),comp_res,by = "Ensembl_Id")
  comp_res[,significant_DE := F]
  comp_res[(abs_log2FoldChange >= lfc_threshold) & (padj < p_value_threshold) & !is.na(padj),significant_DE := T]
  comp_res[,no_filter_significant_DE := F]
  comp_res[(abs_no_filter_log2FoldChange >= lfc_threshold) & (no_filter_padj < p_value_threshold) & !is.na(no_filter_padj),no_filter_significant_DE := T]

  normcounts <- dcast.data.table(count_dt[condition %in% condsToCompare],formula = Ensembl_Id ~ condition + sample_name,value.var = "normcounts")
  names(normcounts) <- gsub(paste(unlist(paste0(condsToCompare,"_")), collapse = "|"), "", names(normcounts))
  setnames(normcounts,names(normcounts)[-1],paste0(names(normcounts)[-1],"_normCounts"))
  comp_res <- merge.data.table(comp_res,normcounts,by = "Ensembl_Id")
  rawcounts <- dcast.data.table(count_dt[condition %in% condsToCompare],formula = Ensembl_Id ~ condition + sample_name,value.var = "rawcounts")
  names(rawcounts) <- gsub(paste(unlist(paste0(condsToCompare,"_")), collapse = "|"), "", names(rawcounts))
  setnames(rawcounts,names(rawcounts)[-1],paste0(names(rawcounts)[-1],"_rawCounts"))
  comp_res <- merge.data.table(comp_res,rawcounts,by = "Ensembl_Id")

  comp_res_summary <- comp_res[, .(test = c("DESeq2","DESeq2_no_filter"),
               total = c(length(Feature_name),length(Feature_name)),
               na = c(sum(is.na(padj)),
                      sum(is.na(no_filter_padj))),
               not_sig = c(sum(!is.na(padj) & padj >= p_value_threshold),
                           sum(!is.na(no_filter_padj) & no_filter_padj >= p_value_threshold)),
               sig = c(sum(!is.na(padj) & padj < p_value_threshold),
                        sum(!is.na(no_filter_padj) & no_filter_padj < p_value_threshold)),
               sig_up = c(sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange > 0),
                           sum(!is.na(no_filter_padj) & no_filter_padj < p_value_threshold & no_filter_log2FoldChange > 0)),
               sig_down = c(sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange < 0),
                             sum(!is.na(no_filter_padj) & no_filter_padj < p_value_threshold & no_filter_log2FoldChange < 0)),
               sig_lfc = c(sum(!is.na(padj) & padj < p_value_threshold & abs_log2FoldChange >= lfc_threshold),
                            sum(!is.na(no_filter_padj) & no_filter_padj < p_value_threshold & abs_no_filter_log2FoldChange >= lfc_threshold)),
               sig_lfc_up = c(sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange >= lfc_threshold),
                               sum(!is.na(no_filter_padj) & no_filter_padj < p_value_threshold & no_filter_log2FoldChange >= lfc_threshold)),
               sig_lfc_down = c(sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange <= -lfc_threshold),
                                sum(!is.na(no_filter_padj) & no_filter_padj < p_value_threshold & no_filter_log2FoldChange <= -lfc_threshold)))]

  fwrite(comp_res_summary, "DESeq2_de_genes_summary.tsv", sep="\t")
  write(paste0("# ",mcols(deseq_obj_comp_res)$description), "DESeq2_de_genes_summary.tsv", ncolumns = 1, append = T)

  # comp_res[,setdiff(names(comp_res),c("no_filter_log2FoldChange","no_filter_pvalue","no_filter_padj","abs_log2FoldChange")),with = F]
  setcolorder(comp_res,c("Ensembl_Id","baseMean","log2FoldChange","lfcSE","pvalue","padj","significant_DE","Feature_name","biotype"))
  setorder(comp_res,padj,pvalue,-abs_log2FoldChange,na.last = T)
  fwrite(comp_res[,setdiff(names(comp_res),c("no_filter_log2FoldChange","no_filter_pvalue","no_filter_padj","no_filter_significant_DE","abs_log2FoldChange","abs_no_filter_log2FoldChange")),with = F], file = "DESeq2.tsv", sep = "\t")
  fwrite(comp_res, file = "detail_results/full_DESeq2.tsv", sep = "\t")

  setwd(orig_dir)

  return(comp_res)
}


create_comparison_specific_DESeq2_results <- function(comp_res,dds,count_dt,condsToCompare,output_dir,paired_samples,TOP,p_value_threshold,lfc_threshold,condition_to_compare_vec){


  ##DESeq2 comparison specific graphic ploting
  ####################################################################################################

  #set output dir but remember where to return
  orig_dir <- getwd()
  dir.create(paste0(output_dir,"/report_data"),showWarnings = F,recursive = T)
  setwd(output_dir)

  TOP_BCKP<-TOP
  TOP <- min(TOP,sum(comp_res$significant_DE))
  RANGE <- seq_len(TOP)

  if(TOP==0){ # Set range for naming the samples - help to avoid naming of samples even if there is no DE gene; THIS ERROR IS OK WHEN WE HAVE 0 DE GENES Error in text.default(res[RANGE, ]$log2FoldChange, -log(res[RANGE, ]$padj,  : zero-length 'labels' specified
    system("touch ZERO_DE_GENES_FOUND")
  }

  if(TOP==1){ # Set range for naming the samples - help to avoid naming of samples even if there is no DE gene; THIS ERROR IS OK WHEN WE HAVE 0 DE GENES Error in text.default(res[RANGE, ]$log2FoldChange, -log(res[RANGE, ]$padj,  : zero-length 'labels' specified
    system("touch ONLY_ONE_DE_GENE_FOUND")
  }

  setorder(comp_res,padj,-abs_log2FoldChange,na.last = T) # Make sure res is ordered by adj.p-value

  if(TOP > 0){
    options(ggrepel.max.overlaps = 2*TOP)
    #volcano plot

    comp_res[,sig := ifelse(padj < p_value_threshold, paste0("padj<", p_value_threshold), "Not Sig")]
    p <- ggplot(comp_res[!is.na(padj)], aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col=sig), size=0.5) +
      scale_color_manual(values=c("black", "red")) +
      geom_text_repel(data=comp_res[RANGE,], aes(label=Feature_name), size=3, max.overlaps = 2*TOP)+
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = c(lfc_threshold, -lfc_threshold), linetype = "longdash", colour="blue") +
      ggtitle(paste("Volcanoplot ",condsToCompare[1], " vs ", condsToCompare[2], " top ", TOP, " genes", sep="")) +
      annotate("text",x=max(comp_res[!is.na(padj)]$log2FoldChange),y=0,vjust=-0.5,label=condsToCompare[1],size=5,fontface = "bold") +
      annotate("text",x=min(comp_res[!is.na(padj)]$log2FoldChange),y=0,vjust=-0.5,label=condsToCompare[2],size=5,fontface = "bold") +
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + theme(plot.title = element_text(face="bold"))

    pdf(file=paste("volcanoplot_", condsToCompare[1], "_vs_", condsToCompare[2],".pdf", sep=""))
    print(p)
    dev.off()

    png(file=paste("report_data/volcanoplot_", condsToCompare[1], "_vs_", condsToCompare[2],".png", sep=""), units = "in",width = 7, height = 7, res = 200)
    print(p)
    dev.off()

    svg(file=paste("report_data/volcanoplot_", condsToCompare[1], "_vs_", condsToCompare[2],".svg", sep=""),width = 7, height = 7)
    print(p)
    dev.off()

    #ma plot

    ma <- ggmaplot(comp_res, main =  paste0("MA plot ", condsToCompare[1], " vs ", condsToCompare[2], " top ", TOP, " genes"),
                   fdr = p_value_threshold, fc = 2^lfc_threshold, size = 0.5,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(comp_res$Feature_name),
                   legend = "top", top = TOP,
                   font.label = c("bold", 11),
                   font.legend = "bold",
                   font.main = "bold",
                   ggtheme = ggplot2::theme_minimal())+
      theme(plot.title = element_text(hjust = 0.5))

    pdf(file=paste("MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],".pdf", sep=""))
    print(ma)
    dev.off()

    png(file=paste("report_data/MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],".png", sep=""), units = "in",width = 7, height = 7, res = 200)
    print(ma)
    dev.off()

    svg(file=paste("report_data/MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],".svg", sep=""),width = 7, height = 7)
    print(ma)
    dev.off()


    #volcano and ma plot for not filtered results

    TOP<-TOP_BCKP
    TOP <- min(TOP,sum(comp_res$no_filter_significant_DE))
    RANGE <- seq_len(TOP)
    options(ggrepel.max.overlaps = 2*TOP)

    #volcano plot

    p <- ggplot(comp_res[!is.na(padj)], aes(no_filter_log2FoldChange, -log10(no_filter_padj))) +
      geom_point(aes(col=sig), size=0.5) +
      scale_color_manual(values=c("black", "red")) +
      geom_text_repel(data=comp_res[RANGE,], aes(label=Feature_name), size=3, max.overlaps = 2*TOP)+
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = c(lfc_threshold, -lfc_threshold), linetype = "longdash", colour="blue") +
      ggtitle(paste("Volcanoplot ",condsToCompare[1], " vs ", condsToCompare[2], " top ", TOP, " genes", sep="")) +
      annotate("text",x=min(comp_res[!is.na(padj)]$log2FoldChange),y=0,vjust=-0.5,label=condsToCompare[1],size=5,fontface = "bold") +
      annotate("text",x=max(comp_res[!is.na(padj)]$log2FoldChange),y=0,vjust=-0.5,label=condsToCompare[2],size=5,fontface = "bold") +
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + theme(plot.title = element_text(face="bold"))

    pdf(file=paste("detail_results/volcanoplot_", condsToCompare[1], "_vs_", condsToCompare[2],"_noIndFilt.pdf", sep=""))
    print(p)
    dev.off()

    png(file=paste("detail_results/report_data/volcanoplot_", condsToCompare[1], "_vs_", condsToCompare[2],"_noIndFilt.png", sep=""), units = "in",width = 7, height = 7, res = 200)
    print(p)
    dev.off()

    svg(file=paste("detail_results/report_data/volcanoplot_", condsToCompare[1], "_vs_", condsToCompare[2],"_noIndFilt.svg", sep=""),width = 7, height = 7)
    print(p)
    dev.off()

    #ma plot

    ma <- ggmaplot(comp_res[,.(Feature_name,baseMean,padj = no_filter_padj,log2FoldChange = no_filter_log2FoldChange)], main =  paste0("MA plot ", condsToCompare[1], " vs ", condsToCompare[2], " top ", TOP, " genes"),
                   fdr = p_value_threshold, fc = 2^lfc_threshold, size = 0.5,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(comp_res$Feature_name),
                   legend = "top", top = TOP,
                   font.label = c("bold", 11),
                   font.legend = "bold",
                   font.main = "bold",
                   ggtheme = ggplot2::theme_minimal())+
      theme(plot.title = element_text(hjust = 0.5))

    pdf(file=paste("detail_results/MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],"_noIndFilt.pdf", sep=""))
    print(ma)
    dev.off()

    png(file=paste("detail_results/report_data/MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],"_noIndFilt.png", sep=""), units = "in",width = 7, height = 7, res = 200)
    print(ma)
    dev.off()

    svg(file=paste("detail_results/report_data/MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],"_noIndFilt.svg", sep=""),width = 7, height = 7)
    print(ma)
    dev.off()

  }

  # Heatmaps of selected genes
  if(TOP >  1){

    if(paired_samples==T){
      print("Using paired design")
      log2.norm.counts <- count_matrix_from_dt(count_dt[condition %in% condsToCompare & Feature_name %in% comp_res[RANGE,]$Feature_name,],"log2counts_batch",condition_to_compare_vec = condition_to_compare_vec)
      df<-unique(count_dt[condition %in% condsToCompare,.(sample_name,condition,patient)])
    }else{
      print("Using simple design")
      log2.norm.counts <- count_matrix_from_dt(count_dt[condition %in% condsToCompare & Feature_name %in% comp_res[RANGE,]$Feature_name,],"log2counts",condition_to_compare_vec = condition_to_compare_vec)
      df<-unique(count_dt[condition %in% condsToCompare,.(sample_name,condition)])
    }
    df <- as.data.frame(df)
    rownames(df)<-df$sample_name
    df$sample_name <- NULL
    log2.norm.counts<-log2.norm.counts[rev(order(rowMeans(log2.norm.counts))),]
    
    
    #pdf(file="heatmap_selected_orderBaseMeanCluster.pdf", onefile=FALSE, height= 2 + (TOP/2))
    pdf(file="heatmap_selected_orderBaseMeanCluster.pdf", onefile=FALSE, height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[1], " vs ", condsToCompare[2], sep=""))
    dev.off()
    
    png(filename = "report_data/heatmap_selected_orderBaseMeanCluster.png", res = 200, units = "in", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[1], " vs ", condsToCompare[2], sep=""))
    dev.off()
    
    svg(filename = "report_data/heatmap_selected_orderBaseMeanCluster.svg", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[1], " vs ", condsToCompare[2], sep=""))
    dev.off()
    
    pdf(file="heatmap_selected_orderBaseMean.pdf", onefile=FALSE, height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[1], " vs ", condsToCompare[2], sep=""))
    dev.off()
    
    png(filename = "report_data/heatmap_selected_orderBaseMean.png", res = 200, units = "in", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[1], " vs ", condsToCompare[2], sep=""))
    dev.off()
    
    svg(filename = "report_data/heatmap_selected_orderBaseMean.svg", height= 7, width = 7)
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df,
             main = paste("Top ", TOP, " significantly DE genes (log2norm)\n", condsToCompare[1], " vs ", condsToCompare[2], sep=""))
    dev.off()
    
  }
  
  #system("for i in heatmap_selected_*; do pdftk $i cat 2-end output tmp.pdf; mv tmp.pdf $i; done") # Cut first empty page from heatmap plot
  tryCatch({
    if(TOP > 1){
      # Plot counts for all significant genes
      dds_plot <- dds[,unique(count_dt[condition %in% condsToCompare,]$sample_name)]
      colData(dds_plot)$condition <- factor(colData(dds_plot)$condition,levels = unique(colData(dds_plot)$condition))
      pdf(file="all_sig_genes_normCounts.pdf")
      for(select_gene_name in comp_res[significant_DE == T]$Feature_name){
        plotCounts(dds_plot, gene=select_gene_name, intgroup="condition", main=select_gene_name)
        mtext(paste0("adj. p-value < ", p_value_threshold, ", logFC >= ", round(lfc_threshold,3)))
        # axis(1, at=seq_along(levels(coldata$condition)), levels(coldata$condition), las=2) # Ugly but works; I am not able to turn off axis() setting in plotCounts function
      }
      dev.off()
    }
  }, error=function(e){})
  
  # Write all normalized counts
  # fwrite(dcast.data.table(count_dt[condition %in% condsToCompare],Feature_name ~ sample_name,value.var = "normcounts"), file="detail_results/norm_counts.tsv", sep="\t")
  
  # # Make gene background/universe
  # background<-as.data.frame(res2[, "gene_name"])
  # colnames(background)<-"gene_name"
  # background$gene_id<-rownames(res2)
  # write.table(background, file="detail_results/background.txt", row.names=F, quote = F, sep="\t") # Write background file
  
  setwd(orig_dir)
}








