
edgeR_computation<- function(txi = NULL,count_dt,experiment_design,analysis_type,condition_to_compare_vec){

  paired_samples <- length(unique(experiment_design$patient)) != nrow(experiment_design)
  
  if(!is.null(txi)){
    analysis_type <- "RSEM"
  } else {
    analysis_type <- "feature_count"
  }
  
  d<-DGEList(counts=count_matrix_from_dt(count_dt,condition_to_compare_vec=condition_to_compare_vec), group=unique(count_dt[,.(sample_name,condition)])$condition) # edgeR DGE object
  
  if(analysis_type == "feature_count"){
    d<-calcNormFactors(d) # Calculate normalization factors
  } else {
    
    normMat <- txi$length
    normMat <- as.data.frame(normMat/exp(rowMeans(log(normMat))))
    cts <- count_matrix_from_dt(count_dt,condition_to_compare_vec=condition_to_compare_vec)
    
    o <- log(calcNormFactors(cts/normMat, na.rm = T)) + log(colSums(cts/normMat, na.rm = T)) # Must not contain NA values
    d <- DGEList(cts, group=experiment_design[match(colnames(cts),sample_name)]$condition)
    d$offset <- t(t(log(normMat)) + o) # d is now ready for estimate dispersion functions see edgeR User's Guide
  }
  
  # Filtering of non-informative genes - very low count-per-million; recommended by edgeR vignette
  # http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf -
  # 4.1.3 Filter low expression tags
  keep<-rowSums(cpm(d)>1) >= 3 # Change to 1/3 of the samples? (length(colnames(d))/3)
  d<-d[keep, , keep.lib.sizes=FALSE]
  
  if(analysis_type == "feature_count"){
    d<-calcNormFactors(d) # Re-calculate normalization factors after filtering
  } else {
    # Recalculate cts and normMat to exclude removed genes
    cts<-cts[rownames(cts) %in% rownames(d$counts), ]
    normMat<-normMat[rownames(txi$length) %in% rownames(d$counts),]
    
    o <- log(calcNormFactors(cts/normMat, na.rm = T)) + log(colSums(cts/normMat, na.rm = T)) # Must not contain NA values
    d <- DGEList(cts, group=experiment_design[match(colnames(cts),sample_name)]$condition)
    d$offset <- t(t(log(normMat)) + o)
  }
  
  if(paired_samples==T){
    print("Using paired design in edgeR")
    design<-model.matrix(~patient+condition, data=experiment_design)
  }else{
    print("Using simple design in edgeR")
    design<-model.matrix(~condition, data=experiment_design)
  }
  
  # The same thing which follows can be done using >d<-estimateDisp(d, design) which should be also used if you plan to use QL methods bellow
  # https://support.bioconductor.org/p/79149/
  # You can use  estimateGLMRobustDisp() is you expect outliers in the data not associated with any sample or particular gene
  d<-estimateGLMCommonDisp(d, design) # Calculate GLM for common dispersion
  d<-estimateGLMTrendedDisp(d, design) # Calculate GLM for trended dispersion
  d<-estimateGLMTagwiseDisp(d, design) # Calculate GLM for tagwise dispersion
  
  # https://support.bioconductor.org/p/76790/
  fit_tgw<-glmFit(d, design, dispersion=d$tagwise.dispersion) # Fit tagwise dispersion;  fit_tgw<-glmQLFit(d, design, dispersion=d$tagwise.dispersion) can be used if the number of replicates is low; QL (glmQLFit and glmQLFTest) is more strict in the assumptions ~ increases adj.pvalues

  return(list(d,fit_tgw))
}


create_normalization_specific_edgeR_results <- function(output_dir,d,count_dt,reduced_plot_design = F){
  
  
  ##edgeR comparison specific graphic ploting
  ####################################################################################################
  
  #set output dir but remember where to return
  orig_dir <- getwd()
  # dir.create(paste0(output_dir,"/report_data"),showWarnings = F,recursive = T)
  dir.create(paste0(output_dir,"/detail_results"),showWarnings = F,recursive = T)
  setwd(output_dir)
  
  # get experiment_design_dt from count_dt
  experiment_design_dt <- unique(count_dt[,.(sample_name,condition,patient)])
  setorder(experiment_design_dt,condition,patient)
  
  #save sample names to control with lib.sizes
  dsamples<-as.data.table(d$samples, keep.rownames = T)
  setnames(dsamples, "rn", "sample_name")
  fwrite(merge(experiment_design_dt, dsamples, by.x=c("sample_name","condition"), by.y=c("sample_name","group"), all=F),"detail_results/edgeR_design_control.txt",sep = "\t")

  experiment_design_dt <- prepare_colors_and_shapes(experiment_design_dt)
  
  #set paired samples if more then one sample belongs to one patient(batch)
  paired_samples <- length(unique(experiment_design_dt$patient)) != nrow(experiment_design_dt)
  
  dfMDS <- as.data.table(plotMDS(d,plot = F)[c("x","y")])
  dfMDS[,sample_name := colnames(d$counts)]
  dfMDS[,condition := d$samples$group]
  dfMDS <- dfMDS[sample_name %in% experiment_design_dt$sample_name]
  
  mds <- ggplot(dfMDS, aes(x, y, color=condition)) +
    geom_point(size = 3) +
    scale_color_manual(values = unique(experiment_design_dt$cond_colours), name="") +
    theme_bw() +
    xlab("BCV distance 1") +
    ylab("BCV distance 2") +
    theme(legend.position="bottom") +
    ggtitle("MDSPlot_BCV_distance") +
    ggrepel::geom_text_repel(aes(x, y, label = sample_name), color="black") +
    theme(plot.title = element_text(face="bold"))
  
  ggsave("MDS_plot.pdf", mds, units = "in", width = 7, height = 7, dpi = 200)
  
  A<-aveLogCPM(d)
  d2<-d[A>1,]
  d2<-calcNormFactors(d2)
  logCPM<-cpm(d2, log=TRUE, prior.count=5)
  logCPM <- logCPM[,experiment_design_dt$sample_name]
  
  if(paired_samples == T){
    # See the effect of batch - should be included in the formula
    logCPMc<-removeBatchEffect(logCPM, experiment_design_dt$patient)
    
    dflogCPM <- as.data.table(plotMDS(logCPM,plot = F)[c("x","y")])
    dflogCPM[,sample_name := experiment_design_dt$sample_name]
    dflogCPM[,condition := experiment_design_dt$condition]

    dflogCPMc <- as.data.table(plotMDS(logCPMc,plot = F)[c("x","y")])
    dflogCPMc[,sample_name := experiment_design_dt$sample_name]
    dflogCPMc[,condition := experiment_design_dt$condition]

    mds.logcpm <- ggplot(dflogCPM, aes(x, y, color=condition)) +
      geom_point(size = 3) +
      scale_color_manual(values = unique(experiment_design_dt$cond_colours), name="") +
      theme_bw() +
      xlab("Leading logFC dim 1") +
      ylab("Leading logFC dim 2") +
      theme(aspect.ratio=1) +
      theme(legend.position="bottom") +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      ggtitle("MDS (logCPM)\nwithout sample pairing") +
      ggrepel::geom_text_repel(aes(x, y, label = sample_name), color="black") +
      theme(plot.title = element_text(face="bold"))
    
    mds.logcpmc <- ggplot(dflogCPMc, aes(x, y, color=condition)) +
      geom_point(size = 3) +
      scale_color_manual(values = unique(experiment_design_dt$cond_colours), name="") +
      theme_bw() +
      xlab("Leading logFC dim 1") +
      ylab("Leading logFC dim 2") +
      theme(aspect.ratio=1) +
      theme(legend.position="bottom") +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      ggtitle("MDS (logCPM)\nwith sample pairing") +
      ggrepel::geom_text_repel(aes(x, y, label = sample_name), color="black") +
      theme(plot.title = element_text(face="bold"))
    
    pmds <- plot_grid(mds.logcpm, mds.logcpmc, ncol = 2,align = "hv")
    
    ggsave("MDS_plot_batchEffect.pdf", pmds, units = "in", width = 7, height = 7, dpi=200)
  }
  
  ### Plot expression profiles
  expcolour <- rep(c(brewer.pal(8, "Dark2"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel2")),10)
  pdens <- ggplot(reshape2::melt(logCPM), aes(value, color=Var2)) +
    geom_density() +
    #theme_bw() +
    scale_color_manual(values = expcolour) +
    ggtitle("Expression profiles") +
    xlab("Log10 of normalized expression per gene (edgeR)") +
    ylab("Density") +
    theme(plot.title = element_text(face="bold")) +
    theme(legend.position="bottom") +
    labs(color = "")
  
  ggsave("detail_results/normalized_gene_expression_check.pdf", units = "in", width = 7, height = 7, dpi=200)
  
  if(!reduced_plot_design){
    # Plot BCV plot
    pdf(file="detail_results/BCV_plot.pdf")
    plotBCV(d, col.tagwise = "black")
    dev.off()
    
    # Plot tagwise mean variation
    pdf(file="detail_results/edgeR_mean_variation_tgwDisp.pdf")
    plotMeanVar(d, show.raw.vars=TRUE, show.tagwise.vars=TRUE, NBline=TRUE, main="Mean Variation Tagwise")
    dev.off()
  }

  setwd(orig_dir)
}


get_comparison_specific_edgeR_table <- function(fit_tgw,d,count_dt,condsToCompare,output_dir,p_value_threshold,lfc_threshold){
  #set output dir but remember where to return
  orig_dir <- getwd()
  dir.create(paste0(output_dir,"/detail_results"),showWarnings = F,recursive = T)
  setwd(output_dir)
  
  compared_samples <- rownames(d$sample)[d$samples$group %in% condsToCompare]
  
  ### Replacement of coef and contrast, should do the same
  if(!length(grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit_tgw$design)))){# If I cannot find condsToCompare[2] (usually intercept) set contrast as 1 for condsToCompare2
    cond_label <<- paste0(colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[1], "\\b"), colnames(fit_tgw$design))])
    my.contrasts <- makeContrasts(postvspre = cond_label, levels=fit_tgw$design) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
    # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
  }else{ # Else make proper contrast
    cond_label <<- paste0(colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[1], "\\b"), colnames(fit_tgw$design))], "-", colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit_tgw$design))])
    # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
    my.contrasts <- makeContrasts(postvspre = cond_label, levels=fit_tgw$design) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
  }
  colnames(my.contrasts)<-"contrast"
  lrt_tgw<-glmLRT(fit_tgw, contrast=my.contrasts[, "contrast"]) # If we have 3 conditions and want to compare 3 vs 2 we set contrast=c(0, -1, 1), if we want to compare 3 vs 1 or 2 vs 1 we just set coef=3 or coef=2, respectively; some more examples of contrast https://www.biostars.org/p/110861/
  
  resultsTbl.tgw<-as.data.table(topTags(lrt_tgw, n=nrow(lrt_tgw$table), adjust.method = "BH", sort.by = "p.value")$table,keep.rownames = T) # Extract all genes
  setnames(resultsTbl.tgw,"rn","Ensembl_Id")
  
  ####################################################################################################
  # Store row number that match between results and raw counts to include raw counts later in the
  #   output
  wh.rows.tgw<-match(resultsTbl.tgw$Ensembl_Id, rownames(d$counts))
  
  # Combine results with extracted DE genes, common dispersion, UpDown values, normalized counts and
  #   raw counts
  
  resultsTbl.tgw[,tgw.Disp := d$tagwise.dispersion[wh.rows.tgw]]
  resultsTbl.tgw <- cbind(resultsTbl.tgw,cpm(d)[wh.rows.tgw,compared_samples])
  setnames(resultsTbl.tgw,compared_samples,paste0(compared_samples,"_normCounts"))
  resultsTbl.tgw <- cbind(resultsTbl.tgw,d$counts[wh.rows.tgw,compared_samples])
  setnames(resultsTbl.tgw,compared_samples,paste0(compared_samples,"_rawCounts"))
  resultsTbl.tgw[,abs_log2FoldChange := abs(logFC)]

  setnames(resultsTbl.tgw,c("logFC","PValue","FDR"),c("log2FoldChange","pvalue","padj"))
  #resultsTbl.tgw[,UpDown := NULL]
  
  resultsTbl.tgw[,significant_DE := padj < p_value_threshold & abs_log2FoldChange > lfc_threshold]
  resultsTbl.tgw<-merge(unique(count_dt[,.(Ensembl_Id,Feature_name,biotype)]),resultsTbl.tgw,by="Ensembl_Id")
  resultsTbl.tgw[,baseMean := exp(logCPM)]
  setcolorder(resultsTbl.tgw,c("Ensembl_Id","baseMean","log2FoldChange","LR","pvalue","padj","significant_DE","Feature_name","biotype","logCPM","tgw.Disp","abs_log2FoldChange"))
  setorder(resultsTbl.tgw,padj,pvalue,-abs_log2FoldChange,-logCPM,na.last = T)

  fwrite(x = resultsTbl.tgw, file = "edgeR.tsv", sep="\t")
  
  setwd(orig_dir)
  
  return(list(resultsTbl.tgw,lrt_tgw))
  ####################################################################################################

}


create_comparison_specific_edgeR_results <- function(edgeR_comp_res,lrt_tgw,condsToCompare,output_dir,paired_samples,TOP,p_value_threshold,lfc_threshold){

  
  ##edgeR comparison specific graphic ploting
  ####################################################################################################
  
  #set output dir but remember where to return
  orig_dir <- getwd()
  dir.create(paste0(output_dir,"/detail_results"),showWarnings = F,recursive = T)
  setwd(output_dir)
  
  TOP <- min(TOP,sum(edgeR_comp_res$significant_DE))
  RANGE <- seq_len(TOP)
  
  if(TOP==0){ # Set range for naming the samples - help to avoid naming of samples even if there is no DE gene; THIS ERROR IS OK WHEN WE HAVE 0 DE GENES Error in text.default(res[RANGE, ]$log2FoldChange, -log(res[RANGE, ]$padj,  : zero-length 'labels' specified
    system("touch ZERO_DE_GENES_FOUND")
  }
  
  if(TOP==1){ # Set range for naming the samples - help to avoid naming of samples even if there is no DE gene; THIS ERROR IS OK WHEN WE HAVE 0 DE GENES Error in text.default(res[RANGE, ]$log2FoldChange, -log(res[RANGE, ]$padj,  : zero-length 'labels' specified
    system("touch ONLY_ONE_DE_GENE_FOUND")
  }

  edgeR_comp_res_summary <- edgeR_comp_res[, .(test = "edgeR",
               total = length(Feature_name),
               na = sum(is.na(padj)),
               not_sig = sum(!is.na(padj) & padj >= p_value_threshold),
               sig = sum(!is.na(padj) & padj < p_value_threshold),
               sig_up = sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange > 0),
               sig_down = sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange < 0),
               sig_lfc = sum(!is.na(padj) & padj < p_value_threshold & abs_log2FoldChange >= lfc_threshold),
               sig_lfc_up = sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange >= lfc_threshold),
               sig_lfc_down = sum(!is.na(padj) & padj < p_value_threshold & log2FoldChange <= -lfc_threshold))]

  fwrite(edgeR_comp_res_summary, "detail_results/edgeR_de_genes_summary.tsv", sep="\t")


  sink("detail_results/edgeR_de_genes_summary.txt")
  print(paste("Number of DE Genes With adj.pval < ", p_value_threshold, " Without LogFC Cut-off", sep=""))
  dt<-decideTestsDGE(lrt_tgw, adjust.method="BH", p.value=p_value_threshold)
  summary(dt)
  print(paste("Number of DE Genes With adj.pval < ", p_value_threshold, " and LogFC >= ", round(lfc_threshold, 2), sep=""))
  dt<-decideTestsDGE(lrt_tgw, adjust.method="BH", p.value=p_value_threshold, lfc=lfc_threshold)
  summary(dt)
  sink()
  
  
  ### EdgeR volcano
  
  if(TOP > 0){
    options(ggrepel.max.overlaps = 2*TOP)

    edgeR_comp_res[,sig := ifelse(padj < p_value_threshold, paste0("padj<", p_value_threshold), "Not Sig")]
    p <- ggplot(edgeR_comp_res[!is.na(padj)], aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col=sig), size=0.5) +
      scale_color_manual(values=c("black", "red")) +
      geom_text_repel(data=edgeR_comp_res[RANGE,], aes(label=Feature_name), size=3, max.overlaps = 2*TOP)+
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = c(lfc_threshold, -lfc_threshold), linetype = "longdash", colour="blue") +
      ggtitle(paste("Volcanoplot ",condsToCompare[1], " vs ", condsToCompare[2], " top ", TOP, " genes", sep="")) +
      annotate("text",x=min(edgeR_comp_res[!is.na(padj)]$log2FoldChange),y=0,vjust=-0.5,label=condsToCompare[2],size=5,fontface = "bold") +
      annotate("text",x=max(edgeR_comp_res[!is.na(padj)]$log2FoldChange),y=0,vjust=-0.5,label=condsToCompare[1],size=5,fontface = "bold") +
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + theme(plot.title = element_text(face="bold"))
    
    pdf(file=paste("edgeR_volcanoplot_", condsToCompare[1], "_vs_", condsToCompare[2],".pdf", sep=""))
    print(p)
    dev.off()
    
    ma <- ggmaplot(edgeR_comp_res, main =  paste0("MA plot ", condsToCompare[1], " vs ", condsToCompare[2], " top ", TOP, " genes (edgeR)"),
                   fdr = p_value_threshold, fc = 2^lfc_threshold, size = 0.5,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(edgeR_comp_res$Feature_name),
                   legend = "top", top = TOP,
                   font.label = c("bold", 11),
                   font.legend = "bold",
                   font.main = "bold",
                   ggtheme = ggplot2::theme_minimal())+
      theme(plot.title = element_text(hjust = 0.5))
    
    
    pdf(file=paste("edgeR_MAplot_", condsToCompare[1], "_vs_", condsToCompare[2],".pdf", sep=""))
    print(ma)
    dev.off()
  }
  
  setwd(orig_dir)
}
