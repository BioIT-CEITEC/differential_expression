comparison_specific_edgeR_DESeq2_overlap <- function(output_dir,comp_res,edgeR_comp_res,p_value_threshold,lfc_threshold){
  orig_dir <- getwd()
  dir.create(paste0(output_dir,"/report_data"),showWarnings = F,recursive = T)
  setwd(output_dir)

  
  
  vennTable <- merge(comp_res[,.(Feature_name,DESeq2 = significant_DE * sign(log2FoldChange),DESeq2_padj = padj,DESeq2_log2FoldChange = log2FoldChange)],
                        edgeR_comp_res[,.(Feature_name,edgeR = significant_DE * sign(log2FoldChange),edgeR_padj = padj,edgeR_log2FoldChange = log2FoldChange)],by = "Feature_name")
  
  upsetTable <- merge(comp_res[,.(Ensembl_Id, DESeq2_up = ifelse(significant_DE * sign(log2FoldChange) == 1,1,0),
                                  DESeq2_down = ifelse(significant_DE * sign(log2FoldChange) == -1,1,0),
                                  DESeq2_no_filter_up = ifelse(no_filter_significant_DE * sign(log2FoldChange) == 1,1,0),
                                  DESeq2_no_filter_down = ifelse(no_filter_significant_DE * sign(log2FoldChange) == -1,1,0))],
                      edgeR_comp_res[,.(Ensembl_Id, edgeR_up = ifelse(significant_DE * sign(log2FoldChange) == 1,1,0),
                                        edgeR_down = ifelse(significant_DE * sign(log2FoldChange) == -1,1,0))],by = "Ensembl_Id")

  upsetTable[,test:=1]
  setcolorder(upsetTable,c("Ensembl_Id","test"))

  if(!is.null(dev.list())) {dev.off()}

  TestResultMatrix <- as.matrix(vennTable[,c("DESeq2","edgeR")],rownames = vennTable$Feature_name)
  TestResultMatrix <- abs(TestResultMatrix)
  
  pdf(file="overlap_DESeq2_edgeR_venn.pdf")
  par(oma=c(0,0,1,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  vennDiagram(vennCounts(TestResultMatrix), counts.col="black", main=paste("Overlap Between DESeq2 and edgeR results \nwith LogFC >= ",
                                                         round(lfc_threshold, 2), " and adj.pval ", p_value_threshold, sep=""))
  dev.off()

  png(file="overlap_DESeq2_edgeR_upset.png",width = 7,height = 7,units="in",res=300)
  UpSetR::upset(upsetTable,
                     nintersects = NA,
                     nsets = 6,
                     sets = c("DESeq2_down","DESeq2_no_filter_down","edgeR_down","edgeR_up","DESeq2_no_filter_up","DESeq2_up"),
                     keep.order = T,
                     text.scale = 1.2,
                     main.bar.color = "black"
                     )
  dev.off()
  
  # Add gene names
  vennTable[,DESeq2_mean_padj := (DESeq2_padj + edgeR_padj) / 2]
  setorder(vennTable,DESeq2_mean_padj,na.last = T)
  setcolorder(vennTable,c(1,2,5,8,3,6,4,7))

  fwrite(vennTable, file = "overlap_DESeq2_edgeR_de.tsv", sep="\t")

  setwd(orig_dir)
}

# ####################################################################################################
# # Overlap between DESeq2 and edgeR and Venn diagrams
# # http://www.ats.ucla.edu/stat/r/faq/venn.htm
# venn_list <- lapply(caller_types,function(x) var_tab[caller == x,paste(chrom,position,alternative,sep = "_")])
# venn_raw <- venn.diagram(venn_list,NULL
#                          ,imagetype = "tiff"
#                          ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
#                          ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
#                          ,cex = 1.5
#                          ,cat.fontface=2
#                          ,category.names=caller_types,main = "Variants all")
# 
# venn_raw_rel <- venn.diagram(venn_list,NULL,print.mode = "percent"
#                              ,imagetype = "tiff"
#                              ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
#                              ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
#                              ,cex = 1.5
#                              ,cat.fontface=2
#                              ,category.names=caller_types,main = "Variants all relative")
# 
# 
# venn_list_pass <- lapply(caller_types,function(x) var_tab[caller == x & is_pass == T,paste(chrom,position,alternative,sep = "_")])
# venn_pass <- venn.diagram(venn_list_pass,NULL
#                           ,imagetype = "tiff"
#                           ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
#                           ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
#                           ,cex = 1.5
#                           ,cat.fontface=2
#                           ,category.names=caller_types,main = "Variants pass")
# 
# venn_pass_rel <- venn.diagram(venn_list_pass,NULL,print.mode = "percent"
#                               ,imagetype = "tiff"
#                               ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
#                               ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
#                               ,cex = 1.5
#                               ,cat.fontface=2
#                               ,category.names=caller_types,main = "Variants pass relative")
# 
# 
# pdf(gsub(".tsv$",".variant_stats.pdf",output_file))
# grid.draw(venn_raw)
# grid.newpage()
# grid.draw(venn_pass)
# grid.newpage()
# grid.draw(venn_raw_rel)
# grid.newpage()
# grid.draw(venn_pass_rel)
# dev.off()
# 
# file.remove(list.files(".",pattern = "VennDiagram.*.log"))

