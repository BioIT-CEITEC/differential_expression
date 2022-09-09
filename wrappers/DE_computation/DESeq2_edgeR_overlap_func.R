comparison_specific_edgeR_DESeq2_overlap <- function(output_dir,comp_res,edgeR_comp_res,p_value_threshold,lfc_threshold){
  ####################################################################################################
  # Overlap between DESeq2 and edgeR and Venn diagrams
  # http://www.ats.ucla.edu/stat/r/faq/venn.htm
  
  vennTable <- merge(comp_res[,.(Feature_name,DESeq2 = significant_DE * sign(log2FoldChange),DESeq2_padj = padj,DESeq2_log2FoldChange = log2FoldChange)],
        edgeR_comp_res[,.(Feature_name,edgeR = significant_DE * sign(log2FoldChange),edgeR_padj = padj,edgeR_log2FoldChange = log2FoldChange)],by = "Feature_name")
  
  TestResultMatrix <- as.matrix(vennTable[,c(2,3)],rownames = vennTable$Feature_name)
  
  pdf(file="overlap_DESeq2_edgeR_venn.pdf")
  par(oma=c(0,0,1,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  vennDiagram(vennCounts(TestResultMatrix), counts.col="black", main=paste("Overlap Between DESeq2 and edgeR results \nwith LogFC >= ",
                                                         round(lfc_threshold, 2), " and adj.pval ", p_value_threshold, sep=""))
  dev.off()
  
  # Add gene names
  vennTable[,DESeq2_mean_padj := (DESeq2_padj + edgeR_padj) / 2]
  setorder(vennTable,DESeq2_mean_padj,na.last = T)
  setcolorder(vennTable,c(1,2,5,8,3,6,4,7))
  
  
  fwrite(vennTable, file = "overlap_DESeq2_edgeR_de.tsv", sep="\t")

}

