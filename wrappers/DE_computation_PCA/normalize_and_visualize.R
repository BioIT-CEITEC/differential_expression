####################################################################################################
# Script to normalize samples and generate PCA/heatmap/MDS plots
####################################################################################################

library(data.table)
library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(RColorBrewer)

args <- commandArgs(trailingOnly = T)
count_data_file <- args[1]
txi_file <- args[2]  # Can be empty for featureCount data
experiment_design_file <- args[3]
analysis_type <- args[4]
output_dir <- args[5]

# Source utilities from parent DE_computation folder
script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir, "/../DE_computation/DESeq2_func.R"))
source(paste0(script.dir, "/../DE_computation/edgeR_func.R"))

# Load count data
count_dt_original <- readRDS(count_data_file)
if(txi_file != "" && file.exists(txi_file)){
  txi <- readRDS(txi_file)
} else {
  txi <- NULL
}
experiment_design <- unique(count_dt_original[,.(sample_name, condition, patient)])
setorder(experiment_design, condition, patient)
condition_to_compare_vec <- unique(experiment_design$condition)

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# DESeq2 normalization and visualization
# ============================================================================
cat("Running DESeq2 normalization and visualization...\n")

res_list <- DESeq2_computation(txi, count_dt_original, experiment_design, condition_to_compare_vec)
dds <- res_list[[1]]
count_dt_normalized <- res_list[[2]]

# Save normalized DESeq2 object
saveRDS(dds, file.path(output_dir, "dds_normalized.RDS"))
saveRDS(count_dt_normalized, file.path(output_dir, "count_data_normalized.RDS"))

# Generate PCA plot using vst counts
pca_data <- plotPCA(DESeq2::vst(dds, blind=FALSE), intgroup=c("condition", "patient"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = patient)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(file.path(output_dir, "DESeq2_PCA.pdf"), pca_plot, width = 8, height = 6)
ggsave(file.path(output_dir, "DESeq2_PCA.png"), pca_plot, width = 8, height = 6, dpi = 300)

# Generate heatmap of top variable genes
topVarGenes <- head(order(rowVars(assay(DESeq2::vst(dds))), decreasing = TRUE), 50)
mat <- assay(DESeq2::vst(dds))[topVarGenes, ]
mat_scaled <- t(scale(t(mat)))

annotation_col <- data.frame(
  condition = factor(experiment_design$condition),
  row.names = experiment_design$sample_name
)
rownames(annotation_col) <- experiment_design$sample_name

pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 Variable Genes (DESeq2 VST)")
dev.off()

# ============================================================================
# edgeR normalization and visualization
# ============================================================================
cat("Running edgeR normalization and visualization...\n")

res_list <- edgeR_computation(txi, count_dt_original, experiment_design, analysis_type, condition_to_compare_vec)
edgeR_DGEList <- res_list[[1]]
fit_tagwise_dispersion_DGEGLM <- res_list[[2]]

# Save edgeR objects
saveRDS(edgeR_DGEList, file.path(output_dir, "edgeR_DGEList_normalized.RDS"))
saveRDS(fit_tagwise_dispersion_DGEGLM, file.path(output_dir, "edgeR_fit_normalized.RDS"))

# Generate MDS plot
pdf(file.path(output_dir, "edgeR_MDS.pdf"), width = 7, height = 6)
plotMDS(edgeR_DGEList,
        col = as.numeric(factor(experiment_design$condition)),
        pch = as.numeric(factor(experiment_design$patient)),
        main = "edgeR MDS Plot")
legend("topright",
       legend = levels(factor(experiment_design$condition)),
       col = 1:length(levels(factor(experiment_design$condition))),
       pch = 16)
dev.off()

# Generate heatmap of top differentially expressed genes
cpm_values <- cpm(edgeR_DGEList, log = TRUE)
topGenes <- order(rowVars(cpm_values), decreasing = TRUE)[1:50]
mat_edgeR <- cpm_values[topGenes, ]
mat_edgeR_scaled <- t(scale(t(mat_edgeR)))

pheatmap(mat_edgeR_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 Variable Genes (edgeR CPM)")

cat("Normalization and visualization complete.\n")
cat("Output saved to:", output_dir, "\n")
