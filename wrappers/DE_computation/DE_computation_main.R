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
# parsed_ensembl - table with gene_id, gene_name and gene_biotype based in Ensembl
# biotypes - table with name and biotype_group where name equals gene biotype
####################################################################################################
# TODO: Add tSNE clustering
# TODO: Add shared DE genes visualization (inspiration at VBCF BioComp)
#####################################################################################################

#devtools::install_github("r-lib/later")
library(data.table)
library(RColorBrewer)
library(biomaRt)
library(ggplot2)
library(DESeq2)
library(cowplot)
library(gplots)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(edgeR)
library(limma)
library(svglite)
library(ashr)
library(plotly)
library(htmlwidgets)
library(rtracklayer)
library(UpSetR)

run_all <- function(args){
  
  experiment_design_file <- args[1]
  counts_file <- args[2]
  relative_output_dir <- dirname(experiment_design_file)
  comparison_vec <- strsplit(args[3],split = ",")[[1]]
  gtf_filename <- args[4]
  analysis_type <- args[5]
  paired_samples <- as.logical(toupper(args[6]))
  normalize_data_per_comparison <- as.logical(toupper(args[7]))
  use_custom_batch_effect_grouping <- as.logical(toupper(args[8]))
  p_value_threshold <- as.numeric(args[9])
  lfc_threshold <- log(as.numeric(args[10]),2)
  TOP <- as.integer(args[11])
  remove_genes_with_mean_read_count_threshold <- as.integer(args[12])
  INTERCEPT<-TRUE
  
  output_dir <- paste0(getwd(),"/",relative_output_dir)
  
  res_list <- read_and_prepare_design_data(comparison_vec,experiment_design_file,paired_samples,use_custom_batch_effect_grouping)
  experiment_design <- res_list[[1]]
  comparison_vec <- res_list[[2]]
  condition_to_compare_vec <- res_list[[3]]
  
  res_list <- read_and_prepare_count_data(counts_file,experiment_design,gtf_filename,analysis_type)
  count_dt <- res_list[[1]]
  txi <- res_list[[2]]
  
  if(!normalize_data_per_comparison){
    #DESeq2 part
    ################
    res_list <- DESeq2_computation(txi,count_dt,experiment_design,remove_genes_with_mean_read_count_threshold,condition_to_compare_vec)
    dds <- res_list[[1]]
    count_dt <- res_list[[2]]
    
    create_normalization_specific_DESeq2_results(paste0(output_dir,"/all_condition_results"),dds,count_dt)
    
    #edgeR part
    ################
    res_list <- edgeR_computation(txi,count_dt,experiment_design,condition_to_compare_vec = condition_to_compare_vec)
    edgeR_DGEList <- res_list[[1]]
    fit_tagwise_dispersion_DGEGLM <- res_list[[2]]
    
    create_normalization_specific_edgeR_results(paste0(output_dir,"/all_condition_results/edgeR"),edgeR_DGEList,count_dt)
  } 

  for(selected_comparison in comparison_vec){
    # to test
    # selected_comparison <- "WT_odd_vs_Mut_odd"
    condsToCompare <- strsplit(selected_comparison,"_vs_")[[1]]
    
    dir.create(paste0(output_dir,"/",selected_comparison), recursive = T,showWarnings = F)
    
    if(normalize_data_per_comparison){
      comparison_experiment_design <- experiment_design[condition %in% condsToCompare]
      
      #DESeq2 part
      ################
      res_list <- DESeq2_computation(txi,count_dt,experiment_design,remove_genes_with_mean_read_count_threshold,condition_to_compare_vec = condition_to_compare_vec)
      dds <- res_list[[1]]
      count_dt <- res_list[[2]]
      create_normalization_specific_DESeq2_results(paste0(output_dir,"/",selected_comparison),dds,count_dt[condition %in% condsToCompare])
      
      #edgeR part
      ################
      res_list <- edgeR_computation(txi,count_dt,experiment_design,condition_to_compare_vec = condition_to_compare_vec)
      edgeR_DGEList <- res_list[[1]]
      fit_tagwise_dispersion_DGEGLM <- res_list[[2]]
      
      create_normalization_specific_edgeR_results(paste0(output_dir,"/",selected_comparison),edgeR_DGEList,count_dt)
    } else {
      
      #DESeq2 part
      ################
      create_normalization_specific_DESeq2_results(paste0(output_dir,"/",selected_comparison,"/detail_results"),dds,count_dt[condition %in% condsToCompare])
      
      #edgeR part
      ################
      create_normalization_specific_edgeR_results(paste0(output_dir,"/",selected_comparison,"/edgeR"),edgeR_DGEList,count_dt[condition %in% condsToCompare],reduced_plot_design = T)
    }
    
    #DESeq2 part
    ################
    comp_res <- get_comparison_specific_DESeq2_table(dds,count_dt,experiment_design,condsToCompare,output_dir = paste0(output_dir,"/",selected_comparison),p_value_threshold,lfc_threshold)
    create_comparison_specific_DESeq2_results(comp_res,dds,count_dt,condsToCompare,output_dir = paste0(output_dir,"/",selected_comparison),paired_samples,TOP,p_value_threshold,lfc_threshold)
    
    #edgeR part
    ################
    res_list <- get_comparison_specific_edgeR_table(fit_tagwise_dispersion_DGEGLM,edgeR_DGEList,count_dt,condsToCompare,output_dir = paste0(output_dir,"/",selected_comparison,"/edgeR"),p_value_threshold,lfc_threshold)
    edgeR_comp_res <- res_list[[1]]
    lrt_tgw <- res_list[[2]]
    
    create_comparison_specific_edgeR_results(edgeR_comp_res,lrt_tgw,condsToCompare,output_dir = paste0(output_dir,"/",selected_comparison,"/edgeR"),paired_samples,TOP,p_value_threshold,lfc_threshold)
    
    #DESeq2 - edgeR overlap
    ################
    comparison_specific_edgeR_DESeq2_overlap(paste0(output_dir,"/",selected_comparison),comp_res,edgeR_comp_res,p_value_threshold,lfc_threshold)
  }
  
}

# develop and test 2
# setwd("/mnt/ssd/ssd_1/workspace/vojta/6303__differential_expression__DE_analysis_refactor__220929")
# args <- c("DE_feature_count_orig/DE_experiment_design.tsv","DE_feature_count_orig/complete_feature_count_table.tsv","CAF_vs_CTRL,NF_O_vs_CTRL,CAF_vs_NF_O","/mnt/ssd/ssd_3/references/homo_sapiens/GRCh38-p10/annot/GRCh38-p10.gtf","feature_count","True","False","False","0.05","2","20","0","")
# script.dir <- "wrappers/DE_computation/"

script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))
source(paste0(script.dir,"/DESeq2_edgeR_overlap_func.R"))
source(paste0(script.dir,"/DESeq2_func.R"))
source(paste0(script.dir,"/edgeR_func.R"))
source(paste0(script.dir,"/load_data_func.R"))

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)





# # Save session info, history and image
# sink("session_info.txt")
# sessionInfo()
# sink()
#
# save.image(file = "image.RData")
# savehistory(file = "history.Rhistory")










