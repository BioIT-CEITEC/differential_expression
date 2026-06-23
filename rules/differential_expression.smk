rule creat_salmon_map_table:
    input:  salmon = expand("qc_reports/{sample}/salmon_map/{sample}.salmon_map.sf",sample=sample_tab.sample_name),
            gtf= config["organism_gtf"],
    output: salmon = "DE_salmon_map/complete_salmon_map_table.RData"
    log:    "logs/DE/create_salmon_map_table.log"
    conda:  "../wrappers/analysis_salmon_table/env.yaml"
    script: "../wrappers/analysis_salmon_table/script.py"

rule creat_salmon_align_table:
    input:  salmon = expand("qc_reports/{sample}/salmon_aln/{sample}.salmon_aln.sf",sample=sample_tab.sample_name),
            gtf= config["organism_gtf"],
    output: salmon = "DE_salmon_align/complete_salmon_align_table.RData"
    log:    "logs/DE/create_salmon_align_table.log"
    conda:  "../wrappers/analysis_salmon_table/env.yaml"
    script: "../wrappers/analysis_salmon_table/script.py"

rule creat_kallisto_table:
    input:  kallisto = expand("qc_reports/{sample}/kallisto/{sample}.kallisto.tsv",sample=sample_tab.sample_name),
            gtf= config["organism_gtf"],
    output: kallisto = "DE_kallisto/complete_kallisto_table.RData"
    log:    "logs/DE/create_kallisto_table.log"
    conda:  "../wrappers/analysis_kallisto_table/env.yaml"
    script: "../wrappers/analysis_kallisto_table/script.py"

rule creat_RSEM_table:
    input:  RSEM = expand("qc_reports/{sample}/RSEM/{sample}.genes.results",sample=sample_tab.sample_name),
            gtf= config["organism_gtf"],
    output: RSEM = "DE_RSEM/complete_RSEM_table.RData"
    log:    "logs/DE/create_RSEM_table.log"
    conda:  "../wrappers/analysis_RSEM_table/env.yaml"
    script: "../wrappers/analysis_RSEM_table/script.py"

rule creat_featureCount_exon_table:
    input:  feature_count = expand("qc_reports/{sample}/featureCount_exon/{sample}.featureCount_exon.tsv",sample=sample_tab.sample_name)
    output: table = "DE_featureCount_exon/complete_featureCount_exon_table.tsv",
    params: type="featureCount",
    log:    "logs/DE/create_featureCount_exon_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_featureCount_gene_table:
    input:  feature_count = expand("qc_reports/{sample}/featureCount_gene/{sample}.featureCount_gene.tsv",sample=sample_tab.sample_name)
    output: table = "DE_featureCount_gene/complete_featureCount_gene_table.tsv",
    params: type="featureCount",
    log:    "logs/DE/create_featureCount_gene_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_featureCount_transcript_table:
    input:  feature_count = expand("qc_reports/{sample}/featureCount_transcript/{sample}.featureCount_transcript.tsv",sample=sample_tab.sample_name)
    output: table = "DE_featureCount_transcript/complete_featureCount_transcript_table.tsv",
    params: type="featureCount"
    log:    "logs/DE/create_featureCount_transcript_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_featureCount_3pUTR_table:
    input:  feature_count = expand("qc_reports/{sample}/featureCount_3pUTR/{sample}.featureCount_3pUTR.tsv",sample=sample_tab.sample_name)
    output: table = "DE_featureCount_3pUTR/complete_featureCount_3pUTR_table.tsv",
    params: type="featureCount"
    log:    "logs/DE/create_featureCount_3pUTR_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_featureCount_5pUTR_table:
    input:  feature_count = expand("qc_reports/{sample}/featureCount_5pUTR/{sample}.featureCount_5pUTR.tsv",sample=sample_tab.sample_name)
    output: table = "DE_featureCount_5pUTR/complete_featureCount_5pUTR_table.tsv",
    params: type="featureCount"
    log:    "logs/DE/create_featureCount_5pUTR_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_HTSeqCount_exon_table:
    input:  feature_count = expand("qc_reports/{sample}/HTSeqCount_exon/{sample}.HTSeqCount_exon.tsv",sample=sample_tab.sample_name)
    output: table = "DE_HTSeqCount_exon/complete_HTSeqCount_exon_table.tsv",
    params: type="HTSeqCount",
    log:    "logs/DE/create_HTSeqCount_exon_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_HTSeqCount_gene_table:
    input:  feature_count = expand("qc_reports/{sample}/HTSeqCount_gene/{sample}.HTSeqCount_gene.tsv",sample=sample_tab.sample_name)
    output: table = "DE_HTSeqCount_gene/complete_HTSeqCount_gene_table.tsv",
    params: type="HTSeqCount",
    log:    "logs/DE/create_HTSeqCount_gene_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_HTSeqCount_transcript_table:
    input:  feature_count = expand("qc_reports/{sample}/HTSeqCount_transcript/{sample}.HTSeqCount_transcript.tsv",sample=sample_tab.sample_name)
    output: table = "DE_HTSeqCount_transcript/complete_HTSeqCount_transcript_table.tsv",
    params: type="HTSeqCount"
    log:    "logs/DE/create_HTSeqCount_transcript_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_HTSeqCount_3pUTR_table:
    input:  feature_count = expand("qc_reports/{sample}/HTSeqCount_3pUTR/{sample}.HTSeqCount_3pUTR.tsv",sample=sample_tab.sample_name)
    output: table = "DE_HTSeqCount_3pUTR/complete_HTSeqCount_3pUTR_table.tsv",
    params: type="HTSeqCount"
    log:    "logs/DE/create_HTSeqCount_3pUTR_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_HTSeqCount_5pUTR_table:
    input:  feature_count = expand("qc_reports/{sample}/HTSeqCount_5pUTR/{sample}.HTSeqCount_5pUTR.tsv",sample=sample_tab.sample_name)
    output: table = "DE_HTSeqCount_5pUTR/complete_HTSeqCount_5pUTR_table.tsv",
    params: type="HTSeqCount"
    log:    "logs/DE/create_HTSeqCount_5pUTR_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

rule creat_mirbase_canonical_table:
    input:  feature_count = expand("qc_reports/{sample}/mirbase_canonical/{sample}.mirbase_canonical.tsv",sample=sample_tab.sample_name)
    output: table = "DE_mirbase_canonical/complete_mirbase_canonical_table.tsv",
    params: type="mirna"
    log:    "logs/DE/create_mirbase_canonical_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

def count_tab_input(wildcards):
    if wildcards.analysis_type.split("_")[0] == "featureCount" or wildcards.analysis_type.split("_")[0] == "HTSeqCount" or wildcards.analysis_type.split("_")[0] == "mirbase":
        suffix = "tsv"
    else:
        suffix = "RData"

    return expand("DE_{analysis_type}/complete_{analysis_type}_table.{suffix}",analysis_type=wildcards.analysis_type,suffix=suffix)[0]

# Original rule commented out - replaced by modular rules below
# rule DE_computation:
#     input:  count_tab = count_tab_input,
#             gtf= config["organism_gtf"], # defined in utilities
#     output: table = expand("DE_{{analysis_type}}/{comparison}/DESeq2.tsv", comparison=comparison_dir_list),
#
#     params: organism = config["organism"],
#             paired_replicates = config["paired_replicates"],
#             sample_tab = sample_tab,
#             experiment_design = "DE_{analysis_type}/DE_experiment_design.tsv",
#             comparison_dir_list = comparison_dir_list,
#             keep_not_compared_samples_for_normalization = config["keep_not_compared_samples_for_normalization"],
#             normalize_data_per_comparison = config["normalize_data_per_comparison"],
#             use_custom_batch_effect_grouping= config["use_custom_batch_effect_grouping"],
#             pvalue_for_viz= config["pvalue_for_viz"],
#             fold_change_threshold= config["fold_change_threshold"],
#             named_in_viz= config["named_in_viz"],
#             remove_genes_with_sum_read_count_threshold=config["remove_genes_with_sum_read_count_threshold"],
#             remove_genes_with_mean_read_count_threshold=config["remove_genes_with_mean_read_count_threshold"],
#             geneList= config["filter_geneList"],
#             keepGene=config["filter_keepGene"],
#             chrmList=config["filter_chrmList"],
#             keepChrm=config["filter_keepChrm"],
#     log:    "logs/DE/DE_{analysis_type}.log"
#     conda:  "../wrappers/DE_computation/env.yaml"
#     script: "../wrappers/DE_computation/script.py"


# ============================================================================
# New modular DE computation rules
# ============================================================================

# Rule: Create experiment design file from sample table
rule create_experiment_design:
    input:
        count_tab = count_tab_input
    output:
        experiment_design = "DE_{analysis_type}/DE_experiment_design.tsv"
    params:
        sample_tab = lambda wildcards: sample_tab,
        paired_replicates = config["paired_replicates"]
    log:
        "logs/DE/create_experiment_design_{analysis_type}.log"
    run:
        import pandas as pd
        df = params.sample_tab.copy()
        if params.paired_replicates:
            # Patient is determined by replicate ID
            if 'batch_group' in df.columns and config.get('use_custom_batch_effect_grouping', False):
                df['patient'] = df['batch_group']
            else:
                df['patient'] = df['replicate']
        else:
            # Each sample is its own patient
            df['patient'] = 'pat' + df['sample_name'].astype(str)
        df.to_csv(output.experiment_design, sep='\t', index=False)


# Rule: Load count data and create count_data_original/txi objects
rule load_count_data:
    input:
        count_tab = count_tab_input,
        gtf = config["organism_gtf"],
        experiment_design = "DE_{analysis_type}/DE_experiment_design.tsv"
    output:
        output_dir = directory("DE_{analysis_type}/loading_data"),
        count_data_original = "DE_{analysis_type}/loading_data/count_data_original.RDS",
        txi = "DE_{analysis_type}/loading_data/txi.RDS"
    params:
        geneList = config["filter_geneList"],
        keepGene = config["filter_keepGene"],
        chrmList = config["filter_chrmList"],
        keepChrm = config["filter_keepChrm"],
        remove_genes_with_sum_read_count_threshold = config["remove_genes_with_sum_read_count_threshold"],
        remove_genes_with_mean_read_count_threshold = config["remove_genes_with_mean_read_count_threshold"]
    log:
        "logs/DE/load_count_data_{analysis_type}.log"
    conda:
        "../wrappers/DE_computation_loading_data/env.yaml"
    script:
        "../wrappers/DE_computation_loading_data/script.py"


# Rule: Normalize samples and generate PCA/heatmap/MDS plots (all samples together)
# Always created, regardless of normalize_data_per_comparison setting
rule normalize_and_visualize:
    input:
        count_data_original = "DE_{analysis_type}/loading_data/count_data_original.RDS",
        txi = "DE_{analysis_type}/loading_data/txi.RDS",
        experiment_design = "DE_{analysis_type}/DE_experiment_design.tsv"
    output:
        output_dir = directory("DE_{analysis_type}/PCA"),
        dds = "DE_{analysis_type}/PCA/dds_normalized.RDS",
        count_data_normalized = "DE_{analysis_type}/PCA/count_data_normalized.RDS",
        edgeR_DGEList = "DE_{analysis_type}/PCA/edgeR_DGEList_normalized.RDS",
        edgeR_fit = "DE_{analysis_type}/PCA/edgeR_fit_normalized.RDS"
    params:
        condition_to_compare_vec = lambda wildcards: "|".join(condition_list),
        analysis_type = lambda wildcards: wildcards.analysis_type
    log:
        "logs/DE/normalize_and_visualize_{analysis_type}.log"
    conda:
        "../wrappers/DE_computation_PCA/env.yaml"
    script:
        "../wrappers/DE_computation_PCA/script.py"


# Rule: Prepare comparison-specific normalized data
# If normalize_data_per_comparison = FALSE: copy from PCA
# If normalize_data_per_comparison = TRUE: recompute normalization with subsetted data
rule prepare_comparison_data:
    input:
        pca_dds = "DE_{analysis_type}/PCA/dds_normalized.RDS",
        pca_count = "DE_{analysis_type}/PCA/count_data_normalized.RDS",
        pca_edger = "DE_{analysis_type}/PCA/edgeR_DGEList_normalized.RDS",
        pca_edger_fit = "DE_{analysis_type}/PCA/edgeR_fit_normalized.RDS",
        count_data_original = "DE_{analysis_type}/loading_data/count_data_original.RDS",
        txi = "DE_{analysis_type}/loading_data/txi.RDS",
        experiment_design = "DE_{analysis_type}/DE_experiment_design.tsv"
    output:
        output_dir = directory("DE_{analysis_type}/{comparison}/normalized"),
        dds = "DE_{analysis_type}/{comparison}/normalized/dds_normalized.RDS",
        count_data_normalized = "DE_{analysis_type}/{comparison}/normalized/count_data_normalized.RDS",
        edgeR_DGEList = "DE_{analysis_type}/{comparison}/normalized/edgeR_DGEList_normalized.RDS",
        edgeR_fit = "DE_{analysis_type}/{comparison}/normalized/edgeR_fit_normalized.RDS"
    params:
        normalize_per_comparison = config["normalize_data_per_comparison"],
        comparison_conditions = lambda wildcards: "_vs_".join(wildcards.comparison.split("_vs_")[::-1]),  # Reverse for contrast
        pvalue_for_viz = config["pvalue_for_viz"],
        fold_change_threshold = config["fold_change_threshold"],
        named_in_viz = config["named_in_viz"]
    log:
        "logs/DE/prepare_comparison_data_{analysis_type}_{comparison}.log"
    conda:
        "../wrappers/DE_computation_PCA/env.yaml"
    script:
        "../wrappers/DE_computation_PCA/prepare_comparison_data.py"


# Rule: Run DESeq2 for a specific comparison
rule deseq2_computation:
    input:
        dds = "DE_{analysis_type}/{comparison}/normalized/dds_normalized.RDS",
        count_data_normalized = "DE_{analysis_type}/{comparison}/normalized/count_data_normalized.RDS",
        experiment_design = "DE_{analysis_type}/DE_experiment_design.tsv"
    output:
        output_dir = directory("DE_{analysis_type}/{comparison}/DESeq2"),
        deseq2_tsv = "DE_{analysis_type}/{comparison}/DESeq2/DESeq2_{comparison}.tsv"
    params:
        pvalue_for_viz = config["pvalue_for_viz"],
        fold_change_threshold = config["fold_change_threshold"],
        named_in_viz = config["named_in_viz"]
    log:
        "logs/DE/deseq2_{analysis_type}_{comparison}.log"
    conda:
        "../wrappers/DE_computation_DESeq2/env.yaml"
    script:
        "../wrappers/DE_computation_DESeq2/script.py"


# Rule: Run edgeR for a specific comparison
rule edger_computation:
    input:
        edgeR_DGEList = "DE_{analysis_type}/{comparison}/normalized/edgeR_DGEList_normalized.RDS",
        edgeR_fit = "DE_{analysis_type}/{comparison}/normalized/edgeR_fit_normalized.RDS",
        count_data_normalized = "DE_{analysis_type}/{comparison}/normalized/count_data_normalized.RDS"
    output:
        output_dir = directory("DE_{analysis_type}/{comparison}/edgeR"),
        edger_tsv = "DE_{analysis_type}/{comparison}/edgeR/edgeR_{comparison}.tsv"
    params:
        pvalue_for_viz = config["pvalue_for_viz"],
        fold_change_threshold = config["fold_change_threshold"],
        named_in_viz = config["named_in_viz"]
    log:
        "logs/DE/edger_{analysis_type}_{comparison}.log"
    conda:
        "../wrappers/DE_computation_edgeR/env.yaml"
    script:
        "../wrappers/DE_computation_edgeR/script.py"


# Rule: Compare DESeq2 and edgeR results
rule deseq2_edger_overlap:
    input:
        deseq2_tsv = "DE_{analysis_type}/{comparison}/DESeq2/DESeq2_{comparison}.tsv",
        edger_tsv = "DE_{analysis_type}/{comparison}/edgeR/edgeR_{comparison}.tsv"
    output:
        output_dir = directory("DE_{analysis_type}/{comparison}/overlap")
    params:
        pvalue_for_viz = config["pvalue_for_viz"],
        fold_change_threshold = config["fold_change_threshold"]
    log:
        "logs/DE/overlap_{analysis_type}_{comparison}.log"
    conda:
        "../wrappers/DE_computation_merge_results/env.yaml"
    script:
        "../wrappers/DE_computation_merge_results/script.py"


# def final_variant_calling_report_input(wildcards):
#     input = {}
#     input['tsv'] = expand("DE_{{analysis_type}}/{comparison}/edgeR.tsv", comparison=comparison_dir_list)
#     return input


rule DE_report:
    input:
        deseq2_tsv = expand("DE_{analysis_type}/{comparison}/DESeq2/DESeq2_{comparison}.tsv", comparison=comparison_dir_list, analysis_type=analysis),
        edger_tsv = expand("DE_{analysis_type}/{comparison}/edgeR/edgeR_{comparison}.tsv", comparison=comparison_dir_list, analysis_type=analysis),
        overlap_dir = expand("DE_{analysis_type}/{comparison}/overlap", comparison=comparison_dir_list, analysis_type=analysis)
    output: html = "final_report.html"
    params: config = "DE_report.json"
    conda: "../wrappers/DE_report/env.yaml"
    log:    "logs/DE/DE_report.log"
    script: "../wrappers/DE_report/script.py"