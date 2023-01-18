rule creat_salmon_map_table:
    input:  salmon = expand("qc_reports/{sample}/salmon_map/{sample}.salmon_map.sf",sample=sample_tab.sample_name),
            gtf= expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: salmon = "DE_salmon_map/complete_salmon_map_table.RData"
    log:    "logs/DE/create_salmon_map_table.log"
    conda:  "../wrappers/analysis_salmon_table/env.yaml"
    script: "../wrappers/analysis_salmon_table/script.py"

rule creat_salmon_align_table:
    input:  salmon = expand("qc_reports/{sample}/salmon_aln/{sample}.salmon_aln.sf",sample=sample_tab.sample_name),
            gtf= expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: salmon = "DE_salmon_align/complete_salmon_align_table.RData"
    log:    "logs/DE/create_salmon_align_table.log"
    conda:  "../wrappers/analysis_salmon_table/env.yaml"
    script: "../wrappers/analysis_salmon_table/script.py"

rule creat_kallisto_table:
    input:  kallisto = expand("qc_reports/{sample}/kallisto/{sample}.kallisto.tsv",sample=sample_tab.sample_name),
            gtf= expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: kallisto = "DE_kallisto/complete_kallisto_table.RData"
    log:    "logs/DE/create_kallisto_table.log"
    conda:  "../wrappers/analysis_kallisto_table/env.yaml"
    script: "../wrappers/analysis_kallisto_table/script.py"

rule creat_RSEM_table:
    input:  RSEM = expand("qc_reports/{sample}/RSEM/{sample}.genes.results",sample=sample_tab.sample_name),
            gtf= expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: RSEM = "DE_RSEM/complete_RSEM_table.RData"
    log:    "logs/DE/create_RSEM_table.log"
    conda:  "../wrappers/analysis_RSEM_table/env.yaml"
    script: "../wrappers/analysis_RSEM_table/script.py"

rule creat_feature_count_table:
    input:  feature_count = expand("qc_reports/{sample}/feature_count/{sample}.feature_count.tsv",sample=sample_tab.sample_name)
    output: table = "DE_feature_count/complete_feature_count_table.tsv",
    log:    "logs/DE/create_feature_count_table.log"
    conda:  "../wrappers/analysis_feature_count_table/env.yaml"
    script: "../wrappers/analysis_feature_count_table/script.py"

def count_tab_input(wildcards):
    if wildcards.analysis_type == "feature_count":
        suffix = "tsv"
    else:
        suffix = "RData"

    return expand("DE_{analysis_type}/complete_{analysis_type}_table.{suffix}",analysis_type=wildcards.analysis_type,suffix=suffix)[0]

rule DE_computation:
    input:  count_tab = count_tab_input,
            gtf= expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0]
    output: table = expand("DE_{{analysis_type}}/{comparison}/DESeq2.tsv", comparison=comparison_dir_list),

    params: organism = config["organism"],
            paired_replicates = config["paired_replicates"],
            sample_tab = sample_tab,
            experiment_design = "DE_{analysis_type}/DE_experiment_design.tsv",
            comparison_dir_list = comparison_dir_list,
            keep_not_compared_samples_for_normalization = config["keep_not_compared_samples_for_normalization"],
            normalize_data_per_comparison = config["normalize_data_per_comparison"],
            use_custom_batch_effect_grouping= config["use_custom_batch_effect_grouping"],
            pvalue_for_viz= config["pvalue_for_viz"],
            fold_change_threshold= config["fold_change_threshold"],
            named_in_viz= config["named_in_viz"],
            remove_genes_with_mean_read_count_threshold=config["remove_genes_with_mean_read_count_threshold"]
    log:    "logs/DE/DE_{analysis_type}.log"
    conda:  "../wrappers/DE_computation/env.yaml"
    script: "../wrappers/DE_computation/script.py"


# def final_variant_calling_report_input(wildcards):
#     input = {}
#     input['tsv'] = expand("DE_{{analysis_type}}/{comparison}/edgeR.tsv", comparison=comparison_dir_list)
#     return input


rule DE_report:
    input: tsv = expand("DE_{analysis_type}/{comparison}/DESeq2.tsv", comparison=comparison_dir_list,analysis_type=selected_analysis_types)
    output: html = "final_report.html"
    params: config = "DE_report.json"
    log:    "logs/DE/DE_report.log"
    script: "../wrappers/DE_report/script.py"