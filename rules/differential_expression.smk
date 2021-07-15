
def final_variant_calling_report_input(wildcards):

    if config['conditions_to_compare'] == "all":
        condition_list = sorted(list(config['condition']))
    else:
        condition_list = config['conditions_to_compare'].split(",")

    comparison_dir_list = list()
    for condition1 in condition_list:
        if ':' in condition1:
            conditions = condition1.split(":")
            comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
        else:
            for condition2 in condition_list[condition_list.index(condition1):]:
                if ':' not in condition2 and condition2 != condition1:
                    comparison_dir_list.append(condition2 + "_vs_" + condition1)

    biotype_dir_list = config['biotypes'].split(",")
    return expand("results/{comparison}/{biotype}/edgeR.tsv", comparison=comparison_dir_list, biotype=biotype_dir_list)

rule final_DE_report:
    input: reports = final_variant_calling_report_input
    output: html = "results/" + PROJECT_NAME + "_finished"
    shell:
        "touch {output.html}"

def DE_computation_input(wildcards):
    if config["feature_count"]:
        analysis_subclass = "feature_count"
        return "results/complete.feature_count.tsv"
    else:
        analysis_subclass = "RSEM_count"
        return "results/complete.RSEM.RData"


rule DE_computation:
    input:  expression_tab=DE_computation_input,
            biotype_groups = expand("{ref_dir}/general/default/annot/biotypes_list_mod.txt",ref_dir=reference_directory)[0],
            sqlite = expand("{ref_dir}/annot/{ref}.sqlite.gz",ref_dir=reference_directory,ref=config["reference"])[0],
    output: table = "results/{comparison}/{biotype}/edgeR.tsv",
    params: config = "config.json",
            count_type = analysis_subclass,
            organism = config["organism"],
            use_tag_to_pair_samples = config["use_tag_to_pair_samples"],
            ref_from_trans_assembly = config["ref_from_trans_assembly"]
    log:    "logs/{comparison}/{biotype}/{comparison}.{biotype}.DE.run.log"
    conda:  "../wrappers/DE_computation/env.yaml"
    script: "../wrappers/DE_computation/script.py"

rule analysis_RSEM_table:
    input:  RSEM = expand("qc_reports/{sample}/rsem_count/{sample}.genes.results",sample=sample_tab.sample_name),
            gtf= expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: RSEM = "results/analysis_RSEM_table/complete.RSEM.RData"
    params: ref_from_trans_assembly = config["ref_from_trans_assembly"]
    log:    "logs/complete.RSEM.log"
    conda:  "../wrappers/analysis_RSEM_table/env.yaml"
    script: "../wrappers/analysis_RSEM_table/script.py"

rule analysis_feature_count_table:
    input:  feature_count = expand("qc_reports/feature_count/{sample}.feature_count.tsv",sample=sample_tab.sample_name)
    output: table = "results/analysis_feature_count_table/complete.feature_count.tsv",
    log:    "logs/complete.feature_count.log"
    conda:  "../wrappers/project_feature_count_table/env.yaml"
    script: "../wrappers/project_feature_count_table/script.py"
