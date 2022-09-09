#############################################################
# wrapper for rule: DE_computation
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: DE_computation \n##\n")
f.close()


snakemake.params.sample_tab.to_csv(snakemake.params.experiment_design, index=False,sep="\t")


command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/DE_computation_main.R "+\
            snakemake.params.experiment_design + " " +\
            snakemake.input.count_tab + " " +\
            ",".join(snakemake.params.comparison_dir_list) + " " +\
            snakemake.params.organism + " " +\
            snakemake.wildcards.analysis_type + " " +\
            str(snakemake.params.paired_replicates) + " " +\
            str(snakemake.params.normalize_data_per_comparison) + " " +\
            str(snakemake.params.use_custom_batch_effect_grouping) + " " +\
            str(snakemake.params.pvalue_for_viz) + " " +\
            str(snakemake.params.fold_change_threshold) + " " +\
            str(snakemake.params.named_in_viz) + " " +\
            str(snakemake.params.remove_genes_with_mean_read_count_threshold) + " " +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
# shell(command)
