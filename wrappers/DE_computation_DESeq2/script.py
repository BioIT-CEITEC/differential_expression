#############################################################
# wrapper for rule: DESeq2_computation
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: DESeq2_computation \n##\n")
f.close()

command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/run_deseq2.R " + \
          snakemake.input.dds + " " + \
          snakemake.input.count_data_normalized + " " + \
          snakemake.params.experiment_design + " " + \
          snakemake.wildcards.comparison + " " + \
          str(snakemake.params.pvalue_for_viz) + " " + \
          str(snakemake.params.fold_change_threshold) + " " + \
          str(snakemake.params.named_in_viz) + " " + \
          snakemake.output.output_dir + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
shell(command)
