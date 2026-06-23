#############################################################
# wrapper for rule: prepare_comparison_data
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: prepare_comparison_data \n##\n")
f.close()

command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/prepare_comparison_data.R " + \
          snakemake.input.pca_dds + " " + \
          snakemake.input.pca_count + " " + \
          snakemake.input.pca_edger + " " + \
          snakemake.input.pca_edger_fit + " " + \
          snakemake.input.count_data_original + " " + \
          snakemake.input.txi + " " + \
          snakemake.input.experiment_design + " " + \
          str(snakemake.params.normalize_per_comparison) + " " + \
          snakemake.wildcards.comparison + " " + \
          snakemake.output.output_dir + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
shell(command)
