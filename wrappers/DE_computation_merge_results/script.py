#############################################################
# wrapper for rule: DESeq2_edgeR_overlap
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: DESeq2_edgeR_overlap \n##\n")
f.close()

command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/run_comparison_overlap.R " + \
          snakemake.input.deseq2_tsv + " " + \
          snakemake.input.edger_tsv + " " + \
          snakemake.wildcards.comparison + " " + \
          str(snakemake.params.pvalue_for_viz) + " " + \
          str(snakemake.params.fold_change_threshold) + " " + \
          snakemake.output.output_dir + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
shell(command)
