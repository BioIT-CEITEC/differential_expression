#############################################################
# wrapper for rule: analysis_RSEM_table
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: analysis_RSEM_table \n##\n")
f.close()

command = "$CONDA_PREFIX/bin/Rscript "+os.path.abspath(os.path.dirname(__file__))+"/combine_RSEM_tabs.R "+\
            snakemake.output.RSEM + " " +\
            " ".join(snakemake.input.RSEM) + " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)