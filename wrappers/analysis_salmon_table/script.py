#############################################################
# wrapper for rule: analysis_salmon_table
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: analysis_salmon_table \n##\n")
f.close()

command = "$CONDA_PREFIX/bin/Rscript "+os.path.abspath(os.path.dirname(__file__))+"/combine_salmon_tabs.R "+\
            snakemake.output.salmon + " " +\
            " ".join(snakemake.input.salmon) + " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)