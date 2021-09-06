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

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/mrna_de_counts.R "+\
            snakemake.input.cfg_tab+ " " +\
            snakemake.input.expression_tab+ " " +\
            os.path.dirname(snakemake.output.table)+ " " +\
            snakemake.input.sqlite+ " " +\
            snakemake.input.biotype_groups+ " " +\
            snakemake.wildcards.comparison+ " " +\
            snakemake.params.count_type + " " +\
            snakemake.params.organism + " " +\
            str(snakemake.params.use_tag_to_pair_samples) + " " +\
            str(snakemake.params.ref_from_trans_assembly) +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
