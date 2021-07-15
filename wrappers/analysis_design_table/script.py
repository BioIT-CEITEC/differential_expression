#############################################################
# wrapper for rule: filter_and_output_annot_germline_variants
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: filter_and_output_annot_germline_variants \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/filter_and_format_variants.R "+\
            snakemake.input.tsv+ " " +\
            snakemake.output.xls+ " " +\
            snakemake.wildcards.format+ " " +\
            snakemake.params.selected_genes+ " " +\
            snakemake.params.custom_annot_DBs+ " " +\
            snakemake.params.min_variant_frequency

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
