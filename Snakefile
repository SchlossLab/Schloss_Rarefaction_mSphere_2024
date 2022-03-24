# Snakefile
# Pat Schloss
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running data processing steps in rarefaction analysis study

rule silva:
    input:
        script="code/get_silva.bash"
    output:
        silva_align="data/references/silva.v4.align",
        silva_tax="data/references/silva.v4.tax",
    resources:  
        ncores=4
    shell:
        "{input.script} {resources.ncores}"

#data/references/trainset16_022016.pds.% :
#	mkdir -p data/references/rdp
#	wget -N -P data/references/ https://mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz; \
#	tar xvzf data/references/Trainset16_022016.pds.tgz -C data/references/rdp;\
#	mv data/references/rdp/trainset16_022016.pds/trainset16_022016.* data/references;\
#	rm -rf data/references/rdp data/references/Trainset*
