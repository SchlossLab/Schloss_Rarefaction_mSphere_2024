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

rule rdp:
    input:
        script="code/get_rdp.bash"
    output:
        rdp_fasta="data/references/trainset18_062020.pds.fasta",
        rdp_tax="data/references/trainset18_062020.pds.tax"
    shell:
        "{input.script}"
        