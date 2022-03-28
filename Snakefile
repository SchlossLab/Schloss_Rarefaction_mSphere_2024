# Snakefile
# Pat Schloss
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running data processing steps in rarefaction analysis study

datasets = ["bioethanol", "human", "lake", "marine","mice", "peromyscus", "rainforest", "rice", "seagrass",
            "sediment", "soil", "stream"]

rule targets:
  input:
    expand("data/{dataset}/data.fasta", dataset=datasets)
  shell:
      '''
      if $(ls | grep -q "mothur.*logfile"); then
          mkdir -p logs/mothur/
          mv mothur*logfile logs/mothur/
      fi
      echo "done."
      '''


rule silva:
  input:
    script="code/get_silva.sh"
  output:
    silva_align="data/references/silva.v4.align",
    silva_tax="data/references/silva.v4.tax"
  resources:  
    ncores=8
  shell:
    """
    {input.script} {resources.ncores}
    """

rule rdp:
  input:
    script="code/get_rdp.sh"
  output:
    rdp_fasta="data/references/trainset18_062020.pds.fasta",
    rdp_tax="data/references/trainset18_062020.pds.tax"
  shell:
    "{input.script}"


################################################################################
#
# Generate pruned versions of datasets based on the original assignments of
# sequences to each sample
#
################################################################################

# Run datasets { mice human soil marine etc. } through mothur pipeline through remove.lineage
rule clean_seqs:
  input: 
    sra="data/{dataset}/sra_info.tsv",
    silva_align=rules.silva.output.silva_align,
    rdp_fasta=rules.rdp.output.rdp_fasta,
    rdp_tax=rules.rdp.output.rdp_tax,
    download="code/datasets_download.sh",
    make_files="code/datasets_make_files.R",
    script="code/datasets_process.sh"
  output:
    "data/{dataset}/data.fasta",
    "data/{dataset}/data.count_table",
    "data/{dataset}/data.taxonomy"
  resources:  
    ncores=8,
    mem_mb=45000,
    time_min=3000
  shell:
    """
    {input.script} data/{wildcards.dataset} {resources.ncores} \
        {input.silva_align} {input.rdp_fasta} {input.rdp_tax}
    """
