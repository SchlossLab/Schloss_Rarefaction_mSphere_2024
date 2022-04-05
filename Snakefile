# Snakefile
# Pat Schloss
# Schloss Lab
# University of Michigan

# Purpose: Snakemake file for running data processing steps in rarefaction analysis study

datasets = ["bioethanol", "human", "lake", "marine","mice", "peromyscus", "rainforest", "rice", "seagrass",
            "sediment", "soil", "stream"]

seeds = list(range(1, 101))

rule targets:
  input:
    # expand("data/{dataset}/data.otu.shared", dataset=datasets),
    # expand("data/{dataset}/data.group_count", dataset=datasets)
    expand("data/{dataset}/data.otu.{seed}.rshared", dataset=datasets, seed=seeds)


rule silva:
  input:
    script="code/get_silva.sh"
  output:
    "data/references/silva.v4.align",
    "data/references/silva.v4.tax"
  resources:  
    job_name="silva",
    ncores=8,
  shell:
    """
    {input.script} {resources.ncores}
    """

rule rdp:
  input:
    script="code/get_rdp.sh"
  output:
    "data/references/trainset18_062020.pds.fasta",
    "data/references/trainset18_062020.pds.tax"
  resources:
    job_name="rdp"
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
    script="code/datasets_process.sh",
    download="code/datasets_download.sh",
    make_files="code/datasets_make_files.R",
    sra="data/{dataset}/sra_info.tsv",
    silva_align="data/references/silva.v4.align",
    rdp_fasta="data/references/trainset18_062020.pds.fasta",
    rdp_tax="data/references/trainset18_062020.pds.tax"
  output:
    "data/{dataset}/data.fasta",
    "data/{dataset}/data.count_table",
    "data/{dataset}/data.taxonomy"
  resources:
    job_name="{dataset}_clean_seqs"  
    ncores=8,
    mem_mb=45000,
    time_min=3000
  shell:
    """
    {input.script} data/{wildcards.dataset} {resources.ncores} \
        {input.silva_align} {input.rdp_fasta} {input.rdp_tax}
    """

# get the number of sequences in each group for each dataset
rule count_seqs:
  input:
    script="code/datasets_count_seqs.sh",
    count_table="data/{dataset}/data.count_table"
  output:
    "data/{dataset}/data.group_count"
  shell:
    """
    {input.script} data/{wildcards.dataset}
    """

# assign sequences to OTUs and generate a shared file
rule cluster_seqs:
  input:
    script="code/datasets_cluster.sh",
    fasta="data/{dataset}/data.fasta",
    count_table="data/{dataset}/data.count_table",
    tax="data/{dataset}/data.taxonomy"
  output:
    "data/{dataset}/data.otu.shared",
  resources:
    ncores=8,
    mem_mb=45000,
    time_min=3000
  shell:
    """
    {input.script} data/{wildcards.dataset} {resources.ncores}
    """


# generate null model shared files
rule null_shared:
  input:
    script="code/get_null_shared.R",
    shared="data/{dataset}/data.otu.shared"
  output:
    "data/{dataset}/data.otu.{seed}.rshared"
  shell:
    """
    {input.script} {input.shared} {seed}
    """

# non-rarefied shannon, sobs, invsimpson
# rarefy shannon, sobs, invsimpson
# estimated sobs with chao1, breakaway

# non-rarefied bray-curtis
# rarefied bray-curtis
# normalized bray-curtis
# relative abundance bray-curtis
# vst bray-curtis
# metagenomeseq bray-curtis