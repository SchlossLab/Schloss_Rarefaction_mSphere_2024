# Snakefile
# Pat Schloss
# Schloss Lab
# University of Michigan

# Purpose:  Snakemake file for running data processing steps in rarefaction
#           analysis study

datasets = ["bioethanol", "human", "lake", "marine","mice", "peromyscus",
            "rainforest", "rice", "seagrass", "sediment", "soil", "stream"]

seeds = list(range(1, 101))

model = ["r"]

rule targets:
  input:
    expand("data/{dataset}/data.otu.shared", dataset=datasets),
    expand("data/{dataset}/data.group_count", dataset=datasets),
    expand("data/{dataset}/data.remove_accnos", dataset=datasets)
    #expand("data/{dataset}/data.otu.{seed}.rshared", dataset=datasets, seed=seeds)


rule silva:
  input:
    script="code/get_silva.sh"
  output:
    "data/references/silva.v4.align",
    "data/references/silva.v4.tax"
  resources:  
    job_name="silva",
    cpus=8,
  shell:
    """
    {input.script} {resources.cpus}
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
# Clean datasets and generate original shared file
#
################################################################################

# Run datasets { mice human soil marine etc. } through mothur pipeline through
# remove.lineage
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
    job_name="{dataset}_clean_seqs",  
    cpus=8,
    mem_mb=45000,
    time_min=3000
  shell:
    """
    {input.script} data/{wildcards.dataset} {resources.cpus} \
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
    cpus=8,
    mem_mb=45000,
    time_min=3000
  shell:
    """
    {input.script} data/{wildcards.dataset} {resources.cpus}
    """

# identify those samples that need to be removed because there are too few
# sequences
rule list_rare_samples:
  input:
    script="code/screen_group_counts.R",
    counts=expand("data/{dataset}/data.group_count", dataset=datasets)
  output:
    accnos=expand("data/{dataset}/data.remove_accnos", dataset=datasets)
  shell:
    """
    {input.script}
    """

################################################################################
#
# Generate simulated community data based on observed shared data
#
################################################################################


# generate null model shared files
# rule null_shared:
#   input:
#     script="code/get_null_shared.R",
#     shared="data/{dataset}/data.otu.shared",
#     accnos="data/{dataset}/data.remove_accnos"
#   output:
#     "data/{dataset}/data.otu.{seed}.rshared"
#   shell:
#     """
#     {input.script} {input.shared} {input.shared} {seed}
#     """


################################################################################
#
# Alpha diversity analysis
#
################################################################################

# non-rarefied nseqs, shannon, sobs, invsimpson, chao, ace, npshannon, coverage
# rule raw_alpha:
#   input:
#     script="code/get_raw_alpha.sh",
#     shared="data/{dataset}/data.otu.{seed}.{model}shared"
#   output:
#     "data/{dataset}/data.otu.{seed}.{model}_raw_alpha"
#   shell:
#     """
#     {input.script} {input.shared}
#     """

# rarefied nseqs, shannon, sobs, invsimpson, chao, ace, npshannon, coverage
# rule rarefy_alpha:
#   input:
#     script="code/get_rarefy_alpha.sh",
#     shared="data/{dataset}/data.otu.{seed}.{model}shared"
#   output:
#     "data/{dataset}/data.otu.{seed}.{model}_rarefy_alpha"
#   shell:
#     """
#     {input.script} {input.shared}
#     """

# estimated sobs with breakaway (multi methods?)
# rule breakaway_alpha:
#   input:
#     script="code/get_breakaway_alpha.R",
#     shared="data/{dataset}/data.otu.{seed}.{model}shared"
#   output:
#     "data/{dataset}/data.otu.{seed}.{model}_breakaway_alpha"
#   shell:
#     """
#     {input.script} {input.shared}
#     """


################################################################################
#
# Beta diversity analysis
#
################################################################################


# non-rarefied bray-curtis / jclass / euclidean
# rarefied bray-curtis / jclass / euclidean
# normalized bray-curtis / jclass / euclidean
# relative abundance bray-curtis / jclass / euclidean
# vst-deseq2 bray-curtis / jclass / euclidean
# vst-metagenomeseq bray-curtis / jclass / euclidean
# aitchison euclidean

# alpha: calculate p-value based on size of sample
# beta: calculate p-value based on size of sample
