# Snakefile
# Pat Schloss
# Schloss Lab
# University of Michigan

# Purpose:  Snakemake file for running data processing steps in rarefaction
#           analysis study

datasets = ["bioethanol", "human", "lake", "marine","mice", "peromyscus",
            "rainforest", "rice", "seagrass", "sediment", "soil", "stream"]

seeds = list(range(1, 101))

models = ["r"]

alpha_process = ["raw", "rarefy", "srs", "breakaway"]

beta_process = ["raw", "rare", "relabund", "srs", "metagenomeseq", "rclr",
                "zclr", "oclr", "nclr", "deseq2"]

beta_calculator = ["bray", "jaccard", "euclidean"]

designs = ["r", "s"]

rule targets:
  input:
    # expand("data/{dataset}/data.otu.shared", dataset=datasets),
    # expand("data/{dataset}/data.group_count", dataset=datasets),
    # expand("data/{dataset}/data.remove_accnos", dataset=datasets),
    # expand("data/{dataset}/data.otu.{seed}.rshared",
    #        dataset=datasets, seed=seeds),
    # expand("data/{dataset}/data.otu.{seed}.{model}_{process}_alpha",
    #        dataset=datasets, seed=seeds, model=models, process=alpha_process),
    # expand("data/{dataset}/data.otu.{seed}.{model}_{process}_{calculator}.dist",
    #        dataset=datasets, seed=seeds, model=models, process=beta_process,
    #        calculator=beta_calculator),
    # expand("data/{dataset}/data.{seed}.{design}design",
    #        dataset=datasets, seed=seeds, design=designs),
    # expand("data/{dataset}/data.{model}_{design}amova",
    #        dataset=datasets, model=models, design=designs)
    expand("data/{dataset}/data.{model}_{design}alpha_kw",
           dataset=datasets, model=models, design=designs)
           


rule silva:
  input:
    script="code/get_silva.sh"
  output:
    "data/references/silva.v4.align",
    "data/references/silva.v4.tax"
  resources:  
    cpus=8
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
rule null_shared:
  input:
    script="code/get_null_shared.R",
    shared="data/{dataset}/data.otu.shared",
    accnos="data/{dataset}/data.remove_accnos"
  resources:
    mem_mb=10000
  output:
    "data/{dataset}/data.otu.{seed}.rshared"
  shell:
    """
    {input.script} {input.shared} {input.accnos} {wildcards.seed}
    """


################################################################################
#
# Alpha diversity analysis
#
################################################################################

# non-rarefied nseqs, shannon, sobs, invsimpson, chao, ace, npshannon, coverage
rule raw_alpha:
  input:
    script="code/get_raw_alpha.sh",
    shared="data/{dataset}/data.otu.{seed}.{model}shared"
  output:
    "data/{dataset}/data.otu.{seed}.{model}_raw_alpha"
  shell:
    """
    {input.script} {input.shared}
    """

# rarefied nseqs, shannon, sobs, invsimpson, chao, ace, npshannon, coverage
rule rarefy_alpha:
  input:
    script="code/get_rarefy_alpha.sh",
    shared="data/{dataset}/data.otu.{seed}.{model}shared"
  output:
    "data/{dataset}/data.otu.{seed}.{model}_rarefy_alpha"
  shell:
    """
    {input.script} {input.shared}
    """

# srs normalized nseqs, shannon, sobs, invsimpson, chao, ace, npshannon, coverage
rule srs_alpha:
  input:
    script="code/get_srs_alpha.R",
    shared="data/{dataset}/data.otu.{seed}.{model}shared"
  output:
    "data/{dataset}/data.otu.{seed}.{model}_srs_alpha"
  resources:
    mem_mb=10000
  shell:
    """
    {input.script} {input.shared}
    """

# observed/estimated sobs with breakaway / poisson
rule breakaway_alpha:
  input:
    script="code/get_breakaway_alpha.R",
    shared="data/{dataset}/data.otu.{seed}.{model}shared"
  output:
    "data/{dataset}/data.otu.{seed}.{model}_breakaway_alpha"
  shell:
    """
    {input.script} {input.shared}
    """


################################################################################
#
# Beta diversity analysis
#
################################################################################

# non-rarefied bray-curtis / jclass / euclidean
rule process_beta:
  input:
    script="code/get_{beta_process}_beta.R",
    shared="data/{dataset}/data.otu.{seed}.{model}shared"
  output:
    dist="data/{dataset}/data.otu.{seed}.{model}_{beta_process}_{beta_calculator}.dist"
  resources:
    mem_mb=16000
  shell:
    """
    {input.script} {input.shared} {output.dist}
    """

################################################################################
#
# Assign samples to treatment groups based on null model or by size
#
################################################################################

rule null_design:
  input:
    script="code/get_null_design.R",
    group_count="data/{dataset}/data.group_count",
    remove_accnos="data/{dataset}/data.remove_accnos"
  output:
    "data/{dataset}/data.{seed}.rdesign"
  resources:
    mem_mb=2000
  shell:
    """
    {input.script} {input.group_count} {input.remove_accnos} {wildcards.seed}
    """


rule samplesize_design:
  input:
    script="code/get_samplesize_design.R",
    group_count="data/{dataset}/data.group_count",
    remove_accnos="data/{dataset}/data.remove_accnos"
  output:
    "data/{dataset}/data.{seed}.sdesign"
  resources:
    mem_mb=2000
  shell:
    """
    {input.script} {input.group_count} {input.remove_accnos} {wildcards.seed}
    """


################################################################################
#
# Test whether alpha and beta diversity data yield significant differences
# for each model approach and for each type of effect size
#
################################################################################

# alpha: calculate p-value based on size of sample
rule run_alpha_kw:
  input:
    script = "code/run_alpha_kw.R",
    alpha_files = expand("data/{dataset}/data.otu.{seed}.{model}_{process}_alpha",
                        seed=seeds, process=alpha_process,
                        allow_missing=True),
    design_files = expand("data/{dataset}/data.{seed}.{design}design",
                           seed=seeds, allow_missing=True)
  output:
    alpha = "data/{dataset}/data.{model}_{design}alpha_kw"
  shell:
    """
    {input.script} data/{wildcards.dataset} {output.alpha}
    """

# beta: calculate p-value based on size of sample
rule run_beta_amova:
  input:
    script="code/run_beta_analysis.R",
    dist_files = expand("data/{dataset}/data.otu.{seed}.{model}_{process}_{calculator}.dist",
                        seed=seeds, process=beta_process,
                        calculator=beta_calculator, allow_missing=True),
    design_files  = expand("data/{dataset}/data.{seed}.{design}design",
                           seed=seeds, allow_missing=True)
  output:
    amova="data/{dataset}/data.{model}_{design}amova"
  resources:
    mem_mb=12000
  shell:
    """
    {input.script} data/{wildcards.dataset} {output.amova}
    """
