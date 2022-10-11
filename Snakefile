# Snakefile
# Pat Schloss
# Schloss Lab
# University of Michigan

# Purpose:  Snakemake file for running data processing steps in rarefaction
#           analysis study

datasets = ["bioethanol", "human", "lake", "marine","mice", "peromyscus",
            "rainforest", "rice", "seagrass", "sediment", "soil", "stream"]

seeds = list(range(1, 101))

indexes = list(range(1, 101))

models = ["r"]

alpha_process = ["raw", "rarefy", "srs", "breakaway"]

beta_process = ["raw", "rare", "relabund", "srs", "metagenomeseq", "rclr",
                "zclr", "oclr", "nclr", "deseq2"]

beta_calculator = ["bray", "jaccard", "euclidean"]

designs = ["r", "s"]

rule targets:
  input:
    ## generate original data files
    # expand("data/{dataset}/data.otu.shared", dataset=datasets),
    # expand("data/{dataset}/data.group_count", dataset=datasets),
    # expand("data/{dataset}/data.remove_accnos", dataset=datasets),
    #
    #
    ## generate null model data files and analyze
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
    #        dataset=datasets, model=models, design=designs),
    # expand("data/{dataset}/data.{model}_{design}alpha_kw",
    #        dataset=datasets, model=models, design=designs),
    #
    ## effect size analysis
    # expand("data/{dataset}/data.otu.{seed}.{model}shared",
    #        dataset=datasets, seed=seeds, model="e"),
    # expand("data/{dataset}/data.otu.{seed}.{model}_{process}_alpha",
    #        dataset=datasets, seed=seeds, model="e", process=alpha_process),
    # expand("data/{dataset}/data.otu.{seed}.{model}_{process}_{calculator}.dist",
    #        dataset=datasets, seed=seeds, model="e", process=beta_process,
    #        calculator=beta_calculator),
    # expand("data/{dataset}/data.{model}_{design}amova",
    #        dataset=datasets, model="e", design="e"),
    # expand("data/{dataset}/data.{model}_{design}alpha_kw",
    #        dataset=datasets, model="e", design="e")
    #
    ## richness depletion analysis
    # expand("data/{dataset}/data.otu.{seed}.{model}shared",
    #        dataset=datasets, seed=seeds, model="c"),
    # expand("data/{dataset}/data.otu.{seed}.{model}_{process}_alpha",
    #        dataset=datasets, seed=seeds, model="c", process=alpha_process),
    ## expand("data/{dataset}/data.otu.{seed}.{model}_{process}_{calculator}.dist",
    ##        dataset=datasets, seed=seeds, model="s", process=beta_process,
    ##        calculator=beta_calculator),
    ## expand("data/{dataset}/data.{model}_{design}amova",
    ##        dataset=datasets, model="s", design="s"),
    # expand("data/{dataset}/data.{model}_{design}alpha_kw",
    #        dataset=datasets, model="c", design="c")
    #
    ## correlation between metric and sample size
    # expand("data/{dataset}/random_alpha_correlation.tsv",
    #        dataset = datasets)
    # expand("data/{dataset}/random_beta_correlation.tsv",
    #        dataset = datasets)
    #
    ## sensitivity between sampling effort and alpha and beta diversity metrics
    # expand("data/{dataset}/data.otu.{seed}.r_rare_alpha.summary",
    #        dataset = datasets, seed = seeds),
    # expand("data/{dataset}/data.otu.alpha_depth.summary",
    #        dataset = datasets),  
    # expand("data/{dataset}/data.otu.{index}.{seed}.r_rare_{calculator}.summary",
    #        dataset = datasets, index = indexes, seed = seeds, calculator = ["bray", "jaccard"]),
    # expand("data/{dataset}/data.otu.beta_depth.summary",
    #        dataset = datasets),
    # expand("data/{dataset}/data.otu.obs_coverage", dataset = datasets),
    # expand("data/{dataset}/data.otu.rarefy_coverage", dataset = datasets),
    #
    ## observed human dataset analysis
    # "data/human/data.otu.oshared",
    # "data/human/data.odesign",
    # expand("data/human/data.otu.o_{process}_alpha",
    #        process = alpha_process),
    # expand("data/human/data.otu.o_{process}_{calculator}.dist",
    #        process = beta_process, calculator= beta_calculator),
    # "data/human/data.o_oalpha_kw",
    # "data/human/data.o_oamova"
    expand("data/mice/data.otu.o_{process}_alpha",
           process = alpha_process),
    expand("data/mice/data.otu.o_{process}_{calculator}.dist",
           process = beta_process, calculator= beta_calculator),

    #
    ## summary data files
    # "data/process/study_summary_statistics.tsv",
    #  
    ## figures
    # "figures/alpha_beta_depth_correlation.tiff",
    # "figures/false_positive_null_model.tiff",
    # "figures/false_positive_null_model_size.tiff",
    # "figures/power_effect_model.tiff",
    # "figures/power_cffect_model.tiff",
    # "figures/intrasample_variation.tiff",
    # "figures/coverage_plot.tiff",
    # "figures/example_alpha_cor.tiff",
    # "figures/example_beta_cor.tiff",
    # "figures/seqs_per_sample.tiff"
    
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


rule effect_shared_design:
  input:
    script="code/get_effect_shared_design.R",
    shared="data/{dataset}/data.otu.shared",
    accnos="data/{dataset}/data.remove_accnos"
  resources:
    mem_mb=10000
  output:
    "data/{dataset}/data.{seed}.edesign",
    "data/{dataset}/data.otu.{seed}.eshared",
  shell:
    """
    {input.script} {input.shared} {input.accnos} {wildcards.seed}
    """


rule richness_shared_design:
  input:
    script="code/get_richness_shared_design.R",
    shared="data/{dataset}/data.otu.shared",
    accnos="data/{dataset}/data.remove_accnos"
  resources:
    mem_mb=10000
  output:
    "data/{dataset}/data.{seed}.cdesign",
    "data/{dataset}/data.otu.{seed}.cshared",
  shell:
    """
    {input.script} {input.shared} {input.accnos} {wildcards.seed}
    """


# rule observed_shared_design:
#   input:
#     script="code/get_observed_shared_design.R",
#     shared="data/human/data.otu.shared",
#     accnos="data/human/data.remove_accnos",
#     sra="data/human/sra_info.tsv"
#   resources:
#   output:
#     "data/human/data.otu.oshared",
#     "data/human/data.odesign"
#   shell:
#     """
#     {input.script}
#     """


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

# rule raw_oalpha:
#   input:
#     script="code/get_raw_alpha.sh",
#     shared="data/{dataset}/data.otu.oshared"
#   output:
#     "data/{dataset}/data.otu.o_raw_alpha"
#   shell:
#     """
#     {input.script} {input.shared}
#     """

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

# rule rarefy_oalpha:
#   input:
#     script="code/get_rarefy_alpha.sh",
#     shared="data/{dataset}/data.otu.oshared"
#   output:
#     "data/{dataset}/data.otu.o_rarefy_alpha"
#   shell:
#     """
#     {input.script} {input.shared}
#     """

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

# rule srs_oalpha:
#   input:
#     script="code/get_srs_alpha.R",
#     shared="data/{dataset}/data.otu.oshared"
#   output:
#     "data/{dataset}/data.otu.o_srs_alpha"
#   resources:
#     mem_mb=10000
#   shell:
#     """
#     {input.script} {input.shared}
#     """

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

# rule breakaway_oalpha:
#   input:
#     script="code/get_breakaway_alpha.R",
#     shared="data/{dataset}/data.otu.oshared"
#   output:
#     "data/{dataset}/data.otu.o_breakaway_alpha"
#   shell:
#     """
#     {input.script} {input.shared}
#     """


################################################################################
#
# Beta diversity analysis
#
################################################################################

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

# rule process_obeta:
#   input:
#     script="code/get_{beta_process}_beta.R",
#     shared="data/{dataset}/data.otu.oshared"
#   output:
#     dist="data/{dataset}/data.otu.o_{beta_process}_{beta_calculator}.dist"
#   resources:
#     mem_mb=16000
#   shell:
#     """
#     {input.script} {input.shared} {output.dist}
#     """

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
rule alpha_kw:
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

# rule oalpha_kw:
#   input:
#     script = "code/run_oalpha_kw.R",
#     alpha_files = expand("data/{dataset}/data.otu.o_{process}_alpha",
#                         process=alpha_process, allow_missing=True),
#     design_file = "data/{dataset}/data.odesign"
#   output:
#     alpha = "data/{dataset}/data.o_oalpha_kw"
#   shell:
#     """
#     {input.script}
#     """


# beta: calculate p-value based on size of sample
rule beta_amova:
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

# rule obeta_amova:
#   input:
#     script="code/run_obeta_analysis.R",
#     dist_files = expand("data/{dataset}/data.otu.o_{process}_{calculator}.dist",
#                         process=beta_process, calculator=beta_calculator,
#                         allow_missing=True),
#     design_file  = "data/{dataset}/data.odesign"
#   output:
#     amova="data/{dataset}/data.o_oamova"
#   resources:
#     mem_mb=12000
#   shell:
#     """
#     {input.script}
#     """



# alpha: calculate correlations with metric based on size of sample
rule alpha_cor:
  input:
    script = "code/alpha_correlation.R",
    alpha_files = expand("data/{dataset}/data.otu.{seed}.r_{process}_alpha",
                         seed=seeds, process=alpha_process, allow_missing=True)
  output:
    output = "data/{dataset}/random_alpha_correlation.tsv"
  shell:
    """
    {input.script} {wildcards.dataset}
    """
    
# beta: calculate correlations with distance based on size of sample
rule beta_cor:
  input:
    script = "code/beta_correlation.R",
    dist_files = expand("data/{dataset}/data.otu.{seed}.r_{process}_{calculator}.dist",
                        seed=seeds, process=beta_process,
                        calculator=beta_calculator, allow_missing=True),
  output:
    output = "data/{dataset}/random_beta_correlation.tsv"
  shell:
    """
    {input.script} {wildcards.dataset}
    """

################################################################################
#
# Test sensitivity between sampling effort and rarefied alpha and beta diversity 
# metrics
#
################################################################################

rule rare_alpha_depth:
  input:
    script = "code/rarefy_alpha_analysis.R",
    shared = "data/{dataset}/data.otu.{seed}.rshared"
  output:
    summary = "data/{dataset}/data.otu.{seed, \d+}.r_rare_alpha.summary"
  resources:  
    mem_mb=8000,
    time_min=7200
  shell:
    """
    {input.script} {input.shared} {output.summary}
    """

rule pool_alpha_depth:
  input:
    script = "code/pool_alpha_summary.R",
    summaries = expand("data/{dataset}/data.otu.{seed}.r_rare_alpha.summary",
                  seed=seeds, allow_missing=True),
  output:
    pool = "data/{dataset}/data.otu.alpha_depth.summary"
  shell:
    """
    {input.script} {input.summaries}
    """

rule rare_beta_depth:
  input:
    script = "code/rarefy_beta_analysis.R",
    shared = "data/{dataset}/data.otu.{seed}.rshared",
    rm_file = "data/{dataset}/data.remove_accnos",
    nseqs = "data/{dataset}/data.group_count"
  output:
    summary = "data/{dataset}/data.otu.{index, \d+}.{seed, \d+}.r_rare_{calculator}.summary"
  wildcard_constraints:
    calculator="jaccard|bray"
  resources:  
    mem_mb=8000,
    time_min=6000,
    cpus=1
  shell:
    """
    {input.script} {input.shared} {input.nseqs} {input.rm_file} {output.summary}
    """

rule pool_beta_depth:
  input:
    script = "code/pool_beta_summary.R",
    summaries = expand("data/{dataset}/data.otu.{index}.{seed}.r_rare_{calculator}.summary",
                  index=indexes, seed=seeds, calculator=["bray","jaccard"], allow_missing=True),
  output:
    pool = "data/{dataset}/data.otu.beta_depth.summary"
  shell:
    """
    {input.script} {input.summaries}
    """

rule observed_coverage:
  input:
    shared = "data/{dataset}/data.otu.shared",
    removal = "data/{dataset}/data.remove_accnos",
    script = "code/get_observed_coverage.R"
  output:
    summary = "data/{dataset}/data.otu.obs_coverage"
  shell:
    """
    {input.script} {input.shared} {input.removal}
    """

rule rarefy_coverage:
  input:
    script = "code/rarefy_coverage.R",
    shared = "data/{dataset}/data.otu.shared",
    remove_file = "data/{dataset}/data.remove_accnos"
  resources:  
    mem_mb=16000
  output:
    "data/{dataset}/data.otu.rarefy_coverage"
  shell:
    """
    {input.script} {input.shared} {input.remove_file}
    """


################################################################################
#
# Generate figures
#
################################################################################

rule plot_alpha_beta_depth_correlation:
  input:
    summary = expand("data/{dataset}/random_{alpha_beta}_correlation.tsv",
                     dataset = datasets, alpha_beta = ["alpha", "beta"]),
    script = "code/plot_alpha_beta_depth_correlation.R"
  output:
    "figures/alpha_beta_depth_correlation.tiff"
  shell:
    """
    {input.script}
    """

rule plot_false_positive_null_model:
  input:
    alpha = expand("data/{dataset}/data.r_ralpha_kw",
                   dataset = datasets),
    beta = expand("data/{dataset}/data.r_ramova",
                  dataset = datasets),
    script = "code/plot_false_positive_null_model.R"
  output:
    "figures/false_positive_null_model.tiff"
  shell:
    """
    {input.script}
    """

rule plot_false_positive_null_model_size:
  input:
    alpha = expand("data/{dataset}/data.r_salpha_kw",
                   dataset = datasets),
    beta = expand("data/{dataset}/data.s_ramova",
                  dataset = datasets),
    script = "code/plot_false_positive_null_model_size.R"
  output:
    "figures/false_positive_null_model_size.tiff"
  shell:
    """
    {input.script}
    """

rule plot_power_effect_model:
  input:
    alpha = expand("data/{dataset}/data.e_ealpha_kw",
                   dataset = datasets),
    beta = expand("data/{dataset}/data.e_eamova",
                  dataset = datasets),
    script = "code/plot_power_effect_model.R"
  output:
    "figures/power_effect_model.tiff"
  shell:
    """
    {input.script}
    """

rule plot_power_cffect_model:
  input:
    alpha = expand("data/{dataset}/data.e_ealpha_kw",
                   dataset = datasets),
    script = "code/plot_power_cffect_model.R"
  output:
    "figures/power_cffect_model.tiff"
  shell:
    """
    {input.script}
    """

rule plot_coverage:
  input:
    script = "code/plot_coverage.R",
    obs_files = expand("data/{dataset}/data.otu.obs_coverage", dataset = datasets),
    rarefy_files = expand("data/{dataset}/data.otu.rarefy_coverage", dataset = datasets),
  output:
    "figures/coverage_plot.tiff"
  shell:
    """
    {input.script}
    """
    
rule plot_example_alpha_cor:
  input:
    script = "code/plot_example_alpha_cor.R",
    raw = "data/human/data.otu.1.r_raw_alpha",
    rare = "data/human/data.otu.1.r_rarefy_alpha",
    srs = "data/human/data.otu.1.r_srs_alpha",
    ba = "data/human/data.otu.1.r_breakaway_alpha"
  output:
    "figures/example_alpha_cor.tiff"
  shell:
    """
    {input.script}
    """

rule plot_example_beta_cor:
  input:
    script = "code/plot_example_beta_cor.R",
    dist_files = [ "data/human/data.otu.1.r_deseq2_euclidean.dist",
                  "data/human/data.otu.1.r_nclr_euclidean.dist",
                  "data/human/data.otu.1.r_oclr_euclidean.dist",
                  "data/human/data.otu.1.r_rclr_euclidean.dist",
                  "data/human/data.otu.1.r_zclr_euclidean.dist",
                  "data/human/data.otu.1.r_rare_bray.dist",
                  "data/human/data.otu.1.r_rare_jaccard.dist",
                  "data/human/data.otu.1.r_raw_bray.dist",
                  "data/human/data.otu.1.r_raw_jaccard.dist",
                  "data/human/data.otu.1.r_relabund_bray.dist",
                  "data/human/data.otu.1.r_relabund_jaccard.dist",
                  "data/human/data.otu.1.r_metagenomeseq_bray.dist",
                  "data/human/data.otu.1.r_metagenomeseq_jaccard.dist",
                  "data/human/data.otu.1.r_srs_bray.dist",
                  "data/human/data.otu.1.r_srs_jaccard.dist" ]
  output:
    "figures/example_beta_cor.tiff"
  shell:
    """
    {input.script}
    """

rule plot_intrasample_variation:
  input:
    removal = expand("data/{dataset}/data.remove_accnos", dataset=datasets),
    nseqs = expand("data/{dataset}/data.group_count", dataset=datasets),
    alpha = expand("data/{dataset}/data.otu.alpha_depth.summary", dataset=datasets),
    beta = expand("data/{dataset}/data.otu.beta_depth.summary", dataset=datasets),
    script = "code/plot_intrasample_variation.R"
  output:
    "figures/intrasample_variation.tiff"
  shell:
    """
    {input.script}
    """

rule plot_seqs_per_sample:
  input:
    removal = expand("data/{dataset}/data.remove_accnos", dataset=datasets),
    nseqs = expand("data/{dataset}/data.group_count", dataset=datasets),
    script = "code/plot_seqs_per_sample.R"
  output:
    "figures/seqs_per_sample.tiff"
  shell:
    """
    {input.script}
    """

rule pool_study_summary_statistics:
  input:
    script = "code/get_sample_summary_statistics.R",
    removal = expand("data/{dataset}/data.remove_accnos", dataset=datasets),
    nseqs = expand("data/{dataset}/data.group_count", dataset=datasets),
    sra = expand("data/{dataset}/sra_info.tsv", dataset=datasets)
  output:
    "data/process/study_summary_statistics.tsv"
  shell:
    """
    {input.script}
    """

################################################################################
#
# Submission related rules
#
################################################################################

rule write_paper:
  input:
    "submission/manuscript.Rmd",
    "submission/references.bib",
    "submission/asm.csl",
    #
    "data/process/study_summary_statistics.tsv",
    "figures/alpha_beta_depth_correlation.tiff",
    "figures/false_positive_null_model.tiff",
    "figures/false_positive_null_model_size.tiff",
    "figures/power_effect_model.tiff",
    "figures/power_cffect_model.tiff",
    "figures/intrasample_variation.tiff",
    "figures/coverage_plot.tiff",
    "figures/example_alpha_cor.tiff",
    "figures/example_beta_cor.tiff",
    "figures/seqs_per_sample.tiff"
  output:
    "submission/manuscript.pdf",
    "submission/manuscript.md",
    "submission/manuscript.tex"
  shell: 
    """
    R -e "library('rmarkdown'); render('submission/manuscript.Rmd', clean=FALSE)"
    mv submission/manuscript.knit.md submission/manuscript.md
    rm -f submission/manuscript.log
    """
    
