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

alpha_process = ["raw", "rarefy", "srs", "inext", "breakaway"]

beta_process = ["raw", "rare", "relabund", "srs", "metagenomeseq", "rclr",
                "zclr", "oclr", "nclr", "deseq2"]

beta_calculator = ["bray", "jaccard", "euclidean"]

rule targets:
  input:
    "submission/manuscript.pdf",
    "submission/manuscript.docx",
    "submission/response_to_reviewers.pdf",
    "submission/supplementary_material.pdf",
    "submission/marked_up.pdf"

################################################################################

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

# estimated sobs, shannon, invsimpson with iNEXT and estimates
rule inext_alpha:
  input:
    script="code/get_inext_alpha.R",
    shared="data/{dataset}/data.otu.{seed}.{model}shared"
  output:
    "data/{dataset}/data.otu.{seed}.{model}_inext_alpha"
  conda:
    "config_files/rare_inext.yml"
  shell:
    """
    {input.script} {input.shared}
    """

################################################################################
#
# Beta diversity analysis
#
################################################################################

ruleorder: process_rclr > process_beta

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

rule process_rclr:
  input:
    r_script = "code/get_rclr_beta.R",
    py_script = "code/run_gemelli.py",
    shared="data/{dataset}/data.otu.{seed}.{model}shared"
  output:
    dist="data/{dataset}/data.otu.{seed}.{model}_rclr_{beta_calculator}.dist"
  resources:
    mem_mb=24000
  conda:
    "config_files/rare_gemelli.yml"
  shell:
    """
    {input.r_script} {input.shared} {output.dist}
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

rule intraquantile_design:
  input:
    script="code/get_intraquantile_design.R",
    group_count="data/{dataset}/data.group_count",
    remove_accnos="data/{dataset}/data.remove_accnos"
  output:
    "data/{dataset}/data.{seed}.idesign"
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

# beta: calculate p-value
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

rule shared_rank:
  input:
    script = "code/get_matrix_rank.R",
    observed = "data/{dataset}/data.otu.shared",
    simulated = expand("data/{dataset}/data.otu.{seed}.{model}shared",
                      seed=seeds, model = ["c", "r", "e"],
                      allow_missing = True)
  output:
    "data/{dataset}/data.rank"
  shell:
    """
    {input.script} {wildcards.dataset}
    """

rule pool_ranks:
  input:
    script = "code/pool_ranks.R",
    ranks = expand("data/{dataset}/data.rank", dataset = datasets)
  output:
    "data/process/rank_fractions.tsv"
  shell:
    """
    {input.script}
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
    alpha = expand("data/{dataset}/data.r_{design}alpha_kw",
                   dataset = datasets, allow_missing = True),
    beta = expand("data/{dataset}/data.r_{design}amova",
                  dataset = datasets, allow_missing = True),
    script = "code/plot_false_positive_null_model_size.R"
  output:
    "figures/false_positive_null_model_{design}size.tiff"
  shell:
    """
    {input.script} {wildcards.design}
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
    alpha = expand("data/{dataset}/data.c_calpha_kw",
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
    ba = "data/human/data.otu.1.r_breakaway_alpha",
    inext = "data/human/data.otu.1.r_inext_alpha"
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

rule fig_1:
  input:
    "figures/alpha_beta_depth_correlation.tiff",
  output:
    "submission/figure_1.tiff",
  shell:
    """
    cp {input} {output}
    """

rule fig_2:
  input:
    "figures/false_positive_null_model.tiff",
  output:
    "submission/figure_2.tiff",
  shell:
    """
    cp {input} {output}
    """

rule fig_3:
  input:
    "figures/false_positive_null_model_isize.tiff",
  output:
    "submission/figure_3.tiff",
  shell:
    """
    cp {input} {output}
    """

rule fig_4:
  input:
    "figures/power_effect_model.tiff",
  output:
    "submission/figure_4.tiff",
  shell:
    """
    cp {input} {output}
    """

rule fig_5:
  input:
    "figures/power_cffect_model.tiff",
  output:
    "submission/figure_5.tiff",
  shell:
    """
    cp {input} {output}
    """

rule fig_6:
  input:
    "figures/intrasample_variation.tiff",
  output:
    "submission/figure_6.tiff",
  shell:
    """
    cp {input} {output}
    """

rule fig_7:
  input:
    "figures/coverage_plot.tiff",
  output:
    "submission/figure_7.tiff",
  shell:
    """
    cp {input} {output}
    """

rule fig_s1:
  input:
    "figures/seqs_per_sample.tiff",
  output:
    "submission/figure_s1.png",
  shell:
    """
    convert {input} {output}
    """

rule fig_s2:
  input:
    "figures/example_alpha_cor.tiff",
  output:
    "submission/figure_s2.png",
  shell:
    """
    convert {input} {output}
    """

rule fig_s3:
  input:
    "figures/example_beta_cor.tiff",
  output:
    "submission/figure_s3.png",
  shell:
    """
    convert {input} {output}
    """

rule fig_s4:
  input:
    "figures/false_positive_null_model_ssize.tiff",
  output:
    "submission/figure_s4.png",
  shell:
    """
    convert {input} {output}
    """


rule write_paper:
  input:
    "code/render_rmd.R",
    "submission/manuscript.Rmd",
    "submission/references.bib",
    "submission/asm.csl",
    #
    "data/process/study_summary_statistics.tsv",
    expand("data/{dataset}/data.otu.obs_coverage", dataset = datasets),
    expand("data/{dataset}/data.otu.rarefy_coverage", dataset = datasets),
    "submission/figure_1.tiff",
    "submission/figure_2.tiff",
    "submission/figure_3.tiff",
    "submission/figure_4.tiff",
    "submission/figure_5.tiff",
    "submission/figure_6.tiff",
    "submission/figure_7.tiff",
    "submission/figure_s1.png",
    "submission/figure_s2.png",
    "submission/figure_s3.png",
    "submission/figure_s4.png"
  output:
    "submission/manuscript.docx",
    "submission/manuscript.pdf",
    "submission/manuscript.md",
    "submission/manuscript.tex"
  shell: 
    """
    code/render_rmd.R submission/manuscript.Rmd
    rm -f submission/manuscript.log
    mv submission/manuscript.knit.md submission/manuscript.md
    """

rule supplemental_text:
  input:
    "submission/supplementary_material.Rmd",
    expand("submission/figure_s{fig_number}.png", fig_number = list(range(1,5))),
    "code/render_rmd.R"
  output:
    "submission/supplementary_material.pdf",
  shell: 
    """
    code/render_rmd.R submission/supplementary_material.Rmd
    mv submission/supplementary_material.knit.md submission/supplementary_material.md 
    rm -f submission/supplementary_material.log submission/supplementary_material.tex
    """

rule response_to_reviewers:
  input:
    rmd="submission/response_to_reviewers.Rmd",
    rscript="code/render_rmd.R"
  output:
    "submission/response_to_reviewers.pdf"
  shell:
    """
    {input.rscript} {input.rmd}
    rm -f submission/response_to_reviewers.knit.md
    """

rule track_changes:
  input:
    "submission/manuscript.tex"
  output:
    "submission/marked_up.pdf"
  shell:
    """
    git cat-file -p 9b21157:submission/manuscript.tex > submission/manuscript_old.tex
    latexdiff submission/manuscript_old.tex submission/manuscript.tex > submission/marked_up.tex
    R -e "tinytex::pdflatex('submission/marked_up.tex')"
    rm marked_up.log
    rm submission/marked_up.tex
    rm submission/manuscript_old.tex
    """
