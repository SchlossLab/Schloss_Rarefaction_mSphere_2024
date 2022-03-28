#!/usr/bin/env bash

# This script assumes that DATA already exists and has the appropriate sra_info.tsv file in it

# Should come in as a directory to dataset like data/marine
DATA=$1
NPROC=$2
SILVA=$3
RDP_FASTA=$4
RDP_TAX=$5

# Download, split, and compress the sequene read files listed in $DATA/sra_info.tsv. This file is described
# in the script's header comment
code/datasets_download.sh $DATA


# Generate $DATA/data.files from downloaded file names and information in $DATA/sra_info.tsv
code/datasets_make_files.R $DATA


# Run data through standard mothur pipeline through the steps prior to clustering
mothur "#set.seed(seed=19760620);
	set.dir(inputdir=$DATA, output=$DATA);
	make.contigs(file=data.files, processors=$NPROC, maxambig=0, maxlength=275, maxhomop=8);
  unique.seqs(fasta=current, count=current);
	align.seqs(fasta=current, reference=$SILVA);
	screen.seqs(fasta=current, count=current, start=9, end=9583);
	filter.seqs(fasta=current, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=T);
	classify.seqs(fasta=current, count=current, reference=$RDP_FASTA, taxonomy=$RDP_TAX);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);"


# Clean up the output data files

# We want to keep $DATA/data.fasta, $DATA/data.count_table, $DATA/data.taxonomy for future steps
mv $DATA/data.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta $DATA/data.fasta
mv $DATA/data.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table $DATA/data.count_table
mv $DATA/data.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.pick.taxonomy $DATA/data.taxonomy

# Remove fastq files and all intermediate mothur files to keep things organized
rm $DATA/*fastq*
rm $DATA/*.contigs.*
rm $DATA/*.filter
rm $DATA/*files
