#!/usr/bin/env bash

NPROC=${1:?ERROR: Need to define NPROC.}

wget -N -P data/references/silva/ https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v132.tgz
tar xvzf data/references/silva/silva.seed_v132.tgz -C data/references/silva

mothur "#get.lineage(fasta=data/references/silva/silva.seed_v132.align, taxonomy=data/references/silva/silva.seed_v132.tax, taxon=Bacteria);
    pcr.seqs(start=13862, end=23445, keepdots=F, processors="${NPROC}");
    degap.seqs();
    unique.seqs()"

grep ">" data/references/silva/silva.seed_v132.pick.pcr.ng.unique.fasta | cut -c 2-  > data/references/silva/silva.seed_v132.pick.pcr.ng.accnos

mothur "#get.seqs(fasta=data/references/silva/silva.seed_v132.pick.pcr.align, accnos=data/references/silva/silva.seed_v132.pick.pcr.ng.accnos);
    screen.seqs(minlength=240, maxlength=275, maxambig=0, maxhomop=8, processors="${NPROC}")"

mv data/references/silva/silva.seed_v132.pick.pcr.pick.good.align data/references/silva.v4.align

grep "^>" data/references/silva.v4.align | cut -c 2- > data/references/silva/silva.v4.accnos

mothur "#get.seqs(taxonomy=data/references/silva/silva.seed_v132.pick.tax, accnos=data/references/silva/silva.v4.accnos)"

mv data/references/silva/silva.seed_v132.pick.pick.tax  data/references/silva.v4.tax

rm -rf data/references/silva mothur*logfile #garbage collection
