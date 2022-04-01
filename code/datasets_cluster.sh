#!/usr/bin/env bash

DATA=$1
PROCESSORSSPLIT=$2
PROCESSORSCLUSTER=1

mothur "#set.seed(seed=19760620);
  set.dir(input=$DATA);
  cluster.split(fasta=data.fasta, count=data.count_table, taxonomy=data.taxonomy, taxlevel=4, 
                processors=$PROCESSORSSPLIT, runsensspec=FALSE, cluster=FALSE);
  cluster.split(file=current, count=current, processors=$PROCESSORSCLUSTER, runsensspec=FALSE);
  make.shared()"

mv $DATA/data.opti_mcc.shared $DATA/data.otu.shared
rm -f $DATA/data.opti_mcc.* $DATA/data.file
