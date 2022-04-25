#!/usr/bin/env bash

SHARED=$1
#SHARED=data/soil/data.otu.1.rshared
mothur "#summary.single(shared=$SHARED,
                        calc=nseqs-sobs-shannon-simpson-chao-ace-npshannon-coverage)"

GROUPS_SUMMARY=`echo $SHARED | sed -E "s/.shared/groups.summary/"`
ALPHA=`echo $SHARED | sed -E "s/shared/_raw_alpha/"`
RABUND=`echo $SHARED | sed -E "s/.shared/*rabund/"`

mv $GROUPS_SUMMARY $ALPHA
rm $RABUND