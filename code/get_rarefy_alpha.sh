#!/usr/bin/env bash

SHARED=$1
#SHARED=data/soil/data.otu.1.rshared
mothur "#summary.single(shared=$SHARED,
                        calc=nseqs-sobs-shannon-simpson-chao-ace-npshannon-coverage,
                        subsample=T)"

GROUPS_SUMMARY=`echo $SHARED | sed -E "s/.shared/groups.ave-std.summary/"`
ALPHA=`echo $SHARED | sed -E "s/shared/_rarefy_alpha/"`
RABUND=`echo $SHARED | sed -E "s/.shared/*rabund/"`

mv $GROUPS_SUMMARY $ALPHA
rm $RABUND