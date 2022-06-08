#!/usr/bin/env bash

SHARED=$1
#SHARED=data/soil/data.otu.1.rshared
TEMP_SHARED=`echo $SHARED | sed -E "s/(.)shared/\1raw.\1shared/"`

cp $SHARED $TEMP_SHARED

mothur "#summary.single(shared=$TEMP_SHARED,
                        calc=nseqs-sobs-shannon-simpson-chao-ace-npshannon-coverage)"

GROUPS_SUMMARY=`echo $TEMP_SHARED | sed -E "s/.shared/groups.summary/"`
ALPHA=`echo $SHARED | sed -E "s/shared/_raw_alpha/"`
RABUND=`echo $TEMP_SHARED | sed -E "s/.shared/*rabund/"`

mv $GROUPS_SUMMARY $ALPHA
rm $RABUND $TEMP_SHARED
