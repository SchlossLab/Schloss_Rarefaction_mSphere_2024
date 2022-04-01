#!/usr/bin/env bash

DATA=$1

mothur "#count.groups(count=$DATA/data.count_table)"

mv $DATA/data.count.summary $DATA/data.group_count
