#!/usr/bin/env bash

mkdir -p data/references/rdp

wget -N -P data/references/ https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset18_062020.pds.tgz

tar xvzf data/references/trainset18_062020.pds.tgz -C data/references/rdp
mv data/references/rdp/trainset18_062020.pds/trainset18_062020.* data/references

#Garbage collection
rm -rf data/references/rdp data/references/trainset18_062020.pds.tgz
