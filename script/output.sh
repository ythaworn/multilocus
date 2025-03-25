#!/bin/bash

wd=$1
minDP=minDP$2
region=$3
dout=$4
aln_filt_flag=$5

mkdir -p $dout/$minDP

mv $wd/out/dset/$region/4.multilocus_$aln_filt_flag/$minDP/*.txt $dout/$minDP
