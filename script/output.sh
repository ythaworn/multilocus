#!/bin/bash

wd=$1
minDP=$2
region=$3
dout=$4
aln_filt_flag=$5

mkdir -p $dout

mv $wd/out/dset/$region/4.multilocus_$aln_filt_flag/minDP$minDP/*.txt $dout/
