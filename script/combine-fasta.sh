#!/bin/bash

wd=$1
mindp=minDP$2
region_full=$3

locus_info=$wd/data/ref/$region_full.bed
dinbase=$wd/out/loci/$region_full
dout=$wd/out/dset/$region_full/1.fasta/$mindp

mkdir -p $dout

# loop over loci
while IFS=$'\t' read -r scaffold begin end ; do
  locus=${scaffold}_${begin}_${end}

  fout=$dout/$locus.fa
  cat $dinbase/*/2.fasta/$mindp/$locus.fa > $fout
done < $locus_info
