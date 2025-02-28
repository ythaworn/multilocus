#!/bin/bash
#
# prepare depth mask for bcftools consensus

wd=$1
refname=$2
dref=$3
d=$4  # depth cutoff


ref=$dref/$refname.fa
chromsize=$dref/$refname.chrom.sizes

ddepth=$wd/out/depth

dout=$wd/out/mask
mkdir -p $dout

cd $ddepth

for f in *.d$d.txt.gz; do
  fout=${f//.txt.gz/.mask.bed}
  
  if [ ! -f $fout.gz ]; then
    echo $fout

    bedtools complement -i $f -g $chromsize > $fout
    gzip $fout
    mv $fout.gz $dout/
  fi
done
