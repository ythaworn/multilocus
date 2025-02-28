#!/bin/bash


wd=$1
refname=$2
dref=$3
pop=$4
fbed=$5  # locus coord file


minDP=$6
maxDP=$7
minGQ=$8


dvcf=$wd/out/vcf/$pop
fin=$dvcf/joint.filt.minDP$minDP.maxDP$maxDP.minGQ$minGQ.vcf.gz

out=${fbed%.bed}

doutbase=$wd/out/loci/$out/$pop
dout=$doutbase/1.vcf_regions/minDP$minDP

mkdir -p $dout

# extract a vcf file for each locus
cat $dref/$fbed | xargs -n 3 -P 4 sh -c 'bcftools view -r $0:$1-$2 '"$fin -Oz -o $dout"'/${0}_${1}_${2}.vcf.gz'
 