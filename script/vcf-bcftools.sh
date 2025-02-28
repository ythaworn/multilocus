#!/bin/bash

wd=$1
refname=$2
dref=$3
pop=$4
sample_info=$5  # for sex info
fploidy=$6

# vcf filter
minDP=$7
maxDP=$8
minGQ=$9

# half depth
minDP2=$(( $minDP / 2 ))

filter1="( FMT/DP >= ${minDP} & FMT/DP <= ${maxDP} & FMT/GQ >= ${minGQ} ) | ( FMT/GT=='"'0'"' & FMT/DP >= ${minDP2} ) | ( FMT/GT=='"'1'"' & FMT/DP >= ${minDP2} ) | ( FMT/GT=='"'2'"' & FMT/DP >= ${minDP2} )"

ref=$dref/$refname.fa
din=$wd/out/bam
dout=$wd/out/vcf 

# individuals from this pop
declare -a ids=( $( grep "$pop" $sample_info | cut -d' ' -f 1 ) )

# output dir
doutpop=$dout/$pop

mkdir -p $doutpop

# input bam files
bamflag=markdup
declare -a ff=()
for id in "${ids[@]}"; do
  ff+=( $din/$id.$bamflag.bam )
done

# ploidy file from joblist
fsample=$doutpop/sample.txt
echo -n "" > $fsample  # clear content if any

for id in "${ids[@]}"; do
  grep $id $sample_info | cut -d' ' -f -3 >> $fsample
done

fvcf=$doutpop/joint.vcf.gz
fmask=$doutpop/joint.mask.minDP$minDP.maxDP$maxDP.minGQ$minGQ.vcf.gz
ffilt=$doutpop/joint.filt.minDP$minDP.maxDP$maxDP.minGQ$minGQ.vcf.gz


# variant sites only (-v)
bcftools mpileup -B -Q $q -Ou -f $ref ${ff[@]} --annotate INFO/AD,FORMAT/AD,FORMAT/DP | bcftools call --annotate GQ -vmO z --ploidy-file $fploidy -S $fsample -o $fvcf


# markfilter
bcftools filter --include "type='snp'" $fvcf | bcftools filter --include "$filter1" --SnpGap 5 --set-GTs . -O z -o $fmask

# apply filter
bcftools view --apply-filters "PASS" -O z -o $ffilt $fmask

# index
tabix -p vcf $fvcf
tabix -p vcf $ffilt

# stats
bcftools stats -F $ref $ffilt > $doutpop/joint.filt.minDP$minDP.maxDP$maxDP.minGQ$minGQ.txt
