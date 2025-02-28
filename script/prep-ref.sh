#!/bin/bash


wd=$1
refname=$2
dref=$3

dscript=$wd/script


# inputs
ref=$dref/$refname.fa
gff=$dref/$refname.gff
frepeat=$dref/$refname.repeat.bed

# outputs
chromSizes=$dref/$refname.chrom.sizes
fcoding=$dref/$refname.coding.bed
fgene=$dref/$refname.gene.bed
ftranscript=$dref/$refname.transcript.bed

fintron=$dref/$refname.intron.bed
fnoncod_unmask=$dref/$refname.noncod.unmask.bed
fnoncod=$dref/$refname.noncod.bed


# 1) scaffold sizes and genome annotations
samtools faidx $ref
cut -f1-2 $ref.fai > $chromSizes


# 2) get coordinates of CDS
grep $'CDS\t' $gff | cut -f1,4,5 | sed 's, ,\t,g' | sort -k1,1 -k2,2n | uniq > $fcoding

grep $'gene\t' $gff | cut -f1,4,5 | sed 's, ,\t,g' | sort -k1,1 -k2,2n | uniq > $fgene

grep $'transcript\t' $gff | cut -f1,4,5 | sed 's, ,\t,g' | sort -k1,1 -k2,2n | uniq > $ftranscript


# 3) get coords for intron [genic non-CDS] & noncoding
bedtools subtract -a $ftranscript -b $fcoding > $fintron

bedtools complement -i $ftranscript -g $chromSizes > $fnoncod_unmask

bedtools subtract -a $fnoncod_unmask -b $frepeat > $fnoncod


# 4) select loci
# 4.1) coding loci (cds)
region=coding
minlen=100; maxlen=1000000; gaplen=2000

Rscript $dscript/DefineLocusToSelect.R \
  -d $dref/ \
  -i $refname.$region.bed \
  -l $minlen -m $maxlen -g $gaplen \
  -o $refname.$region.$minlen.$maxlen.$gaplen.bed
# echo $(wc -l $dref/$refname.$region.$minlen.$maxlen.$gaplen.bed)


# 4.2) intron (genic noncoding) loci
region=intron
minlen=100; maxlen=1000; gaplen=2000

Rscript $dscript/DefineLocusToSelect.R \
  -d $dref/ \
  -i $refname.$region.bed \
  -l $minlen -m $maxlen -g $gaplen \
  -o $refname.$region.$minlen.$maxlen.$gaplen.bed
# echo $(wc -l $dref/$refname.$region.$minlen.$maxlen.$gaplen.bed)


# 4.3) intergenic noncoding loci
region=noncod
minlen=100; maxlen=1000; gaplen=2000;

Rscript $dscript/DefineLocusToSelect.R \
  -d $dref/ \
  -i $refname.$region.bed \
  -l $minlen -m $maxlen -g $gaplen \
  -o $refname.$region.$minlen.$maxlen.$gaplen.bed
# echo $(wc -l $dref/$refname.$region.$minlen.$maxlen.$gaplen.bed)
