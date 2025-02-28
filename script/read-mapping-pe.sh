#!/bin/bash

wd=$1
refname=$2
dref=$3
din=$4
id=$5

ref=$dref/$refname.fa

# create index
bwa index $ref


dout=$wd/out/bam
ddepth=$wd/out/depth

mkdir -p $dout
mkdir -p $ddepth

# number of threads
CPU=4

f1=$din/${id}_1.fastq.gz
f2=$din/${id}_2.fastq.gz


fbam=$dout/$id.bam
fmdk=$dout/$id.markdup.bam
fmdk_txt=$dout/$id.markdup_metrics.txt
fsam=$dout/$id.sam.gz
ffixmate=$dout/$id.fixmate.bam

fcov=$dout/$id.coverage.txt

fd=$ddepth/$id.depth.txt
fd10=$ddepth/$id.d10.txt
fd20=$ddepth/$id.d20.txt
ftmp=$ddepth/$id.tmp.txt

if [ ! -f $fmdk ]; then
  bwa mem -t $CPU -R "@RG\\tID:${id}\\tSM:${id}\\tPL:Illumina" -M $ref $f1 $f2 | gzip -3 > $fsam
  samtools fixmate -O bam $fsam $ffixmate
  samtools sort -O bam -o $fbam -@ $CPU $ffixmate

  # mark duplicate
  gatk --java-options "-Xmx1g" MarkDuplicates -I $fbam -O $fmdk -M $fmdk_txt --REMOVE_DUPLICATES

  samtools index $fmdk
  samtools flagstat $fmdk > $dout/$id.flagstat.txt
  
  # mapped read depth
  genomeCoverageBed -bga -split -ibam $fmdk > $fd
  awk '$4 > 9'  $fd > $ftmp; bedtools merge -i $ftmp > $fd10
  awk '$4 > 19' $fd > $ftmp; bedtools merge -i $ftmp > $fd20

  gzip -f $fd
  gzip -f $fd10
  gzip -f $fd20

  # summary of read depth by chr
  samtools coverage $fmdk > $fcov

  # clean up
  rm $fsam $ffixmate $fbam $ftmp
fi
