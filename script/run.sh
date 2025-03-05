#!/bin/bash


## set working directory ##
wd=PATH/to/working-directory
dscript=$wd/script

## sample name & pop label ##
fsample=$wd/data/sample.txt

## pop list & min number of sequences ##
fpop=$wd/data/pop.txt

## reference genome ##
refname=AatrE4c2Lp1M
dref=$wd/data/ref
fploidy=$dref/$refname.ploidy.txt

## directory containing input fastq files ##
ddata=$wd/data/fastq

## final output ##
dout=$wd/out/multilocus

## variant filtering ##
minDP=20   # min read depth
maxDP=300  # max read depth
minGQ=20   # min GQ score

## alignment filtering options ##
missing_prop_cutoff_row=0.5
missing_prop_cutoff_col=0.5
aln_filt_flag=f50


## START ##
timer=$SECONDS

flog=$wd/log.txt
echo -n "" > $flog


## read mapping ##
echo "read mapping"
while IFS=" " read -r sample_id sex pop; do
  $dscript/read-mapping-pe.sh $wd $refname $dref $ddata $sample_id >> $flog 2>&1
done < $fsample


## genotyping ##
echo "genotyping"
while IFS=" " read -r pop; do
  $dscript/vcf-bcftools.sh $wd $refname $dref $pop $fsample $fploidy $minDP $maxDP $minGQ >> $flog 2>&1
done < $fpop


echo "extracting loci"

## prepare reference ##
$dscript/prep-ref.sh $wd $refname $dref >> $flog 2>&1


## extract individual loci (vcf) ##
while IFS=" " read -r pop; do
  for region_full in $refname.coding.100.1000000.2000 $refname.intron.100.1000.2000 $refname.noncod.100.1000.2000; do
    fbed=$region_full.bed
    $dscript/extract-loci.sh $wd $refname $dref $pop $fbed $minDP $maxDP $minGQ >> $flog 2>&1
  done
done < $fpop


## create depth masks ##
minDP2=$(( $minDP / 2 ))
for d in $minDP2 $minDP; do
  $dscript/depth-mask.sh $wd $refname $dref $d >> $flog 2>&1
done


## convert vcf to fasta ##
ftmp=$wd/tmp.txt
while IFS=" " read -r pop; do
  for region_full in $refname.coding.100.1000000.2000 $refname.intron.100.1000.2000 $refname.noncod.100.1000.2000; do
    $dscript/vcf-to-fasta.sh $wd $refname $dref $pop $minDP $region_full >> $ftmp 2>&1
  done
done < $fpop
rm $ftmp


for region_full in $refname.coding.100.1000000.2000 $refname.intron.100.1000.2000 $refname.noncod.100.1000.2000; do
  ## collect fasta sequences into one file per locus ##
  $dscript/combine-fasta.sh $wd $minDP $region_full

  ## convert fasta to phylip ##
  $dscript/fasta-to-phylip.sh $wd $minDP $region_full

  ## filter phylip ##
  python $dscript/filter-phylip-loci.py $wd $refname $minDP $region_full $dref $missing_prop_cutoff_row $missing_prop_cutoff_col $aln_filt_flag >> $flog 2>&1

  ## OPTIONAL: make blocks ##
  python $dscript/make-blocks.py $wd $refname $minDP $region_full $dref $aln_filt_flag >> $flog 2>&1

  ## move output files to final location  ##
  $dscript/output.sh $wd $minDP $region_full $dout $aln_filt_flag
done

echo "We're done!"

timer=$(($SECONDS-timer))
printf "Time used: %02d:%02d:%02d\n" "$((timer/3600))" "$((timer/60%60))" "$((timer%60))"
