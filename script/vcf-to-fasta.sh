#!/bin/bash

wd=$1
refname=$2
dref=$3
pop=$4
d=$5
region_full=$6

ref=$dref/$refname.fa

# half depth for sex chromosome
d2=$(( $d / 2 ))

din=$wd/out/loci/$region_full/$pop
doutbase=$din/2.fasta

dmask=$wd/out/mask
dvcf=$din/1.vcf_regions/minDP$d

dout=$doutbase/minDP$d
mkdir -p $dout

# loop over loci
for f in $dvcf/*.vcf.gz; do
  fname=$(basename ${f%.vcf.gz})

  ch=$(echo $fname | cut -d'_' -f 1)
  start=$(echo $fname | cut -d'_' -f 2)
  end=$(echo $fname | cut -d'_' -f 3)
  fin=$dvcf/${ch}_${start}_${end}.vcf.gz

  fout=$dout/${ch}_${start}_${end}.fa

  # comment out if overwrite
  if [ ! -s "$fout" ]; then
    echo "${ch}_${start}_${end}"

    # index vcf
    tabix -f $fin

    # create a region file (temporary)
    fregion=$dout/region_${ch}_${start}_${end}.bed
    echo "$ch:$start-$end" > $fregion

    # store screen output from bcftools consensus
    ftmp=$dout/tmp_${ch}_${start}_${end}.txt

    # get a list of samples in the vcf file
    sample_list=($(bcftools query -l $fin))
  
    # create an emppty output file
    echo -n "" > $fout

    # loop over samples in the input vcf
    for sample in "${sample_list[@]}"; do
      #echo $sample

      # get ploidy info from genotype call
      gt=$( bcftools query -f'[%GT\n]' -s $sample $fin | head -1 )

      # low-depth regions to mask
      if [[ $gt != *"/"* && $gt != *"|"* ]]; then
        # half cutoff if haploid
        fmask=$dmask/$sample.d$d2.mask.bed.gz
      else
        fmask=$dmask/$sample.d$d.mask.bed.gz
      fi

      # get sequence from ref
      fout_sample=$dout/${ch}_${start}_${end}_${sample}.fa

      samtools faidx -r $fregion $ref | \
        bcftools consensus $fin \
        -s $sample \
        -o $fout_sample \
        -I \
        -M 'N' \
        -m $fmask --mask-with '-' 2> $ftmp

      # rename sequence
      sed "s/>.*/>${ch}_${start}_${end}^$sample/" $fout_sample >> $fout
      
      rm $fout_sample
    done

    # clean up
    rm $fregion $ftmp
  fi

  # clean up
  if [ -f $fin.tbi ]; then 
    rm $fin.tbi
  fi
done
