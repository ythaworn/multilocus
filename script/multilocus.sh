#!/bin/bash
#
# Generate multilocus datasets from short-read sequence data
# 
# Yuttapong Thawornwattana, 21 Oct 2025


## set working directory ##
wd="PATH/to/working-directory"
dscript="$wd/script"

## sample name & pop label ##
fsample="$wd/data/sample.txt"

## pop list & min number of sequences ##
fpop="$wd/data/pop.txt"

## reference genome ##
refname="AatrE4c2Lp1M"
dref="$wd/data/ref"
ref="$dref/$refname.fa"
fploidy="$dref/$refname.ploidy.txt"

## directory containing input fastq files ##
ddata="$wd/data/fastq"

## final output ##
dout="$wd/out/output"

## optional: make blocks of loci ##
makeblock=true
fblock="$wd/data/block_info.txt"

## genotype filtering parameters ##
minDP=20   # min read depth
maxDP=300  # max read depth
minGQ=20   # min GQ score

## alignment filtering parameters ##
missing_prop_cutoff_row=0.5
missing_prop_cutoff_col=0.5
aln_filt_flag=f50

# number of CPUs to use for read mapping
CPU=4


#######################################
# Read mapping
# Globals:
#   CPU
# Arguments:
#   wd refname dref din id
# Returns:
#   bam and depth files
#######################################
read_mapping_pe() {
  wd="$1"
  refname="$2"
  ref="$3"
  din="$4"
  id="$5"

  # output directories
  local dout="$wd/out/bam"
  local ddepth="$wd/out/depth"

  # input files
  local f1="$din/${id}_1.fastq.gz"
  local f2="$din/${id}_2.fastq.gz"

  # outputs: bam files and summary statistics
  local fmdk="$dout/$id.markdup.bam"
  local fcov="$dout/$id.coverage.txt"
  local ffst="$dout/$id.flagstat.txt"

  # outputs: depth files
  local fd="$ddepth/$id.depth.txt"
  local fd10="$ddepth/$id.d10.txt"
  local fd20="$ddepth/$id.d20.txt"
  local ftmp="$ddepth/$id.tmp.txt"

  # index the reference
  bwa index "$ref"
  
  # create output directories
  mkdir -p "$dout"
  mkdir -p "$ddepth"

  if [[ ! -f "$fmdk" ]]; then
    # read mapping
    bwa mem -t "$CPU" -R "@RG\\tID:$id\\tSM:$id\\tPL:Illumina" -M "$ref" "$f1" "$f2" \
    | samtools sort -n -@ "$CPU" - \
    | samtools fixmate -u -m - - \
    | samtools sort -u -@ "$CPU" -T "$id" - \
    | samtools markdup -r -@ "$CPU" - "$fmdk"

    # index bam file
    samtools index "$fmdk"

    # calculate mapping statistics
    samtools flagstat "$fmdk" > "$ffst"

    # calculate coverage statistics
    samtools coverage "$fmdk" > "$fcov"
    
    # mapped read depth
    bedtools genomecov -bga -split -ibam "$fmdk" > "$fd"
    awk '$4 > 9'  "$fd" > $ftmp; bedtools merge -i $ftmp > "$fd10"
    awk '$4 > 19' "$fd" > $ftmp; bedtools merge -i $ftmp > "$fd20"
    gzip -f "$fd" "$fd10" "$fd20"

    # clean up
    rm "$ftmp"
  fi
}


#######################################
# Genotype calling
# Arguments:
#   wd refname dref pop fsample fploidy minDP maxDP minGQ
# Returns:
#   vcf files
#######################################
genotyping_bcftools() {
  wd="$1"
  refname="$2"
  ref="$3"
  pop="$4"
  sample_info="$5"  # for sex info
  fploidy="$6"
  
  minDP="$7"
  maxDP="$8"
  minGQ="$9"

  # directories
  local din="$wd/out/bam"
  local dout="$wd/out/vcf"
  local dout_pop="$dout/$pop"

  # ploidy file
  local fsample="$dout_pop/sample.txt"

  # output files
  local fvcf="$dout_pop/$pop.vcf.gz"
  local ffilt="$dout_pop/$pop.filt.minDP$minDP.maxDP$maxDP.minGQ$minGQ.vcf.gz"
  local fstat="$dout_pop/$pop.filt.minDP$minDP.maxDP$maxDP.minGQ$minGQ.txt"

  # half depth for haploid genotypes
  local minDP2
  minDP2="$(( $minDP / 2 ))"

  # vcf filter
  local vcf_filter="( FMT/DP >= ${minDP} & FMT/DP <= ${maxDP} & FMT/GQ >= ${minGQ} ) | ( FMT/GT=='"'0'"' & FMT/DP >= ${minDP2} ) | ( FMT/GT=='"'1'"' & FMT/DP >= ${minDP2} ) | ( FMT/GT=='"'2'"' & FMT/DP >= ${minDP2} )"


  if [[ ! -f "$ffilt" ]]; then
    # individuals from this pop
    declare -a ids=( $( grep "$pop" "$sample_info" | cut -d' ' -f 1 ) )

    mkdir -p "$dout_pop"

    # input bam files
    declare -a ff=()
    for id in "${ids[@]}"; do
      ff+=( "$din/$id.markdup.bam" )
    done

    # clear content if any
    echo -n "" > "$fsample"

    for id in "${ids[@]}"; do
      grep "$id" "$sample_info" | cut -d' ' -f -3 >> "$fsample"
    done

    # genotyping using bcftools mpileup & bcftools call
    bcftools mpileup -B -Ou -f "$ref" "${ff[@]}" --annotate INFO/AD,FORMAT/AD,FORMAT/DP \
      | bcftools call --annotate GQ -vmO z --ploidy-file "$fploidy" -S "$fsample" -o "$fvcf"

    # filter
    bcftools filter --include "type='snp'" "$fvcf" \
      | bcftools filter --include "$vcf_filter" --SnpGap 5 --set-GTs . -O z -o "$ffilt"

    # index vcf
    tabix -p vcf "$fvcf"
    tabix -p vcf "$ffilt"

    # optional: vcf stats
    bcftools stats -F "$ref" "$ffilt" > "$fstat"
  fi
}


#######################################
# Prepare reference: select loci to use
# Arguments:
#   wd refname dref pop fsample fploidy minDP maxDP minGQ
# Returns:
#   vcf files
#######################################
prepare_ref() {
  wd="$1"
  refname="$2"
  dref="$3"
  region_list="$4"

  # inputs
  local ref="$dref/$refname.fa"
  local gff="$dref/$refname.gff"
  local frepeat="$dref/$refname.repeat.bed"

  # outputs
  local chrom_size="$dref/$refname.chrom.sizes"
  local fcoding="$dref/$refname.coding.bed"
  local fgene="$dref/$refname.gene.bed"
  local ftranscript="$dref/$refname.transcript.bed"
  local fintron="$dref/$refname.intron.bed"
  local fnoncod_unmask="$dref/$refname.noncod.unmask.bed"
  local fnoncod="$dref/$refname.noncod.bed"

  # 1) scaffold sizes and genome annotations
  samtools faidx "$ref"
  cut -f1-2 "$ref.fai" > "$chrom_size"

  # 2) get coordinates of CDS
  grep $'CDS\t' "$gff" \
    | cut -f1,4,5 \
    | sed 's, ,\t,g' \
    | sort -k1,1 -k2,2n \
    | uniq > "$fcoding"

  grep $'gene\t' "$gff" \
    | cut -f1,4,5 \
    | sed 's, ,\t,g' \
    | sort -k1,1 -k2,2n \
    | uniq > "$fgene"

  grep $'transcript\t' "$gff" \
    | cut -f1,4,5 \
    | sed 's, ,\t,g' \
    | sort -k1,1 -k2,2n \
    | uniq > "$ftranscript"

  # 3) get coords for intron [genic non-CDS] & noncoding
  bedtools subtract -a "$ftranscript" -b "$fcoding" > "$fintron"
  bedtools complement -i "$ftranscript" -g "$chrom_size" > "$fnoncod_unmask"
  bedtools subtract -a "$fnoncod_unmask" -b "$frepeat" > "$fnoncod"

  # 4) select loci
  local region minlen maxlen gaplen

  for region_full in "${region_list[@]}"; do
    IFS="." read -r refname region minlen maxlen gaplen <<< "$region_full"

    Rscript "$dscript"/DefineLocusToSelect.R \
      -d "$dref/" \
      -i "$refname.$region.bed" \
      -l "$minlen" -m "$maxlen" -g "$gaplen" \
      -o "$refname.$region.$minlen.$maxlen.$gaplen.bed"
  done
}


#######################################
# Read mapping
# Globals:
#   CPU
# Arguments:
#   wd refname dref pop fbed minDP maxDP minGQ
# Returns:
#   loci in vcf format
#######################################
extract_loci() {
  wd="$1"
  refname="$2"
  dref="$3"
  pop="$4"
  fbed="$5"  # locus coord file
  minDP="$6"
  maxDP="$7"
  minGQ="$8"

  local dvcf="$wd/out/vcf/$pop"
  local fin="$dvcf/$pop.filt.minDP$minDP.maxDP$maxDP.minGQ$minGQ.vcf.gz"

  local doutbase="$wd/out/loci/${fbed%.bed}/$pop"
  local dout="$doutbase"/1.vcf/minDP"$minDP"

  mkdir -p "$dout"

  # extract a vcf file for each locus
  cat "$dref/$fbed" \
    | xargs -n 3 -P "$CPU" sh -c 'bcftools view -r $0:$1-$2 '"$fin -Oz -o $dout"'/${0}_${1}_${2}.vcf.gz'
}


#######################################
# Prepare depth mask for bcftools consensus
# Arguments:
#   wd refname dref d
# Returns:
#   depth mask
#######################################
depth_mask() {
  wd="$1"
  refname="$2"
  ref="$3"
  chrom_size="$4"
  d="$5"  # depth cutoff

  local ddepth="$wd/out/depth"
  local dout="$wd/out/mask"
  local fout

  mkdir -p "$dout"
  cd "$ddepth"

  for f in *.d$d.txt.gz; do
    fout="${f//.txt.gz/.mask.bed}"
    
    if [[ ! -f "$fout.gz" ]]; then
      echo "$fout"

      bedtools complement -i "$f" -g "$chrom_size" > "$fout"
      gzip "$fout"
      mv "$fout.gz" "${dout}/"
    fi
  done
}


#######################################
# Generate sequence from genotype data
# Arguments:
#   wd refname dref d region_full
# Returns:
#   loci in fasta format
#######################################
vcf_to_fasta() {
  wd="$1"
  refname="$2"
  ref="$3"
  pop="$4"
  d="$5"
  region_full="$6"

  # half depth for haploid genotypes
  local d2="$(( $d / 2 ))"

  local din="$wd/out/loci/$region_full/$pop"
  local dmask="$wd/out/mask"
  local dvcf="$din/1.vcf/minDP$d"
  local doutbase="$din/2.fasta"
  local dout="$doutbase/minDP$d"
  local fname ch start end fin fout locus ftmp sample_list sample gt fmask fout_sample

  mkdir -p "$dout"

  # loop over loci
  for f in "$dvcf"/*.vcf.gz; do
    fname="$(basename ${f%.vcf.gz})"
    ch="$(echo $fname | cut -d'_' -f 1)"
    start="$(echo $fname | cut -d'_' -f 2)"
    end="$(echo $fname | cut -d'_' -f 3)"
    fin="$dvcf/${ch}_${start}_${end}.vcf.gz"
    fout="$dout/${ch}_${start}_${end}.fa"

    # comment out if overwrite
    if [[ ! -s "$fout" ]]; then
      echo "${ch}_${start}_${end}"

      # index vcf
      tabix -f "$fin"

      # locus coordinates
      locus="$ch:$start-$end"

      # store screen output from bcftools consensus
      ftmp="$dout/tmp_${ch}_${start}_${end}.txt"

      # get a list of samples in the vcf file
      sample_list=($(bcftools query -l "$fin"))
    
      # create an emppty output file
      echo -n "" > "$fout"

      # loop over samples in the input vcf
      for sample in "${sample_list[@]}"; do
        # get ploidy info from genotype call
        gt="$( bcftools query -f'[%GT\n]' -s $sample $fin | head -1 )"

        # low-depth regions to mask
        if [[ "$gt" != *"/"* && "$gt" != *"|"* ]]; then
          # half cutoff if haploid
          fmask="$dmask/$sample.d$d2.mask.bed.gz"
        else
          fmask="$dmask/$sample.d$d.mask.bed.gz"
        fi

        # get sequence from ref
        fout_sample="$dout/${ch}_${start}_${end}_${sample}.fa"

        samtools faidx "$ref" "$locus" | \
          bcftools consensus "$fin" \
          -s "$sample" \
          -o "$fout_sample" \
          -I \
          -M 'N' \
          -m "$fmask" --mask-with '-' 2> "$ftmp"

        # rename sequence
        sed "s/>.*/>${ch}_${start}_${end}^$sample/" "$fout_sample" >> "$fout"
        
        # clean up
        rm "$fout_sample"
      done

      # clean up
      rm "$ftmp" "$fin.tbi"
    fi

    # clean up
    if [[ -f "$fin.tbi" ]]; then 
      rm "$fin.tbi"
    fi
  done
}


#######################################
# Combine loci into one file
# Arguments:
#   wd mindp region_full
# Returns:
#   Files containing multiple loci
#######################################
combine_fasta() {
  wd="$1"
  mindp="minDP$2"
  region_full="$3"

  local locus_info="$wd/data/ref/$region_full.bed"
  local dinbase="$wd/out/loci/$region_full"
  local dout="$wd/out/dset/$region_full/1.fasta/$mindp"

  mkdir -p "$dout"

  # loop over loci
  while IFS=$'\t' read -r scaffold begin end ; do
    locus="${scaffold}_${begin}_${end}"
    fout="$dout/$locus.fa"
    cat "$dinbase"/*/2.fasta/"$mindp/$locus.fa" > "$fout"
  done < "$locus_info"
}


#######################################
# Convert fasta to phylip
# Arguments:
#   wd mindp region_full
# Returns:
#   Files in phylip format
#######################################
fasta_to_phylip() {
  wd="$1"
  mindp="minDP$2"
  region_full="$3"

  local dinbase="$wd/out/dset/$region_full"
  local din="$dinbase/1.fasta/$mindp"
  local dout="$dinbase/2.phylip/$mindp"
  local fname fout num_sq name sq hdr

  mkdir -p "$dout"

  # loop over loci
  for f in "$din"/*.fa; do
    fname=$(basename ${f%.fa})
    fout="$dout/$fname.txt"

    # create output file
    echo -n "" > "$fout"

    # create header line
    num_sq=`grep -c '>' $f`

    name=""
    sq=""
    hdr=true

    # go line by line
    while read -r line; do
      if [[ "$line" == ">"* ]]; then
        # fasta header line
        if [[ "$name" != "" ]] && [[ "$sq" != "" ]]; then
          # phylip header line (first line)
          if "$hdr"; then
            num_site="${#sq}"
            echo "$num_sq $num_site" >> "$fout"
            hdr=false
          fi

          echo "$name  $sq" >> "$fout"
        fi

        name="${line/>/}"
        sq=""

      else
        # sequence
        sq+="$line"
      fi
    done < "$f"

    # last line
    echo "$name  $sq" >> "$fout"
  done
}


#######################################
# Move final outputs to one location
# Arguments:
#   wd mindp region dout aln_filt_flag makeblock fblock
# Returns:
#   None
#######################################
transfer_output() {
  wd="$1"
  mindp="minDP$2"
  region_full="$3"
  dout="$4"
  aln_filt_flag="$5"
  makeblock="$6"
  fblock="$7"

  local dsource="$wd/out/dset/${region_full}/4.multilocus_${aln_filt_flag}/$mindp"
  local dtarget="$dout/multilocus/$mindp"

  mkdir -p "$dtarget"
  
  find "$dsource" -maxdepth 1 -type f -name '*.txt' -exec mv -i {} "${dtarget}/" \;

  if "$makeblock" && [[ -f "$fblock" ]]; then
    # list of block sizes
    blocksizes=( $(awk '{print $1}' "$fblock" ) )

    for b in "${blocksizes[@]}"; do
      dsource="$wd/out/dset/${region_full}/4.multilocus_${aln_filt_flag}_block${b}/$mindp"
      dtarget="$dout/multilocus_block${b}/$mindp"
      
      mkdir -p "$dtarget"

      # loop over regions
      for dir in "$dsource"/*/; do
        d=$(basename ${dir%*/})
        mv "${dsource}/$d" "${dtarget}/"
      done
    done
  fi
}


#######################################
# Main pipeline
#######################################
main() {
  timer=$SECONDS

  # log screen outputs
  flog="$wd/log.txt"
  echo -n "" > "$flog"

  # list of regions of loci to use
  # syntax: refname.region.minlen.maxlen.gaplen
  declare -a region_list=(
    "$refname.coding.100.1000000.2000"
    "$refname.intron.100.1000.2000"
    "$refname.noncod.100.1000.2000"
  )


  ## read mapping ##
  echo "read mapping"
  while IFS=" " read -r sample_id sex pop; do
    read_mapping_pe "$wd" "$refname" "$ref" "$ddata" "$sample_id" >> "$flog" 2>&1
  done < "$fsample"


  ## genotyping ##
  echo "genotyping"
  while IFS=" " read -r pop; do
    genotyping_bcftools "$wd" "$refname" "$ref" "$pop" "$fsample" "$fploidy" "$minDP" "$maxDP" "$minGQ" >> "$flog" 2>&1
  done < "$fpop"


  echo "extracting loci"

  ## optional: prepare reference in dref ##
  # prepare_ref "$wd" "$refname" "$dref" "$region_list" >> "$flog" 2>&1
  

  ## extract individual loci (vcf) ##
  echo "extracting vcf loci"
  while IFS=" " read -r pop; do
    for region_full in "${region_list[@]}"; do
      fbed="${region_full}.bed"
      extract_loci "$wd" "$refname" "$dref" "$pop" "$fbed" "$minDP" "$maxDP" "$minGQ" >> "$flog" 2>&1
    done
  done < "$fpop"


  ## create depth masks ##
  echo "creating depth masks"
  
  # require refname.chrom.sizes from prepare_ref
  chrom_size="$dref/$refname.chrom.sizes"
  
  # half depth cutoff for haploid genotypes
  minDP2="$(( $minDP / 2 ))"

  depth_mask "$wd" "$refname" "$ref" "$chrom_size" "$minDP"  >> "$flog" 2>&1
  depth_mask "$wd" "$refname" "$ref" "$chrom_size" "$minDP2" >> "$flog" 2>&1


  ## convert vcf to fasta ##
  echo "converting vcf to fasta"
  ftmp="tmp.txt"
  while IFS=" " read -r pop; do
    for region_full in "${region_list[@]}"; do
      echo " $pop $region_full"
      vcf_to_fasta "$wd" "$refname" "$ref" "$pop" "$minDP" "$region_full" >> "$ftmp" 2>&1
    done
  done < "$fpop"
  rm "$ftmp"


  ## optional: processing of fasta ##
  echo "generating datasets"
  for region_full in "${region_list[@]}"; do
    ## collect fasta sequences into one file per locus ##
    combine_fasta "$wd" "$minDP" "$region_full"

    ## convert fasta to phylip ##
    fasta_to_phylip "$wd" "$minDP" "$region_full"

    ## filter phylip ##
    python $dscript/filter-phylip-loci.py "$wd" "$refname" "$minDP" "$region_full" "$dref" "$missing_prop_cutoff_row" "$missing_prop_cutoff_col" "$aln_filt_flag" >> "$flog" 2>&1


    ## make blocks ##
    python $dscript/make-blocks.py "$wd" "$refname" "$minDP" "$region_full" "$dref" "$aln_filt_flag" "$makeblock" "$fblock" >> "$flog" 2>&1
    
    ## move output files to final location  ##
    transfer_output "$wd" "$minDP" "$region_full" "$dout" "$aln_filt_flag" "$makeblock" "$fblock"
  done

  echo "We're done!"

  timer=$(($SECONDS-timer))
  printf "Time used: %02d:%02d:%02d\n" "$((timer/3600))" "$((timer/60%60))" "$((timer%60))"
}

main "$@"
