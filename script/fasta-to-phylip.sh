#!/bin/bash

wd=$1
mindp=minDP$2
region_full=$3

dinbase=$wd/out/dset/$region_full
din=$dinbase/1.fasta/$mindp
dout=$dinbase/2.phylip/$mindp

mkdir -p $dout

# loop over loci
for f in $din/*.fa; do
  fname=$(basename ${f%.fa})
  fout=$dout/$fname.txt

  # create output file
  echo -n "" > $fout

  # create header line
  num_sq=$( grep '>' $f | wc -l )

  name=""
  sq=""
  hdr=true

  # go line by line
  while read -r line; do
    if [[ $line == ">"* ]]; then
      # fasta header line
      if [[ "$name" != "" ]] && [[ "$sq" != "" ]]; then
        # phylip header line (first line)
        if $hdr; then
          num_site=${#sq}
          echo $num_sq $num_site >> $fout
          hdr=false
        fi

        echo "$name  $sq" >> $fout
      fi

      name=${line/>/}
      sq=""

    else
      # sequence
      sq+="$line"
    fi
  done < "$f"

  # last line
  echo "$name  $sq" >> $fout
done
