#!/usr/local/bin/python

'''
Convert group-level vcf files from gatk genotypeGVCFs into sequence alignment
in phylip format
 - main input: joint_filtered.vcf.gz
 - positions to exclude:
     1) in joint.vcf.gz but not in joint_filtered.vcf.gz (incl low quality SNPs and non-SNP variants)
     2) mapped read depth < 20
 - filling in the remaining non-SNP positions using ref sq
'''

import sys
import re
import os


def make_multilocus_block(dir_in, dir_out, blocksizes, mindp, f_loci, fout_prefix, wd, last_block_size_cutoff, thin, verbose):

  ## STEP 1: # combine loci into a single file
  out = os.path.basename(f_loci).replace('.bed', '')
  print(out)
  
  # optional: skip some loci
  if thin:
    f_loci = f_loci.replace('.bed', '.' + thin + '.bed')
    dir_out += '_' + thin

  din = os.path.join(wd, out, dir_in, 'minDP' + mindp)
  if not os.path.exists(din): sys.exit("Input dir doesn't exist: {}".format(din))

  dout = os.path.join(wd, out, dir_out, 'minDP' + mindp)
  if not os.path.exists(dout): os.makedirs(dout)
  print(dout)

  f = open(f_loci, 'r')
  locus_coord = f.readlines()
  f.close()

  # store sequence data
  sqs = {}  # list of phylip alignment in text format

  # store chrom label of each locus
  chrom_label = {}

  # loop over list of loci in locus_coord
  for i_locus in range(len(locus_coord)):

    chrom, coord_from, coord_to = locus_coord[i_locus].rstrip().split('\t')

    fname = '{}_{}_{}.txt'.format(chrom, coord_from, coord_to)
    fin = os.path.join(din, fname)

    if os.path.isfile(fin):
      if verbose: print(i_locus, chrom, coord_from, coord_to)

      # record chr of this locus
      if chrom not in chrom_label: 
        chrom_label[chrom] = [i_locus]
      else:
        chrom_label[chrom].append(i_locus)
      
      # read in sequence
      f = open(fin, 'r')
      sqs[i_locus] = f.readlines()
      f.close()

  if verbose: print('chrom_label', list(chrom_label.keys()))

  # write output
  for chrom in chrom_label.keys():
    print(chrom + ':', len(chrom_label[chrom]), 'loci')
    
    fname = fout_prefix + chrom + '.txt'
    fout = os.path.join(dout, fname)
    
    with open(fout, 'w') as f:
      for i_locus in chrom_label[chrom]:
        # if verbose: print(*sqs[i_locus])

        for l in sqs[i_locus]:
          f.write(l)
        f.write('\n')


  ## STEP 2: split into blocks
  for b in blocksizes:
    print('block' + str(b))

    dir_out_block = dir_out + '_block' + str(b)
  
    # optional: skip some loci
    if thin:
      f_loci = f_loci.replace('.bed', '.' + thin + '.bed')
      dir_out += '_' + thin
      dir_out_block += '_' + thin
    
    for chrom in chrom_label:
      fname = fout_prefix + chrom + '.txt'
      fin = os.path.join(dout, fname)

      dout_block = os.path.join(wd, out, dir_out_block, 
        'minDP' + mindp, fout_prefix + chrom)
      if not os.path.exists(dout_block): os.makedirs(dout_block)

      # count #blocks
      num_blocks = len(chrom_label[chrom]) // b + 1
      last_block_size = len(chrom_label[chrom]) % b

      # combine last block with the preceding one
      if num_blocks > 1 and last_block_size > 0 and last_block_size < last_block_size_cutoff[b]:
        num_blocks -= 1
      
      print(' {}: {:3d} block{} (last block: {})'.format(chrom, num_blocks, 's' if num_blocks > 1 else '', last_block_size))

      for k in range(num_blocks):
        first_locus = k * b
        last_locus = (k + 1) * b if k < num_blocks - 1 else len(chrom_label[chrom])
        
        # print(k, first_locus, last_locus)

        fout = os.path.join(dout_block, '{:03d}.txt'.format(k + 1))
        f = open(fout, 'w')

        for i_locus in chrom_label[chrom][first_locus:last_locus]:
          for l in sqs[i_locus]:
            f.write(l)
          f.write('\n')

        f.close()


def main():

  wd = sys.argv[1]
  refname = sys.argv[2]
  mindp = sys.argv[3]
  region_full = sys.argv[4]
  dref = sys.argv[5]  # locus data
  flag = sys.argv[6]

  blocksizes = [100]
  last_block_size_cutoff = {100: 40}

  verbose = False
  # verbose = True

  dout = wd + '/out/dset'

  dir_in  = '3.phylip_filtered_' + flag
  dir_out = '4.multilocus_' + flag
  
  # optional: increase locus spacing by thinning out loci
  thin = ''
    
  f_loci = os.path.join(dref, '{}.bed'.format(region_full))
  fout_prefix = os.path.basename(f_loci).replace('.bed', '') + '_'

  print('\n' 'minDP' + mindp, dir_in, '->', dir_out, thin)
  make_multilocus_block(dir_in, dir_out, blocksizes, mindp, f_loci, fout_prefix, dout, last_block_size_cutoff, thin, verbose)


if __name__ == '__main__':
  main()
