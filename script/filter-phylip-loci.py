#!/usr/local/bin/python3

'''
Loop over loci
- remove sequences with missing data > missing_prop_cutoff_row
- remove columns with all missing data
- discard loci if too few sequences left after above filtering
'''

import sys
import re
import os
import math


def filter_loci(dir_in, dir_out, mindp, f_loci, wd, missing_prop_cutoff_row, missing_prop_cutoff_col, min_num_sqs, verbose):

  f = open(f_loci, 'r')
  locus_coord = f.readlines()
  f.close()

  # out = region
  out = os.path.basename(f_loci).replace('.bed', '')
  print(f_loci, 'minDP' + str(mindp), out, dir_in, '->', dir_out)

  din = os.path.join(wd, out, dir_in, 'minDP' + str(mindp))
  if not os.path.exists(din): sys.exit("Input dir doesn't exist: {}".format(din))

  dout = os.path.join(wd, out, dir_out, 'minDP' + str(mindp))
  if not os.path.exists(dout): os.makedirs(dout)

  # loop over list of loci in locus_coord
  for i_locus in range(len(locus_coord)):
    chrom, coord_from, coord_to = locus_coord[i_locus].rstrip().split('\t')
    
    # convert to int
    # use float() first to avoid 'ValueError: invalid literal for int() with base 10'
    coord_from = int(float(coord_from))
    coord_to = int(float(coord_to))

    if verbose: print(chrom, coord_from, coord_to)

    fname = '{}_{}_{}.txt'.format(chrom, str(coord_from), str(coord_to))

    fin = os.path.join(din, fname)
    fout = os.path.join(dout, fname)

    if os.path.isfile(fin):
      # sequence data
      f = open(fin, 'r')
      lines = [line for line in f.readlines() if line.strip()]
      f.close()

      if not lines:
        print('data error: ' + fname)
        continue

      # header line
      hdr = lines[0].rstrip()
      num_sq, sq_len = hdr.split(' ')
      num_sq = int(num_sq)
      sq_len = int(sq_len)
      if verbose: print(hdr)

      # remaining lines
      lines = [lines[i].rstrip() for i in range(1, len(lines))]
      
      # store output
      sqs = {}
      missing = {}  # counts of '-' and 'N'
      for line in lines:
        try:
          name, sq = line.split('  ')
          num_missing = sq.count('-') + sq.count('N')
          
          if verbose: print(name, str(num_missing), 'EXCLUDE' if num_missing >= missing_prop_cutoff_row * sq_len else '')

          # only include sequences if proportion of missing nts are below the cutoff
          if num_missing < missing_prop_cutoff_row * sq_len:
            sqs[name] = sq
            missing[name] = num_missing

        except:
          print('data error: ' + fname)

      # number of remaining sequences
      num_sq_new = len(sqs)

      # skip if num_sq_new < min_num_sqs
      if num_sq_new >= min_num_sqs:

        # remove columns with missing data > missing_prop_cutoff_col
        cols_remove = []
        for i_col in range(sq_len):
          # print(i_col, ''.join([sqs[name][i_col] for name in sqs]))
          sq_col = ''.join([sqs[name][i_col] for name in sqs])
          num_missing_col = sq_col.count('-') + sq_col.count('N')

          if (num_missing_col > 0) and (num_missing_col / num_sq_new >= missing_prop_cutoff_col):
            cols_remove.append(i_col)
        # print(len(cols_remove))

        sq_len_new = sq_len - len(cols_remove)
        hdr_new = str(num_sq_new) + ' ' + str(sq_len_new)
        if verbose: print(hdr_new)

        # write to output file
        if sq_len_new > 0:
          # columns to keep
          cols_keep = [k for k in range(sq_len) if k not in cols_remove]

          with open(fout, 'w') as f:
            f.write(hdr_new + '\n')

            for name in sqs:
              f.write(name + '  ' + ''.join(sqs[name][k] for k in cols_keep) + '\n')

      if verbose: print('\n')


def main():

  wd = sys.argv[1]
  refname = sys.argv[2]
  mindp = sys.argv[3]
  region_full = sys.argv[4]
  dref = sys.argv[5]  # locus data

  # row filter
  missing_prop_cutoff_row = 0.5

  # column filter
  missing_prop_cutoff_col = 0.5
  flag = 'f50'

  # minimun number of sequences per locus
  min_num_sqs = 2

  verbose = False
  # verbose = True

  dout = wd + '/out/dset'
  
  dir_in = '2.phylip'
  dir_out = '3.phylip_filtered_' + flag

  f_loci = os.path.join(dref, '{}.bed'.format(region_full))

  filter_loci(dir_in, dir_out, mindp, f_loci, dout, missing_prop_cutoff_row, missing_prop_cutoff_col, min_num_sqs, verbose)


if __name__ == '__main__':
  main()
