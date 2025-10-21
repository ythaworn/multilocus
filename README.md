# multilocus
Generate multilocus sequence data from short-read sequencing data.

# Inputs
- paired-end short reads: `<sample>_{1,2}.fastq.gz`
- sample metadata
  * sample list: `sample.txt`, each line is `<sample id> <sex> <population label>`; `<sex>` should match with ploidy file (below)
  * population list: `pop.txt`, list of population labels, one per line
- reference genome
  * sequence: `<ref>.fasta`
  * CDS annotation: `<ref>.gff`
  * repeat coordinates: `<ref>.repeat.bed`
  * [ploidy file](https://samtools.github.io/bcftools/bcftools.html#ploidy): `<ref>.ploidy.txt`

Optional
- info for grouping loci into blocks: `block_info.txt>`, two columns: block size, minimum size for the last block to be included

# Outputs
Multilocus sequence data in phylip format suitable for [bpp](https://github.com/bpp/bpp) analysis.

# Requirements
- [`bwa`](https://github.com/lh3/bwa)
- [`samtools`, `bcftools`, `tabix`](https://www.htslib.org/download/)
- [`bedtools`](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- `python`
- optional: `R` with package `optparse`

# Usage
See [tutorial](https://github.com/ythaworn/multilocus/wiki/Tutorial)
