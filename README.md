# multilocus
Generate multilocus sequence data from short-read whole genome sequencing data.

# inputs
- paired-end short reads: `<sample>_{1,2}.fastq.gz`
- reference genome
  * sequence: `<ref>.fasta`
  * CDS annotation: `<ref>.gff`
  * repeat coordinates: `<ref>.repeat.bed`
  * [ploidy file](https://samtools.github.io/bcftools/bcftools.html#ploidy): `<ref>.ploidy.txt`

# outputs
Multilocus sequence data in phylip format for [bpp](https://github.com/bpp/bpp) analysis.

# requirements
- [`bwa`](https://github.com/lh3/bwa)
- [`samtools`, `bcftools`](https://www.htslib.org/download/)
- [`bedtools`](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- [`gatk`](https://github.com/broadinstitute/gatk) for `MarkDuplicates`
optional: `R` with package `optparse`
