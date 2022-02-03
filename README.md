# glopp (name to be decided) : polyploid phasing from read sequencing

## Introduction

**glopp** is a software package for single individual haplotype phasing of polyploid organisms from read sequencing. 

Given 

1. a list of variants in .vcf format
2. a set of reads mapped to a reference in .bam format

**glopp** outputs a set of phased haplotypes.

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
2. [cmake](https://cmake.org/download/) version > 3.12 is required. It's sufficient to download the binary from the link and do `PATH="/path/to/cmake-3.xx.x-linux-x86_64/bin/:$PATH"` before installation. 
### Install

```
git clone https://github.com/bluenote-1577/glopp
cd glopp
git checkout flow
cargo build --release
./target/release/glopp -h
```

`cargo build --release` builds the **glopp** binary, which is found in the ./target/release/ directory. 

## Using glopp

```
glopp -b bamfile.bam -c vcffile.vcf -o output_dir (Estimate ploidy using heuristic)
```
For a quick test, we provide a VCF and BAM files in the tests folder. Run ``./target/release/glopp -b tests/test_bams/pds_ploidy3.bam -c tests/test_vcfs/pds.vcf -o results`` to run glopp on a 3 Mb section of a simulated 3x ploidy potato chromosome with 30x read coverage.

### BAM + VCF
The standard mode of usage is to specify a bam file using the option **-b** and a vcf file using the option **-c**. The output is written to a text file with value of option **-o**. 

glopp currently only uses SNP information and does not take into account indels. 

The bam file may contain multiple contigs/references which the reads are mapped to as long as the corresponding contigs also appear in the vcf file.

## Output

glopp outputs the sequence of SNPs (i.e. the phasing) and the partition of reads associated to each haplotype. 

The results are found in the output directory. This directory is specified by the `-o` option, or `glopp_out_dir` by default.

### Phased haplotype output
For each contig, glopp outputs a phased haplotype file `(contig_name)_phasing.txt` in the following format:

1. Column 1 is (variant) : (genome position) where (variant) is the i-th variant, and the genome position is the the position of the genome on the reference.
2. The next k columns are the k phased haplotypes for an organism of ploidy k. 0 represents the reference allele, 1 the first alternate, and so forth. 
3. The next k columns are of the form (allele):(# calls)|(allele):(# calls) where (allele) = 0,1,... and (# calls) is the number of reads assigned to a specific haplotype which contain that allele. For example, 0:10|1:5 indicates that 10 reads assigned to this haplotype have allele 0 at this position, and 5 reads have allele 1. 

If using a bam file with multiple contigs being mapped to, the output file contains multiple phased haplotypes of the above format which are delimited by `**(contig name)**`.

Glopp will automatically detect regions that can not be phased across unambiguously. These regions are delimited by `--------` in this file. 

### Read partition output 
For each contig, glopp outputs a partition of reads in the file `(contig_name)_partition.txt` in the following format:

```
#1 (partition #1)
(read_name1) (first SNP position covered) 
(read_name2) (first SNP position covered)
...
#2 (partition #2)
...
```
## Misc.

### VCF requires contig headers
We found that some variant callers don't put contig headers in the VCF file. In this situation, run `python scripts/write_contig_headers_vcf.py (vcf_file)` to get a new VCF with contig headers.

### Output BAM partition
To get a set of BAM files which correspond to the output read partition (i.e. the haplotypes), use

``python scripts/get_bam_partition.py (-P output file) (original BAM file) (prefix name for output)``

This will output a set of bams labelled `prefix_name1.bam`, `prefix_name2.bam` and so forth. This script requires pysam.

### Manually consensus for testing

Suppose you already have a partitioning of reads. That is, you have bam files `bam_file1, bam_file2, bam_file3` and you want to use this partioning for the phasing. Use the `consensus` binary to get a phasing from the .bam files by `consensus -v (vcf_file) -b (bam_file1) (bam_file2) (bam_file3) -o (consensus_file.txt)`. This is useful for visualizing error rates.


