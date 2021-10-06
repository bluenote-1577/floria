# glopp : global polyploid phasing from long read sequencing.

## Introduction

**glopp** is a software package for single individual haplotype phasing of polyploid organisms from long read sequencing. 

Given 

1. a list of variants in .vcf format
2. a set of reads mapped to a reference in .bam format

**glopp** outputs a set of phased haplotypes.

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
### Install

```
git clone https://github.com/bluenote-1577/glopp
cd glopp
cargo build --release
./target/release/glopp -h
```

`cargo build --release` builds the **glopp** binary, which is found in the ./target/release/ directory. 

## Using glopp

```
flopp -b bamfile.bam -c vcffile.vcf -p (ploidy) -o results.txt -P output_partition_dir (Fix the ploidy)
flopp -b bamfile.bam -c vcffile.vcf -o results.txt (Estimate ploidy using heuristic)
```
For a quick test, we provide a VCF and BAM files in the tests folder. Run ``glopp -b tests/test_bams/pds_ploidy3.bam -c tests/test_vcfs/pds.vcf -p 3 -o results.txt -P test_results`` to run glopp on a 3 Mb section of a simulated 3x ploidy potato chromosome with 30x read coverage.

### BAM + VCF
The standard mode of usage is to specify a bam file using the option **-b** and a vcf file using the option **-c**. The output is written to a text file with value of option **-o**. 

glopp currently only uses SNP information and does not take into account indels. However, the user may define their own fragments which can be indexed by other types of variants. 

The bam file may contain multiple contigs/references which the reads are mapped to as long as the corresponding contigs also appear in the vcf file.

## Output
### Phased haplotype output (-o option)
glopp outputs a phased haplotype file in the following format:

1. Column 1 is (variant) : (genome position) where (variant) is the i-th variant, and the genome position is the the position of the genome on the reference.
2. The next k columns are the k phased haplotypes for an organism of ploidy k. 0 represents the reference allele, 1 the first alternate, and so forth. 
3. The next k columns are of the form (allele):(# calls)|(allele):(# calls) where (allele) = 0,1,... and (# calls) is the number of reads assigned to a specific haplotype which contain that allele. For example, 0:10|1:5 indicates that 10 reads assigned to this haplotype have allele 0 at this position, and 5 reads have allele 1. 

If using a bam file with multiple contigs being mapped to, the output file contains multiple phased haplotypes of the above format which are delimited by `**(contig name)**`.

### Read partition output (-P option)
If also using `-P` option, glopp outputs the read partition obtained by glopp. That is, set of reads corresponding to each haplotype. The format looks like:
```
#1 (partition #1)
(read_name1) (first SNP position covered) (last SNP position covered)
(read_name2) (first SNP position covered) (last SNP position covered)
...
#2 (partition #2)
...
```

To get a set of BAM files which correspond to the output read partition (i.e. the haplotypes), use

``python scripts/get_bam_partition.py (-P output file) (original BAM file) (prefix name for output)``

This will output a set of bams labelled `prefix_name1.bam`, `prefix_name2.bam` and so forth. This script requires pysam.


