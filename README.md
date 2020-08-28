# flopp : fast local polyploid phasing from long read sequencing.

## Introduction

**flopp** is a software package for single individual haplotype phasing of polyploid organisms from long read sequencing. flopp is a extremely fast due to being multithreaded and written entirely in the rust programming language, offering an order of magnitude speedup and better accuracies compared to current polyploid haplotype phasing algorithms. 

Given 

1. a list of variants in .vcf format
2. a set of reads mapped to a reference in .bam format

for a single individual, **flopp** outputs a set of phased haplotypes.

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
### Install

```
git clone https://github.com/bluenote-1577/flopp
cd flopp
cargo build --release
./target/release/flopp -h
```

`cargo build --release` builds the **flopp** binary, which is found in the ./target/release/ directory. 

## Using flopp

```
flopp -b bamfile.bam -v vcffile.vcf -p (ploidy) -o results.txt (with VCF and BAM)
flopp -f fragfile.frags -p (ploidy) -o unpolished_results.txt (with fragment file)
flopp -f fragfile.frags -v vcffile.vcf -p (ploidy) -o polished_results.txt (polished output with fragment file and VCF)
```
The ploidy of the organism must be specified. The number of threads (default 10) can be specified using the -t option. See `flopp -h` for more information.  

### BAM + VCF
The standard mode of usage is to specify a bam file using the option **-b** and a vcf file using the option **-v**. The output is written to a text file with value of option **-o**. 

### Fragment file
A user can also input a fragment file using the option **-f**. The fragment file is a line separated file of reads indexed by SNP positions; see https://github.com/MinzhuXie/H-PoPG or https://github.com/realabolfazl/AltHap for more details about the fragment file specifcaiton (called the *input snp matrix* by H-PoP). Specifying a compatible vcf file with a fragment file uses genotyping information to produce a higher quality output. 

For testing purposes and compatibility with other haplotype phasing algorithms, the binary **frag-dump** is provided in the same folder as the **flopp** binary. 

`frag-dump -b bamfile.bam -v vcffile.vcf -o frags.txt` gives a fragment file a.k.a input snp matrix which is compatible with H-PoP and other haplotype phasing algorithms. 

### Output

TODO

## Citing flopp

Shaw, J. Yu, Y.W. flopp - fast polyploid phasing by long read sequencing using a local assembly approach. (In preparation)


