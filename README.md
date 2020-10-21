# flopp : fast local polyploid phasing from long read sequencing.

## Introduction

**flopp** is a software package for single individual haplotype phasing of polyploid organisms from long read sequencing. flopp is extremely fast, multithreaded, and written entirely in the rust programming language. flopp offers an order of magnitude speedup and better accuracies compared to current polyploid haplotype phasing algorithms. 

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
A user can also input a fragment file using the option **-f**. The fragment file is a file where each line is a read which is indexed by SNP positions; see https://github.com/MinzhuXie/H-PoPG or https://github.com/realabolfazl/AltHap for more details about the fragment file specifcaiton (called the *input snp matrix* by H-PoP). Specifying a compatible vcf file with a fragment file uses genotyping information to produce a higher quality output. 

For testing purposes and compatibility with other haplotype phasing algorithms, the binary **frag-dump** is provided in the same folder as the **flopp** binary. 

`frag-dump -b bamfile.bam -v vcffile.vcf -o frags.txt` gives a fragment file a.k.a input snp matrix which is compatible with H-PoP and other haplotype phasing algorithms. 

### Output

flopp outputs a phased haplotype file in the following format:

1. Column 1 is (variant) : (genome position) where (variant) is the i-th variant, and the genome position is the the position of the genome on the reference.
2. The next k columns are the k phased haplotypes for an organism of ploidy k. 0 represents the reference allele, 1 the first alternate, and so forth. 
3. The next k columns are of the form (allele):(# calls)|(allele):(# calls) where (allele) = 0,1,... and (# calls) is the number of reads in cluster for a specific haplotype which contain that allele. For example, 0:10|1:5 indicates that 10 reads assigned to this haplotype have allele 0 at this position, and 5 reads have allele 1.

## Citing flopp

Jim Shaw and Yun William Yu; Practical probabilistic and graphical formulations of long-read polyploid haplotype phasing. (In preparation)


