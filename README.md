# glopp (name to be decided) 

## Introduction

**glopp** is a software package for single individual haplotype phasing of polyploid organisms from read sequencing. 

Given 

1. a list of variants in .vcf format
2. a set of reads mapped to contigs/references in .bam format (**sorted and indexed**)

**glopp** performs strain/haplotype phasing.

### Requirements 

A relatively recent toolchain is needed, but no other dependencies. 

1. [rust](https://www.rust-lang.org/tools/install) **version > 1.63.0** and associated tools such as cargo are required and assumed to be in PATH.
2. [cmake](https://cmake.org/download/) **version > 3.12** is required. It's sufficient to download the binary from the link and do `PATH="/path/to/cmake-3.xx.x-linux-x86_64/bin/:$PATH"` before installation. 
3. make 
4. GCC 

Alternatively, we offer a statically compiled binary for x86-64 linux if the above requirements can't be met. 

### Install

#### Option 1 - compile from scratch

If you're using an x86-64 architecture with avx2 instructions: 

```sh
git clone https://github.com/bluenote-1577/glopp
cd glopp

# Option 1) AVX2 instructions are available
cargo install --path . --root ~/.cargo 
glopp -h # binary is available in PATH

# OR IF ~/.cargo is unavailable for some reason

cargo build --release
./target/release/glopp -h # binary built in ./target/release instead.
```
If you're using an ARM architecture with NEON instructions (e.g. Mac M1): 

```sh

# Option 2) If using ARM architecture with NEON instructions
cargo install --path . --root ~/.cargo --features=neon --no-default-features
glopp -h # binary is available in PATH

```

#### Option 2 - precompiled static binary on **x86-64-linux**

If you're using ARM then you'll have to compile. See above. 

```sh
wget https://github.com/bluenote-1577/glopp/releases/download/latest/glopp
chmod +x glopp
./glopp -h
```

## Using glopp

``` sh

# Options for glopp
glopp -h

# Basic minimal run -- automatic estimation of parameters, 10 threads by default
glopp -b bamfile.bam -v vcffile.vcf -r references.fa

# Phase only contigs listed in the -G option to output_dir, 30 threads
glopp -b bamfile.bam -v vcffile.vcf -r references.fa -o output_dir -G contig_1 contig_2 contig_3 -t 30

# Options to filter bam using MAPQ, outputting reads, filtering contigs. See manual for more information.
glopp -b bamfile.bam -v vcffile.vcf -r references.fa -m 60 --output-reads --snp-count-filter 1000

```
For a quick test, we provide a VCF and BAM files in the tests folder. Run
```
glopp -b tests/test_long.bam -v tests/test.vcf -r tests/MN-03.fa 
```
to run glopp on a 100 kb section of a simulated 3-strain Klebsiella Pneumoniae sample. Look through the resulting `glopp_out_dir` folder to get a sense of glopp's output format. 

### Parameters for best performance

4. **-e** controls how sensitive your blocks are during phasing. Blocks with error rate less than **-e** with not be phased further, so haplotypes that differ less than **-e** may be combined. Use a higher value for more contiguous but less sensitive phasings, and a lower value if you want a more sensitive but broken phasings. If using hybrid correction, maybe try setting this lower to 0.02 (default is 0.04). 

## Output

```
results
│      
└───contig1_in_bam
│   │   all_part.txt
│   │   pet_graph.dot
|   |   (other debug files)
│   │
│   └───haplotypes
│   |   │   0_hap.txt
│   |   │   1_hap.txt
│   |   │   ...
|   |
│   └───long_reads
│   |   │   0_part.fastq
│   |   |   ...
|   |
│   └───short_reads
│   |   │   0_part_paired1.fastq
│   |   |   0_part_paired2.fastq
|   |
|   └───(debug_folders)
└───contig2_in_bam
    │   ...
    │   ...
```
glopp outputs a set of **haplotigs**. We define a haplotig to be a set of reads that belong to the same strain. The collection of all haplotigs is found in the `results/contig1/all_part.txt` file. 

The following information is output for each haplotig:

1. a **haplotype** which is the sequence of SNPs on each haplotig in the `haplotypes` folder. 
2. trimmed long-reads (if using long-reads) corresponding to each haplotig are found in the `long_reads` folder. 
3. trimmed short-reads (if using short-reads) corresponding to each haplotig are found in the `short_reads` folder. 

### Haplotigs ``results/contig/all_part.txt`` 

Each haplotig corresponds to a cluster of reads and is presented in the following format:

```
#0,(coverage for haplotig 0),(error_rate for haplotig 0) 
(read_name1) (first SNP position covered) 
(read_name2) (first SNP position covered)
...
#1,(coverage for haplotig 1),(error_rate for haplotig 1)
...
```

### Haplotype output ``results/contig/haplotypes/``
For each haplotig, glopp outputs a haplotype file `#_hap.txt` in the following format:

```
>(haplotig number),(left snp cutoff position),(right snp cutoff position)
(snp #1):(genome position)     (consensus allele #: 0/1/2...)    (allele #1):(support)|(allele #2):(support)|...
(snp #2):(genome position)     (consensus allele #: 0/1/2...)    (allele #1):(support)|(allele #2):(support)|...
...

```

The haplotig is valid for the region of the genome that lies within the left snp to the right snp cutoff positions. 

1. Col. 1 indicates the # and position of the SNP. 
2. Col. 2 states which allele is the consensus allele.
3. Col. 3 describes how many reads support each allele (i.e. how many reads in the haplotig have the allele at the SNP position). 

### Read output ``results/contig/*_reads/``

The reads in each haplotig can be found in either the `long_reads` or `short_reads` folder, depending on which type of read is used. Note that fastq files in these folders are trimmed to lie within an interval and thus differ from the original reads. This is done so that all reads in a haplotig fall within an interval on the genome and do not extend past the interval. 

### Debugging

Extra debug files in `local_parts` and `debug_paths` show the local partitions and the path corresponding to the haplotigs and local partitions. To visualize the flow-graph constructed, a graphviz `pet_graph.dot` file is included. If graphviz is installed, this can be visualized by running `dot -Tps results/contig/pet_graph.dot -o outfile.ps` and looking at the resulting `outfile.ps`. 

## Assembling output reads in `results/contig/*_reads/`

If you want to assemble the haplotigs (in the same way strainberry does) then the utility scripts `strains_phase_scripts/assemble_from_glopp_out.py` or `assemble_shortreads_from_glopp_out.py` for long and short reads respectively allow you to do so. Ensure that 

1. wtdbg2 (for long-reads)
2. minimap2
3. abyss (for short-reads)

and associated binaries are present in path or edit the files to include paths to the binaries. Run 

```strains_phase_scripts/assemble_from_glopp_out.py results_dir results_dir reference.fa```

where `results_dir` is the output of glopp (i.e. the `-o` output) for long-reads and the other script for short reads. 

All haplotigs will be assembled in `results_dir/intermediate` and all assembled haplotigs will be mapped onto `reference.fa` in the `results_dir/all_assemblies.bam` bam file. 

**NOTE**: this only works for single contig phasings, i.e. your bam file is only aligned to a single reference. I will test this later for multi-contig phasings.


## Extra scripts

### VCF requires contig headers
We found that some variant callers don't put contig headers in the VCF file. In this situation, run `python scripts/write_contig_headers_vcf.py (vcf_file)` to get a new VCF with contig headers.

### Haplotagging bams for visualization

-- requires natsort
-- requires pysam 

To generate a new bam file that is tagged with the `HP:i` tag, we offer a python script called `haplotag_bam.py` in the `scripts` directory. 

`python glopp/scripts/haplotag_bam.py glopp_out_dir/CONTIG/all_part.txt ORIGINAL_BAM_FILE.bam NEW_BAM_FILE.bam CONTIG`

For a single contig CONTIG, this generates a new bam file with `HP:i` tags called NEW_BAM_FILE.bam. This new bam file contains all reads in the ORIGINAL_BAM_FILE.bam mapped to CONTIG. You can then index NEW_BAM_FILE.bam and visualize it with [igv](https://igv.org/).  
