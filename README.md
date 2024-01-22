# floria - metagenomic long or short-read strain haplotype phasing

## Introduction

**floria** is a software package for recovering microbial haplotypes and clustering reads at the strain level from metagenomic sequencing data. 

After calling SNPs against reference genomes/a metagenomic assembly, floria produce strain-level clusters of reads and their haplotypes **in minutes**. 

<p align="center">
  <img width="460" height="400" src="https://github.com/bluenote-1577/vartig-utils/blob/main/visualize-vartig-example.png", caption="asdf">
</p>
<p align="center">
  <i>
    A 1Mbp contig (Brevefilum fermentans) was automatically phased into two strains by floria in minutes. Only two strains are present with high HAPQ; spurious "haplosets" are given low HAPQ.
  </i>
</p>

### Inputs and outputs

Given 

1. a list of variants in .vcf format
2. a set of reads mapped to contigs/references in .bam format

**floria** performs strain/haplotype phasing on short-read or long-read shotgun metagenomic samples. For more an introduction to floria and its inputs/outputs, see the documentation below. 

## Full documentation

See https://phase-doc.readthedocs.io/en/latest/index.html for more information.

## End-to-end phasing pipeline

See the **"Production"** pipeline [here](https://github.com/jsgounot/Floria_analysis_workflow) for reads-to-haplotype pipelines if you do not know how to get started with generating VCFs or BAMs. 

## Install 

Compiling floria from scratch should be relatively simple. Otherwise, a static binary is provided. 

#### Option 1 - compile from scratch

A relatively recent standard toolchain is needed.

1. [rust](https://www.rust-lang.org/tools/install) **version > 1.63.0** and associated tools such as cargo are required and assumed to be in PATH.
2. [cmake](https://cmake.org/download/) **version > 3.12** is required. It's sufficient to download the binary from the link and do `PATH="/path/to/cmake-3.xx.x-linux-x86_64/bin/:$PATH"` before installation. 
3. make 
4. GCC 

If you're using an **x86-64 architecture with SSE instructions (most linux systems)**: 

```sh
git clone https://github.com/bluenote-1577/floria
cd floria

cargo install --path . 
floria -h # binary is available in PATH

# OR IF ~/.cargo is unavailable for some reason

#cargo build --release
#./target/release/floria -h # binary built in ./target/release instead.
```

If you're using an **ARM architecture with NEON instructions** (e.g. Mac M1): 

```sh

# If using ARM architecture with NEON instructions
cargo install --path . --root ~/.cargo --features=neon --no-default-features
floria -h # binary is available in PATH

```

#### Option 2 - precompiled static binary on **x86-64-linux**

The static binary only for x86-64 linux with SSE instructions currently. 

```sh
wget https://github.com/bluenote-1577/floria/releases/download/latest/floria
chmod +x floria
./floria -h
```
