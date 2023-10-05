# floria

## Introduction

**floria** is a software package for phasing metagenomic communities at the strain level.

Given 

1. a list of variants in .vcf format
2. a set of reads mapped to contigs/references in .bam format

**floria** performs strain/haplotype phasing on short-read or long-read whole genome metagenomic samples. For more an introduction to floria and its inputs/outputs, see the documentation below. 

## Full documentation

See https://phase-doc.readthedocs.io/en/latest/index.html for more information.

## Install 

Compiling floria from scratch should be relatively simple. Otherwise, a static binary is provided. Conda install is forthcoming.

#### Option 1 - compile from scratch

A relatively recent standard toolchain is needed.

1. [rust](https://www.rust-lang.org/tools/install) **version > 1.63.0** and associated tools such as cargo are required and assumed to be in PATH.
2. [cmake](https://cmake.org/download/) **version > 3.12** is required. It's sufficient to download the binary from the link and do `PATH="/path/to/cmake-3.xx.x-linux-x86_64/bin/:$PATH"` before installation. 
3. make 
4. GCC 

If you're using an **x86-64 architecture with avx2 instructions (e.g. most linux systems)**: 

```sh
git clone https://github.com/bluenote-1577/floria
cd floria

cargo install --path . --root ~/.cargo 
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
If you don't have AVX2 but have **SSE2 only instead (e.g. older systems)** 

```sh

# If using AVX2 not available
cargo install --path . --root ~/.cargo --features=sse2 --no-default-features
floria -h # binary is available in PATH
```

#### Option 2 - precompiled static binary on **x86-64-linux**

The static binary only for x86-64 linux with AVX2 instructions currently. 

```sh
wget https://github.com/bluenote-1577/floria/releases/download/latest/floria
chmod +x floria
./floria -h
```
