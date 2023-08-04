import re
from pyfaidx import Fasta, FastaRecord
import pysam
import argparse


parser = argparse.ArgumentParser(description='Generate new contigs from vartigs by swapping SNPs. ')

parser.add_argument("-f", "--fasta", help="fasta with .fai index as well.", type=str, required=True)
parser.add_argument("-o", "--output-name", help="output name of the output contigs file", type=str,required=True)
parser.add_argument("-v", "--vcf", help="vcf file with tabix index.", type=str, required=True)
parser.add_argument("-t", "--vartigs", help = "vartigs file.", type = str, required=True)
args = parser.parse_args()

vartig_file = args.vartigs
fasta_file = args.fasta
vcf_file = args.vcf
out = args.output_name

# read the file containing contigs and ranges
with open(vartig_file, 'r') as f:
    lines = [line.strip() for line in f]

# read the fasta file
fasta = Fasta(fasta_file)

# read the vcf file
vcf = pysam.VariantFile(vcf_file)

open(out,'w')
for i in range(0,len(lines),2):
    # parse the line
    hapid, contig_name, snprange, baserange, cov, err, hapq, rel_err = lines[i].split('\t')

    # get contig and ranges
    contig_name = contig_name.split(':')[1]
    snprange = list(map(int, snprange.split(':')[1].split('-')))
    baserange = list(map(int, baserange.split(':')[1].split('-')))

    # get fasta sequence
    sequence = fasta.get_seq(contig_name, baserange[0], baserange[1]).seq
    sequence = list(sequence)  # convert string to list for mutability

    # get vcf records
    snp_positions = {}
    for j, rec in enumerate(vcf.fetch(contig_name, baserange[0]-1, baserange[1]), start=1):
        if j in range(snprange[0], snprange[1] + 1):
            # storing the SNP positions with the alternate allele
            snp_positions[rec.pos + 1 - baserange[0]] = [rec.ref] + list(rec.alts)  # assuming first alt is the desired one

    # replace the bases in sequence with alternate alleles
    vartig = lines[i+1]
    for j, (pos, alt) in enumerate(snp_positions.items()):
        if vartig[j] == '?':
            sequence[pos - 1] = 'N' # adjusting for 0-based indexing in Python
        else:
            sequence[pos - 1] = alt[int(vartig[j])]  # adjusting for 0-based indexing in Python

    # write the modified sequence to new fasta record
    with open(out,'a') as f:
        # truncate sequence to lie between the first and fourteenth SNP
        truncated_sequence = ''.join(sequence)
        f.write(hapid)
        f.write('\n')
        f.write(truncated_sequence)
        f.write('\n')
print("Completed and written to output file")
