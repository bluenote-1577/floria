import pysam
import glob
import os
import subprocess
import sys
import argparse
from natsort import natsorted
import re

parser = argparse.ArgumentParser(description='Generate a new bam which contains all haplotagging information over all contigs in the result folder. Consider the haplotag_bam.py script if you only care about a single contig.')

parser.add_argument("-d", "--result-directory", help="the haploset file to haplotag.", type=str, required=True)
parser.add_argument("-b", "--bam", help="bam file to haplotag.", type=str, required=True)
parser.add_argument("-o", "--output-name", help="output name of the bam file. '.bam' is appended as a file extension. ", type=str,required=True)
parser.add_argument("-q", "--min-hapq", help = "minimum HAPQ threshold for haplotagging (default = 0)", type = int, default = 0) 
args = parser.parse_args()


bam = pysam.AlignmentFile(args.bam)
new_bam_name = args.output_name + '.bam'
new_bam_file = pysam.AlignmentFile(new_bam_name, "wb", template=bam)
min_hapq = args.min_hapq

p = re.compile('COV:(\d*\.?\d+)')
snp_p = re.compile('BASERANGE:(\d+)-(\d+)')
hapq_p = re.compile('HAPQ:(\d+)')
index_p = re.compile('HAP(\d+)')

for dir_object in natsorted(glob.glob(args.result_directory + '/*')):
    if os.path.isdir(dir_object):

        contig_name = os.path.basename(dir_object)

        haploset_g = glob.glob(dir_object + '/*_haploset*')
        if len(haploset_g) == 0:
            haploset_g = glob.glob(dir_object + '/all_part.txt')
        if len(haploset_g) == 0:
            print(f"ERROR: could not find haploset file for contig {contig_name}. Skipping ...")
            continue

        haploset_file = haploset_g[0]

        #Start haplotagging
        read_part = dict()
        query_name_to_read_part = dict()
        hapq_good = False;
        for line in open(haploset_file,'r'):
            if '>' in line:
                index = int(index_p.findall(line)[0])
                hapq = int(hapq_p.findall(line)[0])
                if hapq >= min_hapq:
                    hapq_good = True
                    read_part[index] = set()
                else:
                    hapq_good = False
            else:
                if hapq_good:
                    qname = line.split()[0]
                    read_part[index].add(qname)
                    query_name_to_read_part[qname] = index

        print(f"Tagging for contig {contig_name} ...")
        for b in bam.fetch(until_eof=True,contig=contig_name):
            not_frag = True;
            if b.query_name in query_name_to_read_part:
                i = query_name_to_read_part[b.query_name]
                b.set_tag('HP',i,value_type='i')
                new_bam_file.write(b)
                not_frag=False
            #for i in read_part.keys():
            #    qnames = read_part[i]
            #    if b.query_name in qnames:
            #        b.set_tag('HP',i,value_type='i')
            #        new_bam_file.write(b)
            #        not_frag=False
            if not_frag:
                    new_bam_file.write(b)

new_bam_file.close()
print("Indexing output bam file...")
print(f"Done! HP:i tags are now added to {new_bam_name}")
pysam.index(new_bam_name)





