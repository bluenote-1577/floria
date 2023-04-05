import pysam
import subprocess
import sys
import argparse

parser = argparse.ArgumentParser(description='Generate a new bam file with haplotagging information for a single contig. ')

parser.add_argument("-t", "--haploset", help="the haploset file to haplotag.", type=str, required=True)
parser.add_argument("-b", "--bam", help="bam file to haplotag.", type=str, required=True)
parser.add_argument("-o", "--output-name", help="output name of the bam file. '.bam' is appended as a file extension. ", type=str,required=True)
parser.add_argument("-n", "--name-contig", help="name of the contig to haplotag.", type=str, required=True)
parser.add_argument("-q", "--min-hapq", help = "minimum HAPQ threshold for haplotagging (default = 0)", type = int, default = 0) 
args = parser.parse_args()

query_name_to_read_part = dict()
read_part_file = args.haploset
bam_file = args.bam
new_name = args.output_name
contig_name = args.name_contig
min_hapq = args.min_hapq

bam = pysam.AlignmentFile(bam_file)

read_part = dict()

index = 0
hapq_good = False;
for line in open(read_part_file,'r'):
    if '#' in line:
        split = line.split();
        index = int(split[0][1:])
        hapq = int(split[-2].split(':')[-1])
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


new_bam_name = new_name + '.bam'
new_bam_file = pysam.AlignmentFile(new_bam_name, "wb", template=bam)

print(f"Tagging {bam_file} and outputting new bam to {new_name}.bam ...")

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

bam.close()
new_bam_file.close()

print(f"Indexing {new_bam_name} ...")
pysam.index(new_bam_name)

print(f"Done! HP:i tags are now added to {new_name}.bam")
#cmd = f"samtools index {new_name}.bam"
#stream = subprocess.Popen(cmd, shell = True)
