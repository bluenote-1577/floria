import pysam
import subprocess
import sys

if len(sys.argv) < 4:
    print("usage: haplotag_bam.py contig_part.txt original_bam.bam new_haplotagged_bam_name.bam (contig_name)")
    exit()

read_part_file = sys.argv[1]
bam_file = sys.argv[2]
new_name = sys.argv[3]
contig_name = None
if len(sys.argv) == 5:
    contig_name = sys.argv[4]

bam = pysam.AlignmentFile(bam_file)

read_part = dict()

index = 0
for line in open(read_part_file,'r'):
    if '#' in line:
        split = line.split(',');
        index = int(split[0][1:])
        read_part[index] = set()
    else:
        read_part[index].add(line.split()[0])


new_bam_name = new_name
new_bam_file = pysam.AlignmentFile(new_bam_name, "wb", template=bam)

for b in bam.fetch(until_eof=True,contig=contig_name):
    not_frag = True;
    for i in read_part.keys():
        qnames = read_part[i]
        if b.query_name in qnames:
            b.set_tag('HP',i,value_type='i')
            new_bam_file.write(b)
            not_frag=False
    if not_frag:
            new_bam_file.write(b)

bam.close()
