import pysam
import subprocess
import sys

if len(sys.argv) < 2:
    print("USAGE: read_part.txt bamfile.bam prefix_name contig_name")
    exit()

read_part_file = sys.argv[1]
bam_file = sys.argv[2]
pref_nam = sys.argv[3]
bam = pysam.AlignmentFile(bam_file)

read_part = []

count_i = -1 
for line in open(read_part_file,'r'):
    if '#' in line:
        read_part.append(set())
        count_i += 1
    else:
        read_part[count_i].add(line.split()[0])

ploidy = len(read_part)

obam_files = []
file_names = []

for i in range(ploidy):
    file_name = pref_nam+str(i)+".bam"
    file_names.append(file_name)
    obam_files.append(pysam.AlignmentFile(file_name, "wb", template=bam))

file_names.append(pref_nam+"-not_mapped.bam")
obam_files.append(pysam.AlignmentFile(pref_nam+"-not_mapped.bam", "wb", template=bam))

if len(sys.argv) < 4:
    fetch = bam.fetch(until_eof=True)
else:
    contig_name = sys.argv[4]
    fetch = bam.fetch(contig_name)
for b in fetch:
    not_frag = True;
    for i in range(ploidy):
        qnames = read_part[i]
        if b.query_name in qnames:
            obam_files[i].write(b)
            not_frag=False
    if not_frag:
        obam_files[-1].write(b);


for obam in obam_files:
    obam.close()

for file_name in file_names:
    subprocess.run("samtools index " + file_name, shell=True, check=True)

bam.close()
