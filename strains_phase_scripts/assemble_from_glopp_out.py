import sys
import subprocess
import os
import glob

import re
minimap2_bin = 'minimap2'
wtdbg2_bin = 'wtdbg2'
wtpoa_bin  = 'wtpoa-cns'

def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.

    Required arguments:
    l -- The iterable to be sorted.

    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

partition_folder = sys.argv[1]
folder= sys.argv[2]
ref = sys.argv[3]

fastq_files = glob.glob(partition_folder+'/*/long_reads/*.fastq')
fastq_files = [x for x in fastq_files if 'paired' not in x]
fastq_files = sorted_nicely(fastq_files)
print(fastq_files)

rm_command = f'rm -r {folder}/intermediate/'
subprocess.run(rm_command, shell=True)
mkdir_string = f'mkdir {folder}; mkdir {folder}/intermediate'
subprocess.run(mkdir_string, shell=True)

for (i,fastq_file) in enumerate(fastq_files):
    mkdir_command = f'mkdir {folder}/intermediate/p{i}_output'
    subprocess.run(mkdir_command, shell=True)
    wtdbg2_command = f'{wtdbg2_bin} -x preset2 -i {fastq_file} -o {folder}/intermediate/p{i}_output/p{i} -t 20'
    subprocess.run(wtdbg2_command, shell=True, check=True)
    final_contigs = f'{folder}/intermediate/p{i}_output/p{i}_contigs.fa'
    cns_command = f'{wtpoa_bin} -t 20 -i {folder}/intermediate/p{i}_output/p{i}.ctg.lay.gz -fo {final_contigs}'
    subprocess.run(cns_command, shell=True, check=True)
    rename_contig_cmd = f'sed -i \'s/ctg/p{i}_ctg/g\' {final_contigs}'
    subprocess.run(rename_contig_cmd, shell=True, check=True)
    minimap_cmd = f'{minimap2_bin} -a {ref} {final_contigs} | samtools sort -o {folder}/intermediate/assembly_p{i}.bam'
    subprocess.run(minimap_cmd, shell=True, check=True)
    subprocess.run(f'samtools index {folder}/intermediate/assembly_p{i}.bam', shell=True, check=True)

merge_index_cmd = f'samtools merge -f {folder}/all_assemblies.bam {folder}/intermediate/assembly_*.bam; samtools index {folder}/all_assemblies.bam'
subprocess.run(merge_index_cmd, shell=True, check=True)
