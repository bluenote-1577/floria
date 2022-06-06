import sys
import subprocess
import os
import glob
from joblib import Parallel, delayed


import re
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

minimap2_bin = 'minimap2'
abyss_pe_bin = 'abyss-pe'

fastq_files = glob.glob(partition_folder+'/*/short_reads/*.fastq')
fastq_files = [x for x in fastq_files if 'paired' in x]
fastq_files1 = [x for x in fastq_files if 'paired1' in x]
fastq_files2 = [x for x in fastq_files if 'paired2' in x]
fastq_files1 = sorted_nicely(fastq_files1)
fastq_files2 = sorted_nicely(fastq_files2)

rm_command = f'rm -r {folder}/intermediate/'
subprocess.run(rm_command, shell=True)
mkdir_string = f'mkdir {folder}; mkdir {folder}/intermediate'
subprocess.run(mkdir_string, shell=True)

def process(i):
    f1 = fastq_files1[i]
    f2 = fastq_files2[i]
    mkdir_command = f'mkdir {folder}/intermediate/p{i}_output'
    subprocess.run(mkdir_command, shell=True)
    
    abyss_command = f'{abyss_pe_bin} contigs k=25 name={i} in=\'../../../{f1} ../../../{f2} \' -C {folder}/intermediate/p{i}_output -j 10'
    subprocess.run(abyss_command, shell=True)
    output_bam = f'{folder}/assembly_p{i}.bam'
    abyss_contigs = f'{folder}/intermediate/p{i}_output/{i}-contigs.fa'
    abyss_unitigs = f'{folder}/intermediate/p{i}_output/{i}-3.fa'
    cont = True
    if os.path.isfile(abyss_contigs):
        final_contigs = abyss_contigs
    elif os.path.isfile(abyss_unitigs):
        final_contigs = abyss_unitigs
    else:
        cont = False
    if cont:
        rename_contig_cmd = f'sed -i \'s/>/>{i}_/g\' {final_contigs}'
        subprocess.run(rename_contig_cmd, shell=True, check=True)
        minimap_cmd = f'{minimap2_bin} -a {ref} {final_contigs} | samtools sort -o {folder}/intermediate/assembly_p{i}.bam'
        subprocess.run(minimap_cmd, shell=True, check=True)
        subprocess.run(f'samtools index {folder}/intermediate/assembly_p{i}.bam', shell=True, check=True)
results = Parallel(n_jobs=5)(delayed(process)(i) for i in range(len(fastq_files1)))
merge_index_cmd = f'samtools merge -f {folder}/all_assemblies.bam {folder}/intermediate/assembly_*.bam; samtools index {folder}/all_assemblies.bam'
subprocess.run(merge_index_cmd, shell=True, check=True)
