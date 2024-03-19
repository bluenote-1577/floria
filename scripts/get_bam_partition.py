import argparse
import re
import pysam

p = re.compile('COV:(\d*\.?\d+)')
snp_p = re.compile('BASERANGE:(\d+)-(\d+)')
hapq_p = re.compile('HAPQ:(\d+)')
index_p = re.compile('HAP(\d+)')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Partition a BAM file into multiple files based on a partition file.')
    parser.add_argument('-t','--haplosets', help='1 or more haploset files', required=True, type=str, nargs='+')
    parser.add_argument('-b','--bam_file', help='Input BAM file',required=True, type=str)
    parser.add_argument('-p', '--prefix_name', help='Prefix for output BAM files', default='split-bam', type=str)
    parser.add_argument("-q", "--min-hapq", help = "minimum HAPQ threshold for haplotagging (default = 0)", type = int, default = 0) 
    return parser.parse_args()

hapq_good = True
def read_partitions(read_part_file, min_hapq):
    read_parts = []
    with open(read_part_file, 'r') as file:
        for line in file:
            if '>' in line:
                hapq = int(hapq_p.findall(line)[0])
                index = int(index_p.findall(line)[0])
                if hapq >= min_hapq:
                    hapq_good = True
                    read_parts.append({'index':index,'hapq':hapq,'set':set()})
                else:
                    hapq_good = False
            else:
                if hapq_good:
                    read_parts[-1]['set'].add(line.split()[0])
    return read_parts

def create_output_files(prefix_name, indices, template_bam):
    obam_files = [(i,pysam.AlignmentFile(f"{prefix_name}{i}.bam", "wb", template=template_bam)) for i in indices]
    #create dictionary from this tuple
    obam_files = dict(obam_files)
    return obam_files

def partition_bam(read_parts, bam_file, obam_files, haploset_file):
    with pysam.AlignmentFile(bam_file) as bam:
        #get contig name EXAMPLE.test.test from haploset file folder/EXAMPLE.test.test.haploset
        #contig_name = haploset_file.split('/')[-1].split('.')[0:-1].join('.')
        contig_name = '.'.join(haploset_file.split('/')[-1].split('.')[0:-1])
        fetch = bam.fetch(contig_name) if contig_name else bam.fetch(until_eof=True)
        for b in fetch:
            not_frag = True
            for d in read_parts:
                qnames = d['set']
                i = d['index']
                if b.query_name in qnames:
                    obam_files[i].write(b)
                    not_frag = False
                    break

def main():
    args = parse_arguments()
    for haploset in args.haplosets:
        print("Splitting bam file for", haploset, "with", args.bam_file)
        read_parts = read_partitions(haploset,args.min_hapq)
        obam_files = create_output_files(args.prefix_name, [x['index'] for x in read_parts], pysam.AlignmentFile(args.bam_file))

        partition_bam(read_parts, args.bam_file, obam_files, haploset)

        for obam in obam_files.values():
            obam.close()
            pysam.index(obam.filename.decode())  # Indexing using pysam

        print("Splitting complete")

if __name__ == '__main__':
    main()

