from sys import argv

vcf_file = argv[1]
refs = set()
for line in open(vcf_file,'r'):
    if line == "":
        continue
    if line[0] == '#':
        continue
    ref_chrom = line.split()[0]
    print(ref_chrom)
    refs.add(ref_chrom)

refs = list(refs)
refs.sort()
print(refs)

new_vcf =  open(vcf_file+"c_header",'w')
count = 0
for line in open(vcf_file,'r'):
    if count != 2:
        new_vcf.write(line)
    else:
        for ref in refs:
            new_vcf.write("##contig=<ID="+ref+">\n")
        new_vcf.write(line)
    count += 1
    





