use crate::types_structs::{build_frag, update_frag, Frag, HapBlock};
use fxhash::{FxHashMap, FxHashSet};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::{bam, bam::Read as DUMMY_NAME1};
use rust_htslib::{bcf, bcf::Read as DUMMY_NAME2};
use crate::utils_frags;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::LineWriter;
use std::io::Write;
use std::io::{self, BufRead};
use std::path::Path;

// The output is wrapped in a Result to allow matching on errors
// returns an Iterator to the Reader of the lines of the file.
//
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

// Given a frags.txt file specified as in H-PoP, we return a collection
// (vector) of fragments after processing it.
//
pub fn get_frags_container<P>(filename: P) -> Vec<Frag>
where
    P: AsRef<Path>,
{
    let mut all_frags = Vec::new();
    let mut counter = 0;

    //Make sure file is able to be read
    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(l) = line {
                let v: Vec<&str> = l.split('\t').collect();

                //First column is the # of blocks
                if let Ok(num_blocks) = v[0].parse::<i32>() {
                    //                    println!("{}",num_blocks);
                    let mut seqs = FxHashMap::default();
                    let mut quals = FxHashMap::default();
                    let mut positions = FxHashSet::default();
                    let mut list_of_positions = Vec::new();
                    let mut first_position = 1;
                    let mut last_position = 1;

                    // For each block, read it into a dictionary with corresp. base
                    for i in 0..num_blocks {
                        let index = i as usize;
                        let start_pos = v[2 * index + 2].parse::<usize>().unwrap();
                        if i == 0 {
                            first_position = start_pos;
                        }
                        for (j, c) in v[2 * index + 3].chars().enumerate() {
                            seqs.insert(start_pos + j, c.to_digit(10).unwrap() as usize);
                            list_of_positions.push(start_pos + j);
                            positions.insert(start_pos + j);
                            last_position = start_pos + j
                        }
                    }

                    let qual_string = v.last().unwrap().as_bytes();
                    for (i, key) in list_of_positions.iter().enumerate() {
                        //We usually have a 33 offset for phred qualities. Rust should throw an
                        //error here if this result is negative.
                        quals.insert(*key, qual_string[i] - 33);
                    }

                    let new_frag = Frag {
                        id: v[1].to_string(),
                        counter_id: counter,
                        seq_dict: seqs,
                        qual_dict: quals,
                        positions: positions,
                        first_position: first_position,
                        last_position: last_position,
                    };

                    all_frags.push(new_frag);
                    counter += 1
                } else {
                    panic!("Not a number found in first column");
                }
            }
        }
    }

    all_frags
}


//Write a vector of blocks into a file.
pub fn write_blocks_to_file<P>(filename: P, blocks: &Vec<HapBlock>, lengths: &Vec<usize>, snp_to_genome : &Vec<usize>, part : &Vec<FxHashSet<&Frag>>)
where
    P: AsRef<Path>,
{
    let ploidy = blocks[0].blocks.len();
    let file = File::create(filename).expect("Can't create file");
    let mut file = LineWriter::new(file);
    let mut length_prev_block = 1;
    let emptydict = FxHashMap::default();
    let unpolished_block = utils_frags::hap_block_from_partition(part);
    //dbg!(snp_to_genome.len(),lengths[0] + 1);

    for (i, block) in blocks.iter().enumerate() {
        file.write_all(b"**BLOCK**\n").unwrap();
        for pos in length_prev_block..length_prev_block + lengths[i] {
            write!(file, "{}:{}\t", pos,snp_to_genome[pos-1]).unwrap();
            //Write haplotypes
            for k in 0..ploidy {
                let allele_map = block.blocks[k].get(&pos).unwrap_or(&emptydict);
                //If a block has no coverage at a position, we write -1.
                if *allele_map == emptydict {
                    file.write_all(b"-1\t").unwrap();
                } else {
                    let best_allele = allele_map.iter().max_by_key(|entry| entry.1).unwrap().0;
                    write!(file, "{}\t", best_allele).unwrap();
                }
            }
            
            //Write stats
            for k in 0..ploidy{
                let allele_map_unpolish = unpolished_block.blocks[k].get(&pos).unwrap_or(&emptydict);
                if *allele_map_unpolish == emptydict{
                    write!(file, "NA\t").unwrap();
                }
                else{
                    let mut first = true;
                    for (site,count) in allele_map_unpolish{
                        if !first{
                            write!(file,"|").unwrap();
                        }
                        if first{
                            first = false;
                        }
                       write!(file, "{}:{}",site,count).unwrap();
                    }
                    write!(file,"\t").unwrap();
                }
            }
            write!(file, "\n").unwrap();
        }
        write!(file, "*****").unwrap();
        length_prev_block += lengths[i]
    }
}

//Given a vcf file and a bam file, we get a vector of frags.
pub fn get_frags_from_bamvcf<P>(vcf_file: P, bam_file: P) -> Vec<Frag>
where
    P: AsRef<Path>,
{
    //Get which SNPS correspond to which positions on the genome.
    let mut vcf = bcf::Reader::from_path(vcf_file).unwrap();
    let mut snp_counter = 1;
    let mut set_of_pos = FxHashSet::default();
    let mut pos_allele_map = FxHashMap::default();
    let mut pos_to_snp_counter_map = FxHashMap::default();
    let header = vcf.header();

    if header.sample_count() > 1 {
        panic!("More than 1 sample detected in header of vcf file; please use only 1 sample");
    }

    if header.contig_count() > 1 {
        panic!("More than 1 contig detected in header of vcf file; please use only 1 contig/reference per vcf file.");
    }

    for rec in vcf.records() {
        let unr = rec.unwrap();
        let alleles = unr.alleles();
        let mut is_snp = true;
        let mut al_vec = Vec::new();

        for allele in alleles.iter() {
            if allele.len() > 1 {
                is_snp = false;
                break;
            }
            al_vec.push(allele[0]);
        }

        if !is_snp {
            println!(
                "BAM : Variant at position {} is not a snp. Ignoring.",
                unr.pos()
            );
            continue;
        }

        set_of_pos.insert(unr.pos());
        pos_to_snp_counter_map.insert(unr.pos(), snp_counter);
        snp_counter += 1;
        pos_allele_map.insert(unr.pos(), al_vec);
    }

    let mut bam = bam::Reader::from_path(bam_file).unwrap();
    //This may be important : We assume that distinct reads have different names. I can see this
    //being a problem in some weird bad cases, so be careful.
    let mut id_to_frag = FxHashMap::default();
    let mut counter_id = 0;

    //Scan the pileup table for every position on the genome which contains a SNP to get the aligned reads corresponding to the SNP. TODO : There should be a way to index into the bam.pileup() object so we don't have to iterate through positions which we already know are not SNPs.
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let pos_genome = pileup.pos();

        if !set_of_pos.contains(&(pos_genome as i64)) {
            continue;
        }

        let snp_id = pos_to_snp_counter_map.get(&(pos_genome as i64)).unwrap();

        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {
                let flags = alignment.record().flags();
                let errors_mask = 1796;
                let secondary_mask = 256;

                //Erroneous alignment, skip
                if flags & errors_mask > 0 {
                    continue;
                }

                //Secondary alignment, skip
                if flags & secondary_mask > 0 {
                    continue;
                }

                let id_string = String::from_utf8(alignment.record().qname().to_vec()).unwrap();
                let id_string2 = id_string.clone();

                if !id_to_frag.contains_key(&id_string) {
                    counter_id += 1;
                }

                let frag_to_ins = build_frag(id_string, counter_id);

                let readbase = alignment.record().seq()[alignment.qpos().unwrap()];
                let qualbase = alignment.record().qual()[alignment.qpos().unwrap()];

                for (i, allele) in pos_allele_map
                    .get(&(pos_genome as i64))
                    .unwrap()
                    .iter()
                    .enumerate()
                {
                    //Only build the frag if the base is one of the SNP alleles.
                    if readbase == *allele {
                        let mut frag = id_to_frag.entry(id_string2).or_insert(frag_to_ins);
                        update_frag(&mut frag, i, qualbase, *snp_id);
                        break;
                    }
                }
            }
        }
    }

    let mut vec_frags = Vec::new();
    for (_id, frag) in id_to_frag.into_iter() {
        if frag.positions.len() > 1 {
            vec_frags.push(frag);
        }
    }

    vec_frags
}

//Read a vcf file to get the genotypes. We read genotypes into a dictionary of keypairs where the
//keys are positions, and the values are dictionaries which encode the genotypes. E.g. the genotype
//1 1 0 0 at position 5 would be (5,{1 : 2, 0 : 2}).
pub fn get_genotypes_from_vcf_hts<P>(vcf_file: P) -> (Vec<usize>, FxHashMap<usize, FxHashMap<usize, usize>>, usize)
where
    P: AsRef<Path>,
{
    let mut vcf = bcf::Reader::from_path(vcf_file).unwrap();
    let mut positions_vec = Vec::new();
    let mut genotype_dict = FxHashMap::default();
    let header = vcf.header();
    let mut snp_counter = 1;
    let mut vcf_ploidy = 0;

    if header.sample_count() > 1 {
        panic!("More than 1 sample detected in header of vcf file; please use only 1 sample");
    }

    if header.contig_count() > 1 {
        panic!("More than 1 contig detected in header of vcf file; please use only 1 contig/reference per vcf file.");
    }

    for rec in vcf.records() {
        let mut unr = rec.unwrap();
        let alleles = unr.alleles();
        let mut is_snp = true;

        for allele in alleles.iter() {
            if allele.len() > 1 {
                is_snp = false;
                break;
            }
        }

        if !is_snp {
            println!(
                "VCF : Variant at position {} is not a snp. Ignoring.",
                unr.pos()
            );
            continue;
        }

        let genotypes = unr.genotypes().unwrap().get(0);
        vcf_ploidy = genotypes.len();
        let mut genotype_counter = FxHashMap::default();
        for allele in genotypes.iter() {
            match allele {
                GenotypeAllele::Unphased(x) => {
                    let count = genotype_counter.entry(*x as usize).or_insert(0);
                    *count += 1
                }
                GenotypeAllele::Phased(x) => {
                    let count = genotype_counter.entry(*x as usize).or_insert(0);
                    *count += 1
                }
                GenotypeAllele::UnphasedMissing => {
                    for i in 0..4 {
                        let count = genotype_counter.entry(i).or_insert(0);
                        *count += 1;
                    }
                }
                GenotypeAllele::PhasedMissing => {
                    for i in 0..4 {
                        let count = genotype_counter.entry(i).or_insert(0);
                        *count += 1;
                    }
                }
            }
        }

        genotype_dict.insert(snp_counter,genotype_counter);
        //+1 because htslib is 0 index by default
        positions_vec.push(unr.pos() as usize + 1); 
        snp_counter += 1;
    }

    (positions_vec,genotype_dict,vcf_ploidy)
}

//Convert a fragment which stores sequences in a dictionary format to a block format which makes
//writing to frag files easier.
fn convert_dict_to_block(frag: Frag) -> (Vec<usize>, Vec<Vec<usize>>, Vec<u8>) {
    let d = frag.seq_dict;
    let vec_d: BTreeMap<usize, usize> = d.into_iter().collect();
    let vec_q: BTreeMap<usize, u8> = frag.qual_dict.into_iter().collect();
    let mut prev_pos = 0;
    let mut block_start_pos = Vec::new();
    let mut blocks = Vec::new();
    let mut block = Vec::new();
    let mut qual_block = Vec::new();

    for (pos, var) in &vec_d {
        if prev_pos == 0 {
            prev_pos = *pos;
            block.push(*var);
            block_start_pos.push(*pos);
        } else if pos - prev_pos > 1 {
            blocks.push(block);
            block = vec![*var];
            block_start_pos.push(*pos);
            prev_pos = *pos;
        } else if pos - prev_pos == 1 {
            block.push(*var);
            prev_pos = *pos;
        }
    }

    for (_pos, q) in &vec_q {
        qual_block.push(*q);
    }

    blocks.push(block);
    (block_start_pos, blocks, qual_block)
}

//Write a vector of sorted fragment files by first position (no guarantees on end position) to a
//file in the same format as H-PoP and other haplotypers.
pub fn write_frags_file(frags: Vec<Frag>, filename: String) {
    let file = File::create(filename).expect("Can't create file");
    let mut file = LineWriter::new(file);
    for frag in frags.into_iter() {
        let frag_id = frag.id.clone();
        let (start_vec, blocks, qual_block) = convert_dict_to_block(frag);
        if start_vec.len() != blocks.len() {
            dbg!(start_vec.len(), blocks.len());
            panic!("Block length diff");
        }

        write!(file, "{}\t", blocks.len()).unwrap();
        write!(file, "{}\t", frag_id).unwrap();
        for i in 0..blocks.len() {
            write!(file, "{}\t", start_vec[i]).unwrap();
            for var in blocks[i].iter() {
                write!(file, "{}", *var).unwrap();
            }
            write!(file, "\t").unwrap();
        }

        for q in qual_block.iter() {
            write!(file, "{}", (*q + 33) as char).unwrap();
        }

        write!(file, "\n").unwrap();
    }
}

