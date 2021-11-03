use crate::types_structs::{build_frag, update_frag, Frag, HapBlock};
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use rust_htslib::bam::header::Header;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::HeaderView as HeaderViewBam;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::{bam, bam::Read as DUMMY_NAME1};
use rust_htslib::{bcf, bcf::Read as DUMMY_NAME2};
use std::collections::BTreeMap;
use std::fs;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::LineWriter;
use std::io::Write;
use std::io::{self, BufRead};
use std::mem;
use std::path::Path;
use std::str;

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
pub fn get_frags_container<P>(filename: P) -> FxHashMap<String, Vec<Frag>>
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
                        supp_aln: None,
                    };

                    all_frags.push(new_frag);
                    counter += 1
                } else {
                    panic!("Not a number found in first column");
                }
            }
        }
    }

    let mut frags_map = FxHashMap::default();
    frags_map.insert(String::from("frag_contig"), all_frags);
    frags_map
}

//Write a vector of blocks into a file.
pub fn write_blocks_to_file<P>(
    out_part_dir: P,
    blocks: &Vec<HapBlock>,
    lengths: &Vec<usize>,
    snp_to_genome: &Vec<usize>,
    part: &Vec<FxHashSet<&Frag>>,
    _first_iter: bool,
    contig: &String,
    break_positions: &FxHashMap<usize,FxHashSet<usize>>
) where
    P: AsRef<Path>,
{
    let ploidy = blocks[0].blocks.len();
    let filename = out_part_dir
        .as_ref()
        .join(format!("{}_phasing.txt", contig));

    let file;
    file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(filename)
        .unwrap();
    //let file = File::create(filename).expect("Can't create file");
    let mut file = LineWriter::new(file);
    let mut length_prev_block = 1;
    let emptydict = FxHashMap::default();
    let unpolished_block = utils_frags::hap_block_from_partition(part);
    //dbg!(snp_to_genome.len(),lengths[0] + 1);

    for (i, block) in blocks.iter().enumerate() {
        let title_string = format!("**{}**\n", contig);
        file.write_all(title_string.as_bytes()).unwrap();
        for pos in length_prev_block..length_prev_block + lengths[i] {
            if break_positions.contains_key(&pos){
                write!(file, "--------\n").unwrap();
            }
            if snp_to_genome.len() == 0 {
                write!(file, "{}:NA\t", pos).unwrap();
            } else {
                write!(file, "{}:{}\t", pos, snp_to_genome[pos - 1]).unwrap();
            }
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
            for k in 0..ploidy {
                let allele_map_unpolish =
                    unpolished_block.blocks[k].get(&pos).unwrap_or(&emptydict);
                if *allele_map_unpolish == emptydict {
                    write!(file, "NA\t").unwrap();
                } else {
                    let mut first = true;
                    for (site, count) in allele_map_unpolish {
                        if !first {
                            write!(file, "|").unwrap();
                        }
                        if first {
                            first = false;
                        }
                        write!(file, "{}:{}", site, count).unwrap();
                    }
                    write!(file, "\t").unwrap();
                }
            }
            write!(file, "\n").unwrap();
        }
        write!(file, "*****\n").unwrap();
        length_prev_block += lengths[i]
    }
}

//Given a vcf file and a bam file, we get a vector of frags.
pub fn get_frags_from_bamvcf<P>(vcf_file: P, bam_file: P, filter_supplementary: bool) -> FxHashMap<String, Vec<Frag>>
where
    P: AsRef<Path>,
{
    //Get which SNPS correspond to which positions on the genome.
    let mut vcf = match bcf::Reader::from_path(vcf_file) {
        Ok(vcf) => vcf,
        Err(_) => panic!("rust_htslib had an error reading the VCF file. Exiting."),
    };
    let mut snp_counter = 1;
    let mut vcf_set_of_pos = FxHashMap::default();
    let mut vcf_pos_allele_map = FxHashMap::default();
    let mut vcf_pos_to_snp_counter_map = FxHashMap::default();
    let mut all_set_of_pos = FxHashSet::default();
    //    let mut set_of_pos = FxHashSet::default();
    //    let mut pos_allele_map = FxHashMap::default();
    //    let mut pos_to_snp_counter_map = FxHashMap::default();
    let vcf_header = vcf.header().clone();

    //    if header.contig_count() > 1 {
    //        panic!("More than 1 contig detected in header of vcf file; please use only 1 contig/reference per vcf file.");
    //    }

    let mut last_ref_chrom: &[u8] = &[];
    for rec in vcf.records() {
        let unr = rec.unwrap();
        let alleles = unr.alleles();
        let mut al_vec = Vec::new();
        let mut is_snp = true;

        let record_rid = unr.rid().unwrap();
        let ref_chrom_vcf = vcf_header.rid2name(record_rid).unwrap();
        //dbg!(String::from_utf8_lossy(ref_chrom_vcf));
        if last_ref_chrom != ref_chrom_vcf {
            snp_counter = 1;
            last_ref_chrom = ref_chrom_vcf;
        }
        let set_of_pos = vcf_set_of_pos
            .entry(ref_chrom_vcf)
            .or_insert(FxHashSet::default());
        let pos_allele_map = vcf_pos_allele_map
            .entry(ref_chrom_vcf)
            .or_insert(FxHashMap::default());
        let pos_to_snp_counter_map = vcf_pos_to_snp_counter_map
            .entry(ref_chrom_vcf)
            .or_insert(FxHashMap::default());

        for allele in alleles.iter() {
            if allele.len() > 1 {
                is_snp = false;
                break;
            }
            al_vec.push(allele[0]);
        }

        if !is_snp {
            //            println!(
            //                "BAM : Variant at position {} is not a snp. Ignoring.",
            //                unr.pos()
            //            );
            continue;
        }

        set_of_pos.insert(unr.pos());
        all_set_of_pos.insert(unr.pos());
        pos_to_snp_counter_map.insert(unr.pos(), snp_counter);
        snp_counter += 1;
        pos_allele_map.insert(unr.pos(), al_vec);
    }

//    let mut bam = bam::Reader::from_path(bam_file).unwrap();
    let mut bam = match bam::Reader::from_path(bam_file) {
        Ok(bam) => bam,
        Err(_) => panic!("rust_htslib had an error while reading the BAM file. Exiting"),
    };

    //Check the headers to see how many references there are.
    let header = Header::from_template(bam.header());
    let bam_header_view = HeaderViewBam::from_header(&header);
//    let mut number_references = 0;
//    for (id, content) in header.to_hashmap() {
//        if id == "SQ" {
//            number_references = content.len();
//        }
//    }

    //    if number_references > 1{
    //        dbg!(vcf_pos_to_snp_counter_map.keys());
    //    }

    //This may be important : We assume that distinct reads have different names. I can see this
    //being a problem in some weird bad cases, so be careful.
    let mut ref_id_to_frag = FxHashMap::default();
    let mut counter_id = 0;

    //Scan the pileup table for every position on the genome which contains a SNP to get the aligned reads corresponding to the SNP. TODO : There should be a way to index into the bam.pileup() object so we don't have to iterate through positions which we already know are not SNPs.
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let pos_genome = pileup.pos();

        if !all_set_of_pos.contains(&(pos_genome as i64)) {
            continue;
        }

        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {
                let aln_record = alignment.record();
                let tid = aln_record.tid();
                let ref_chrom = bam_header_view.tid2name(tid as u32);
                //dbg!(String::from_utf8_lossy(ref_chrom));
                let get_ref_chrom = vcf_pos_to_snp_counter_map.get(ref_chrom);

                let pos_to_snp_counter_map = match get_ref_chrom {
                    Some(pos_to_snp_counter_map) => pos_to_snp_counter_map,
                    None => continue,
                };

                if pos_to_snp_counter_map.contains_key(&(pos_genome as i64)) == false {
                    continue;
                }
                let snp_id = pos_to_snp_counter_map.get(&(pos_genome as i64)).unwrap();
                let flags = aln_record.flags();
                let errors_mask = 1796;
                let secondary_mask = 256;
                let supplementary_mask = 2048;
                let mapq_supp_cutoff = 30;
                let mapq_normal_cutoff = 5;
                let id_to_frag = ref_id_to_frag
                    .entry(ref_chrom)
                    .or_insert(FxHashMap::default());
                let mapq = aln_record.mapq();

                let id_string = String::from_utf8(aln_record.qname().to_vec()).unwrap();
                //Erroneous alignment, skip
                if flags & errors_mask > 0 {
                    //dbg!(&flags,&id_string);
                    continue;
                }

                //Secondary alignment, skip
                if flags & secondary_mask > 0 {
                    //dbg!(&flags,&id_string);
                    continue;
                }

                if flags & supplementary_mask > 0{
                    if filter_supplementary{
                        if mapq < mapq_supp_cutoff{
                            continue;
                        }
                    }
                }

                let id_string2 = id_string.clone();

                let mut supp_aln = None;
                if !id_to_frag.contains_key(&id_string) {
                    if let Ok(tag) = aln_record.aux(b"SA") {
                        if let Aux::String(v) = tag {
                            let str_v = v;
//                            let str_v = str::from_utf8(v).unwrap();
                            let split: Vec<&str> = str_v.split(';').collect();
                            let best_SA: Vec<&str> = split[0].split(',').collect();
                            let best_contig = best_SA[0];
                            supp_aln = Some(String::from(best_contig));
                        }
                    }
                    counter_id += 1;
                }

                let frag_to_ins = build_frag(id_string, counter_id, supp_aln);

                let readbase = alignment.record().seq()[alignment.qpos().unwrap()];
                let qualbase = alignment.record().qual()[alignment.qpos().unwrap()];

                for (i, allele) in vcf_pos_allele_map
                    .get(ref_chrom)
                    .unwrap()
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

    let mut ref_vec_frags = FxHashMap::default();
    let mut keys = FxHashSet::default();
    for ref_chrom in ref_id_to_frag.keys() {
        ref_vec_frags.insert(String::from_utf8(ref_chrom.to_vec()).unwrap(), Vec::new());
        keys.insert(ref_chrom.clone());
    }
    for ref_chrom in keys {
        let vec_frags = ref_vec_frags
            .get_mut(&String::from_utf8(ref_chrom.to_vec()).unwrap())
            .unwrap();
        let id_to_frag = ref_id_to_frag.get_mut(ref_chrom).unwrap();
        let id_to_frag = mem::replace(id_to_frag, FxHashMap::default());
        for (_id, frag) in id_to_frag.into_iter() {
            //IMPORTANT: I'm turning this off for metagenomics because some fragments may only
            //index one read. However, this is still useful because we don't know ploidy info. 
            if frag.positions.len() > 0 {
                vec_frags.push(frag);
            }
        }
    }

    ref_vec_frags
}

//Read a vcf file to get the genotypes. We read genotypes into a dictionary of keypairs where the
//keys are positions, and the values are dictionaries which encode the genotypes. E.g. the genotype
//1 1 0 0 at position 5 would be (5,{1 : 2, 0 : 2}).
pub fn get_genotypes_from_vcf_hts<P>(
    vcf_file: P,
) -> (
    FxHashMap<String, Vec<usize>>,
    FxHashMap<String, FxHashMap<usize, FxHashMap<usize, usize>>>,
    usize,
)
where
    P: AsRef<Path>,
{
    let mut vcf = match bcf::Reader::from_path(vcf_file) {
        Ok(vcf) => vcf,
        Err(_) => panic!("rust_htslib had an error while reading the BAM file. Exiting."),
    };
    let mut map_positions_vec = FxHashMap::default();
    let mut map_genotype_dict = FxHashMap::default();
    //let mut positions_vec = Vec::new();
    //let mut genotype_dict = FxHashMap::default();
    let header = vcf.header().clone();
    let mut snp_counter = 1;
    let mut vcf_ploidy = 0;

    if header.sample_count() > 1 {
        panic!("More than 1 sample detected in header of vcf file; please use only 1 sample");
    }

    //    if header.contig_count() > 1 {
    //        panic!("More than 1 contig detected in header of vcf file; please use only 1 contig/reference per vcf file.");
    //    }

    let mut last_ref_chrom: &[u8] = &[];

    for rec in vcf.records() {
        let mut unr = rec.unwrap();
        let alleles = unr.alleles();
        let mut is_snp = true;
        let record_rid = unr.rid().unwrap();
        let ref_chrom_vcf = header.rid2name(record_rid).unwrap();
        if last_ref_chrom != ref_chrom_vcf {
            last_ref_chrom = ref_chrom_vcf;
            snp_counter = 1;
        }

        for allele in alleles.iter() {
            if allele.len() > 1 {
                is_snp = false;
                break;
            }
        }

        if !is_snp {
            //            println!(
            //                "VCF : Variant at position {} is not a snp. Ignoring.",
            //                unr.pos()
            //            );
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

        let genotype_dict = map_genotype_dict
            .entry(String::from_utf8(ref_chrom_vcf.to_vec()).unwrap())
            .or_insert(FxHashMap::default());
        let positions_vec = map_positions_vec
            .entry(String::from_utf8(ref_chrom_vcf.to_vec()).unwrap())
            .or_insert(Vec::new());
        genotype_dict.insert(snp_counter, genotype_counter);
        //+1 because htslib is 0 index by default
        positions_vec.push(unr.pos() as usize + 1);
        snp_counter += 1;
    }

    (map_positions_vec, map_genotype_dict, vcf_ploidy)
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

pub fn write_output_partition_to_file<P>(
    part: &Vec<FxHashSet<&Frag>>,
    out_bam_part_dir: P,
    contig: &String,
    break_positions: &FxHashMap<usize,FxHashSet<usize>>
) where
    P: AsRef<Path>,
{
    fs::create_dir_all(&out_bam_part_dir).unwrap();
    let contig_path = out_bam_part_dir
        .as_ref()
        .join(format!("{}_part.txt", contig));
    let file = File::create(contig_path).expect("Can't create file");
    let mut file = LineWriter::new(file);

    for (i, set) in part.iter().enumerate() {
        let mut vec_part: Vec<&&Frag> = set.into_iter().collect();
        vec_part.sort_by(|a, b| a.first_position.cmp(&b.first_position));
        write!(file, "#{}\n", i).unwrap();
        for frag in vec_part {
            write!(file, "{}\t{}\t{}\n", frag.id.clone(), frag.first_position, frag.last_position).unwrap();
        }
    }
}
