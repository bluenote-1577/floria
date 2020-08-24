use crate::types_structs::{Frag, HapBlock,build_frag,update_frag};
use fxhash::{FxHashMap, FxHashSet};
use rust_htslib::{bam, bam::Read as DUMMY_NAME1};
use rust_htslib::{bcf, bcf::Read as DUMMY_NAME2};
use std::fs::File;
use std::io::LineWriter;
use std::io::Write;
use std::io::{self, BufRead};
use std::path::Path;
use std::collections::BTreeMap;

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

// Given a SORTED!!! frags.txt file specified as in H-PoP, we return a collection
// (vector) of fragments after processing it
pub fn get_frags_container<P>(filename: P) -> Vec<Frag>
where
    P: AsRef<Path>,
{
    let mut all_frags = Vec::new();
    let mut counter = 0;
    let mut prev_first_pos = 0;

    //Make sure file is able to be read
    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(l) = line {
                let v: Vec<&str> = l.split('\t').collect();
                //                println!("{:?}",v);

                //First column is the # of blocks
                if let Ok(num_blocks) = v[0].parse::<i32>() {
                    //                    println!("{}",num_blocks);
                    let mut seqs = FxHashMap::default();
                    let mut quals = FxHashMap::default();
                    let mut positions = FxHashSet::default();
                    let mut list_of_positions = Vec::new();
                    let mut first_position = 1;
                    let mut last_position = 1;

                    // For each block, read it into a dictionary with corresp.
                    // base
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
                        //Offset for phred qualities
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

                    if first_position < prev_first_pos {
                        panic!("Frags file is not sorted.")
                    }

                    prev_first_pos = first_position;
                    counter += 1
                } else {
                    panic!("Not a number found in first column");
                }
            }
        }
    }
    all_frags
}

//Read a vcf file to get the genotypes : TODO may have to refactor this method into a more general
//vcf method
pub fn get_genotypes_from_vcf<P>(filename: P) -> (FxHashMap<usize, FxHashMap<usize, usize>>, usize)
where
    P: AsRef<Path>,
{
    let mut genotype_dict = FxHashMap::default();
    let mut ploidy = 0;
    if let Ok(lines) = read_lines(filename) {
        let mut counter = 1;
        for line in lines {
            if let Ok(l) = line {
                //Skip the first few files of the VCF
                if l.contains('#') {
                    continue;
                }

                let v: Vec<&str> = l.split('\t').collect();
                let v: Vec<&str> = v.last().unwrap().split(':').collect();
                let genotypes = v[0];
                let mut split_genotypes = Vec::new();

                if genotypes.contains('|') {
                    split_genotypes = genotypes.split('|').collect();
                } else if genotypes.contains('/') {
                    split_genotypes = genotypes.split('/').collect();
                } else {
                    panic!("Genotype column not processed correctly : {}", genotypes);
                }
                ploidy = split_genotypes.len();
                let mut genotype_value = FxHashMap::default();
                for allele in split_genotypes.iter() {
                    let num_allele: usize = allele.parse().unwrap();
                    let val = genotype_value.entry(num_allele).or_insert(0);
                    *val += 1
                }
                //                println!("{:?}",genotype_value);
                genotype_dict.insert(counter, genotype_value);
                counter += 1;
            }
        }
    }
    (genotype_dict, ploidy)
}

//Write a vector of blocks into a file
pub fn write_blocks_to_file<P>(filename: P, blocks: &Vec<HapBlock>, lengths: &Vec<usize>)
where
    P: AsRef<Path>,
{
    let ploidy = blocks[0].blocks.len();
    let file = File::create(filename).expect("Can't create file");
    let mut file = LineWriter::new(file);
    let mut length_prev_block = 1;
    let emptydict = FxHashMap::default();
    for (i, block) in blocks.iter().enumerate() {
        file.write_all(b"**BLOCK**\n").unwrap();
        for pos in length_prev_block..length_prev_block + lengths[i] {
            write!(file, "{}\t", pos).unwrap();
            for k in 0..ploidy {
                let allele_map = block.blocks[k].get(&pos).unwrap_or(&emptydict);
                if *allele_map == emptydict {
                    file.write_all(b"-1\t").unwrap();
                } else {
                    let best_allele = allele_map.iter().max_by_key(|entry| entry.1).unwrap().0;
                    write!(file, "{}\t", best_allele).unwrap();
                }
            }
            write!(file, "\n").unwrap();
        }
        write!(file, "*****").unwrap();
        length_prev_block += lengths[i]
    }
}

pub fn get_frags_from_bamvcf<P>(vcf_file: P, bam_file: P) -> Vec<Frag>
where
    P: AsRef<Path>,
{
    //Get which SNPS correspond to which positions
    let mut vcf = bcf::Reader::from_path(vcf_file).unwrap();
    let mut snp_counter = 1;
    let mut set_of_pos = FxHashSet::default();
    let mut pos_allele_map = FxHashMap::default();
    let mut pos_to_snp_counter_map = FxHashMap::default();
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
            println!("Variant at position {} is not a snp. Ignoring.", unr.pos());
            continue;
        }

        set_of_pos.insert(unr.pos());
        pos_to_snp_counter_map.insert(unr.pos(),snp_counter);
        snp_counter +=1;
        pos_allele_map.insert(unr.pos(),al_vec);
    }

    let mut bam = bam::Reader::from_path(bam_file).unwrap();
    let mut id_to_frag = FxHashMap::default();
    let mut counter_id = 0;

    for p in bam.pileup() {
        let pileup = p.unwrap();
        let pos_genome = pileup.pos();

        if !set_of_pos.contains(&(pos_genome as i64)){
            continue;
        }

        let snp_id = pos_to_snp_counter_map.get(&(pos_genome as i64)).unwrap();

        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {

                //Chimeric or secondary alignments do this with minimap
                if alignment.record().seq().len() == 0{
                    continue
                }

                let id_string = String::from_utf8(alignment.record().qname().to_vec()).unwrap();
                let id_string2 = id_string.clone();

                if !id_to_frag.contains_key(&id_string){
                    counter_id += 1;
                }

                let frag_to_ins = build_frag(id_string,counter_id);


                let readbase = alignment.record().seq()[alignment.qpos().unwrap()];
                let qualbase = alignment.record().qual()[alignment.qpos().unwrap()];

                for (i,allele) in pos_allele_map.get(&(pos_genome as i64)).unwrap().iter().enumerate(){
                    if readbase == *allele{
                        let mut frag = id_to_frag.entry(id_string2).or_insert(frag_to_ins);
                        update_frag(&mut frag,i,qualbase,*snp_id);
                        break
                    }
                }
            }
        }
    }

    let mut vec_frags = Vec::new();
    for (_id,frag) in id_to_frag.into_iter(){
        if frag.positions.len() > 1{
            vec_frags.push(frag);
        }
    }

    vec_frags.sort_by(|a,b| a.first_position.cmp(&b.first_position));
    vec_frags

}

fn convert_dict_to_block(frag :&Frag) -> (Vec<usize>, Vec<Vec<usize>>) {
    let d = frag.seq_dict;
    let vec_d : BTreeMap<usize,usize> = d.into_iter().collect();
    let mut prev_pos = 0;
    let mut block_start_pos = Vec::new();
    let mut blocks = Vec::new();
    let mut block = Vec::new();
    for (pos,var) in &vec_d{
        if prev_pos == 0{
            prev_pos = *pos;
            block.push(*var);
            block_start_pos.push(pos);
        }
        
        else if pos-prev_pos > 1{
            blocks.push(block);
            block = vec![*var];
            block_start_pos.push(pos);
            prev_pos = *pos;
        }

        else if pos - prev_pos == {
            block.push(*var);
        }

        else{
            panic!("Something went wrong here");
        }
    }

    blocks.push(block);
    blocks
}

pub fn write_frags_file(frags : &Vec<Frag>,filename : String){
    let file = File::create(filename).expect("Can't create file");
    let mut file = LineWriter::new(file);
    for frag in frags.iter(){
    }
//    let emptydict = FxHashMap::default();
//    for (i, block) in blocks.iter().enumerate() {
//        file.write_all(b"**BLOCK**\n").unwrap();
//        for pos in length_prev_block..length_prev_block + lengths[i] {
//            write!(file, "{}\t", pos).unwrap();
//            for k in 0..ploidy {
//                let allele_map = block.blocks[k].get(&pos).unwrap_or(&emptydict);
//                if *allele_map == emptydict {
//                    file.write_all(b"-1\t").unwrap();
//                } else {
//                    let best_allele = allele_map.iter().max_by_key(|entry| entry.1).unwrap().0;
//                    write!(file, "{}\t", best_allele).unwrap();
//                }
//            }
//            write!(file, "\n").unwrap();
//        }
//        write!(file, "*****").unwrap();
//        length_prev_block += lengths[i]
//    }
}
