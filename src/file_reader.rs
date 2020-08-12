use std::fs::File;
use fnv::FnvHashMap;
use fnv::FnvHashSet;
use std::io::{self, BufRead};
use std::path::Path;
use crate::types_structs::Frag;

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

// Given a SORTED!!! frags.txt file specified as in H-PoP, we return a collection
// (vector) of fragments after processing it
pub fn get_frags_container<P>(filename : P) -> Vec<Frag> where P: AsRef <Path>,{
    let mut all_frags = Vec::new();
    let mut counter = 0;
    let mut prev_first_pos = 0;
    
    //Make sure file is able to be read
    if let Ok(lines) = read_lines(filename){
        for line in lines{
            if let Ok(l) = line{
                let v : Vec<&str> = l.split('\t').collect();
//                println!("{:?}",v);

                //First column is the # of blocks
                if let Ok(num_blocks) = v[0].parse::<i32>(){

//                    println!("{}",num_blocks);
                    let mut seqs = FnvHashMap::default();
                    let mut quals = FnvHashMap::default();
                    let mut positions = FnvHashSet::default();
                    let mut list_of_positions = Vec::new();
                    let mut first_position = 1;
                    let mut last_position = 1;

                    // For each block, read it into a dictionary with corresp.
                    // base
                    for i in 0..num_blocks {
                        let index = i as usize;
                        let start_pos = v[2*index + 2].parse::<usize>().unwrap();
                        if i == 0 {
                            first_position = start_pos;
                        }
                        for (j,c) in v[2*index+3].chars().enumerate(){
                            seqs.insert(start_pos+j,c.to_digit(10).unwrap() as usize);
                            list_of_positions.push(start_pos+j);
                            positions.insert(start_pos+j);
                            last_position = start_pos+j
                        }
                    }

                    let qual_string = v.last().unwrap().as_bytes();
                    for (i,key) in list_of_positions.iter().enumerate(){
                        //Offset for phred qualities
                        quals.insert(*key,qual_string[i] - 33);
                    }

                    let new_frag = Frag{
                        id:v[1].to_string(),
                        counter_id:counter,
                        seq_dict : seqs,
                        qual_dict : quals,
                        positions: positions,
                        first_position : first_position,
                        last_position : last_position,};

                    all_frags.push(new_frag);

                    if first_position < prev_first_pos{
                        panic!("Frags file is not sorted.")
                    }

                    prev_first_pos = first_position;
                    counter+=1
                }

                else{
                    panic!("Not a number found in first column");
                }
            } 
        }
    }
    all_frags
}

pub fn get_genotypes_from_vcf<P>(filename : P) -> FnvHashMap<usize,FnvHashMap<usize,usize>> where P: AsRef <Path>, {
    let mut genotype_dict = FnvHashMap::default();
    if let Ok(lines) = read_lines(filename){
        let mut counter = 1;
        for line in lines{
            if let Ok(l) = line{
                //Skip the first few files of the VCF
                if l.contains('#'){
                    continue
                }

                let v : Vec<&str> = l.split('\t').collect();
                let v : Vec<&str> = v.last().unwrap().split(':').collect();
                let genotypes = v[0];
                let mut split_genotypes =  Vec::new();

                if genotypes.contains('|'){
                    split_genotypes = genotypes.split('|').collect();
                }
                else if genotypes.contains('/'){
                    split_genotypes = genotypes.split('/').collect();
                }
                else{
                    panic!("Genotype column not processed correctly : {}",genotypes);
                }
                let mut genotype_value = FnvHashMap::default();
                for allele in split_genotypes.iter(){
                    let num_allele : usize = allele.parse().unwrap();
                    let val = genotype_value.entry(num_allele).or_insert(0);
                    *val += 1
                }
//                println!("{:?}",genotype_value);
                genotype_dict.insert(counter,genotype_value);
                counter+=1;
            }
        }
    }
    genotype_dict
}
