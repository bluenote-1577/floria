use crate::types_structs::{Frag,HapBlock};
use fnv::{FnvHashMap,FnvHashSet};

pub fn polish_using_vcf(genotype_dict : &FnvHashMap<usize,FnvHashMap<usize,usize>>,
                        hap_block : &HapBlock,
                        positions_to_polish : Vec<usize>) -> HapBlock{

    let ploidy = hap_block.blocks.len();
//    let mut polished_block = Vec::new();
//    for i in 0..ploidy{
//        polished_block.push(FnvHashMap::default());
//    }
    for pos in positions_to_polish.iter(){

        //Get types of alleles : needs this for polyallelic case
        
        let mut set_of_variants = FnvHashSet::default();
        let emptydict = FnvHashMap::default();

        for _i in 0..ploidy{
            for var in genotype_dict.get(pos).unwrap_or(&emptydict).keys(){
                set_of_variants.insert(*var);
            }
        }
        println!("{:?}",set_of_variants);

        //Get a sorted vector of objects [ith_hap,allele_call,error]
        let mut best_calls_vec = Vec::new();
        for i in 0..ploidy{
            let mut total_cov= 0;
            for val in hap_block.blocks[i].get(pos).unwrap_or(&emptydict).values(){
                total_cov += val;
            }
            let total_cov = total_cov as i32;
            for allele in set_of_variants.iter(){
               let allele_cov = *hap_block.blocks[i].get(pos).unwrap_or(&emptydict).get(allele).unwrap_or(&0) as i32;
               //Can make htis error whatever you like
               let call_err = allele_cov* 2 - total_cov;
               let i_i32 = i as i32;
               let allele_i32 = *allele as i32;
               best_calls_vec.push(vec![i_i32,allele_i32, call_err])
            }
        }

        best_calls_vec.sort_by(|a,b| b[2].cmp(&a[2]));
        println!("{:?}",best_calls_vec)
    }

    HapBlock{blocks : Vec::new()}
}

