use crate::types_structs::{Frag,HapBlock};
use fnv::{FnvHashMap,FnvHashSet};

pub fn polish_using_vcf(genotype_dict : &FnvHashMap<usize,FnvHashMap<usize,usize>>,
                        hap_block : &HapBlock,
                        positions_to_polish : Vec<usize>) -> HapBlock{

    let ploidy = hap_block.blocks.len();
    let mut polished_block = Vec::new();
    for i in 0..ploidy{
        polished_block.push(FnvHashMap::default());
    }
    for pos in positions_to_polish.iter(){

        //Get types of alleles : needs this for polyallelic case
        let mut set_of_variants = FnvHashSet::default();
        let emptydict = FnvHashMap::default();

        for _i in 0..ploidy{
            for var in genotype_dict.get(pos).unwrap_or(&emptydict).keys(){
                set_of_variants.insert(*var);
            }
        }

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
//        println!("{:?}",best_calls_vec);
        let mut vec_pos_maps = Vec::new();
        let mut genotype_counter = FnvHashMap::default();
        for _i in 0..ploidy{
            vec_pos_maps.push(usize::MAX);
        }

        //Very ambiguous calls, dont' do anything
        if best_calls_vec[0][2] == 0{
            println!("Ambiguous position at {}, skipping polishing",pos);
            continue
        }

        for call in best_calls_vec.iter(){
            let ith_hap = call[0] as usize;

            //Already finished calling for this haplotype
            if vec_pos_maps[ith_hap] != usize::MAX{
                continue;
            }

            let allele = call[1] as usize;
            let geno_count = genotype_counter.entry(allele).or_insert(0);
            let geno_count_truth =  genotype_dict.get(pos).unwrap().get(&allele).unwrap_or(&0);
            if (*geno_count_truth as i32) - ((*geno_count) as i32 + 1) >= 0{
                vec_pos_maps[ith_hap] = allele;
                *geno_count += 1;
//                println!("{:?} vec_pos_maps, {:?} geno_count",&vec_pos_maps,genotype_counter);
                if genotype_counter == *genotype_dict.get(pos).unwrap(){
                    break;
                }

            }
        }

        for i in 0..ploidy{
            let var_to_count = polished_block[i].entry(*pos).or_insert(FnvHashMap::default());
            var_to_count.insert(vec_pos_maps[i],1);
        }
//        println!("{:?} vec_pos_maps, {:?} geno_count_truth",vec_pos_maps,genotype_dict.get(pos).unwrap());
    }

    HapBlock{blocks : polished_block}
}

