use crate::local_clustering;
use std::io::LineWriter;
use std::cell::RefCell;
use crate::types_structs::Frag;
use crate::types_structs::HapBlock;
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use permute::permute;
use rayon::prelude::*;
use std::mem;
use std::sync::Mutex;
use std::fs::File;
use std::io::Write;

///This function takes polishes a haplotype block by using genotyping information.
///The algorithm used is simple; we sort the calls for what a haplotype should be
///based on the minimum number of errors to correct for what allele to call. We then
///do the least erroneous calls subject to constraints.
pub fn polish_using_vcf(
    genotype_dict: &FxHashMap<usize, FxHashMap<usize, usize>>,
    hap_block: &HapBlock,
    positions_to_polish: &Vec<usize>,
) -> HapBlock {
    let ploidy = hap_block.blocks.len();
    let mut polished_block = Vec::new();
    for _i in 0..ploidy {
        polished_block.push(FxHashMap::default());
    }
    for pos in positions_to_polish.iter() {
        //Get types of alleles : needs this for polyallelic case
        let mut set_of_variants = FxHashSet::default();
        let emptydict = FxHashMap::default();

        for _i in 0..ploidy {
            for var in genotype_dict.get(pos).unwrap_or(&emptydict).keys() {
                set_of_variants.insert(*var);
            }
        }

        //Get a sorted vector of objects [ith_hap,allele_call,error]
        let mut best_calls_vec = Vec::new();
        for i in 0..ploidy {
            let mut total_cov = 0;
            for val in hap_block.blocks[i].get(pos).unwrap_or(&emptydict).values() {
                total_cov += val;
            }
            let total_cov = total_cov as i32;
            for allele in set_of_variants.iter() {
                let allele_cov = *hap_block.blocks[i]
                    .get(pos)
                    .unwrap_or(&emptydict)
                    .get(allele)
                    .unwrap_or(&0) as i32;
                //Can make htis error whatever you like
                //               let call_err = (allele_cov* 2 - total_cov) as f64;
                let call_err = (allele_cov as f64) / (total_cov as f64 - allele_cov as f64 + 1.0);
                //               dbg!(call_err,allele_cov,total_cov);
                let i_i32 = i as i32;
                let allele_i32 = *allele as i32;
                best_calls_vec.push((i_i32, allele_i32, call_err))
            }
        }

        best_calls_vec.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
        //        println!("{:?}",best_calls_vec);
        let mut vec_pos_maps = Vec::new();
        let mut genotype_counter = FxHashMap::default();
        for _i in 0..ploidy {
            vec_pos_maps.push(usize::MAX);
        }

        //Very ambiguous calls, dont' do anything
        
//        dbg!(&best_calls_vec);
        if best_calls_vec[0].2 == 0.0 {
            //println!("Ambiguous position at {}, skipping polishing", pos);
            continue;
        }

        for call in best_calls_vec.iter() {
            let ith_hap = call.0 as usize;

            //Already finished calling for this haplotype
            if vec_pos_maps[ith_hap] != usize::MAX {
                continue;
            }

            let allele = call.1 as usize;
            let geno_count = genotype_counter.entry(allele).or_insert(0);
            let geno_count_truth = genotype_dict.get(pos).unwrap().get(&allele).unwrap_or(&0);
            if (*geno_count_truth as i32) - ((*geno_count) as i32 + 1) >= 0 {
                vec_pos_maps[ith_hap] = allele;
                *geno_count += 1;
                //                println!("{:?} vec_pos_maps, {:?} geno_count",&vec_pos_maps,genotype_counter);
                if genotype_counter == *genotype_dict.get(pos).unwrap() {
                    break;
                }
            }
        }

        for i in 0..ploidy {
            let var_to_count = polished_block[i]
                .entry(*pos)
                .or_insert(FxHashMap::default());
            var_to_count.insert(vec_pos_maps[i], 1);
        }
        //        println!("{:?} vec_pos_maps, {:?} geno_count_truth",vec_pos_maps,genotype_dict.get(pos).unwrap());
    }

    HapBlock {
        blocks: polished_block,
    }
}


//Link two partitions by best MEC score permutation. This doesn't help much
//on the simulated datasets. 
fn get_best_perms_mec(part1: &Vec<FxHashSet<&Frag>>, part2: &Vec<FxHashSet<&Frag>>) -> Vec<Vec<usize>> {
    let ploidy = part1.len();
    let mut best_perms = Vec::new();
    let rangevec: Vec<usize> = (0..ploidy).collect();
    let perms = permute(rangevec);
    for perm in perms {
        let mut test_part = Vec::new();
        for i in 0..ploidy {
            let mut new_set = FxHashSet::default();
            let j = perm[i];

            let set1 = &part1[i];
            let set2 = &part2[j];

            for read in set1.iter(){
                new_set.insert(*read);
            }
            for read in set2.iter(){
                new_set.insert(*read);
            }

            test_part.push(new_set);
        }
        let score = get_mec_from_part(&test_part);
        best_perms.push((-1*score,perm));
    }
    best_perms.sort_by(|a,b| b.0.cmp(&a.0));
    let best_perms : Vec<Vec<usize>> = best_perms.into_iter().map(|x| x.1).collect();
    best_perms
}

//Link two partitions by which permutation gives the most amount ofintersections between the sets.
fn get_best_perms(part1: &Vec<FxHashSet<&Frag>>, part2: &Vec<FxHashSet<&Frag>>) -> Vec<Vec<usize>> {
    let ploidy = part1.len();
    let mut best_perms = Vec::new();
    let rangevec: Vec<usize> = (0..ploidy).collect();
    let perms = permute(rangevec);
    for perm in perms {
        let mut score = 0;
        for i in 0..ploidy {
            let j = perm[i];

            let set1 = &part1[i];
            let set2 = &part2[j];

            let int_vec: Vec<&Frag> = set1.intersection(set2).copied().collect();
            score += int_vec.len();
        }
        best_perms.push((score,perm));
    }
    best_perms.sort_by(|a,b| b.0.cmp(&a.0));
    let best_perms : Vec<Vec<usize>> = best_perms.into_iter().map(|x| x.1).collect();
    best_perms
}

///Link all partitions in a vector using get_best_perm.
pub fn link_blocks<'a>(all_parts: &Vec<Vec<FxHashSet<&'a Frag>>>) -> Vec<FxHashSet<&'a Frag>> {
    //Multithreaded version -- not super useful unless ploidy > 6. Might as well though.
    let mut final_part = all_parts[0].clone();
    let perm_vector: Mutex<Vec<(usize, Vec<usize>)>> = Mutex::new(vec![]);
    let ploidy = final_part.len();
    let mut current_perm: Vec<usize> = (0..ploidy).collect();

    //Get the best permutation for each consecutive pair.
    (1..all_parts.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
//            let best_perm = get_best_perms_mec(&all_parts[i - 1], &all_parts[i])[0].clone();
            let best_perm = get_best_perms(&all_parts[i - 1], &all_parts[i])[0].clone();
            let mut locked_perm = perm_vector.lock().unwrap();
            locked_perm.push((i, best_perm));
        });

    let mut perm_vector = perm_vector.lock().unwrap().to_vec();
    perm_vector.sort_by(|a, b| a.0.cmp(&b.0));
    let perm_vector: Vec<Vec<usize>> = perm_vector.into_iter().map(|x| x.1).collect();

    //We have to do compose permutations here to keep track of
    //what state the permutations between the pairs needs to be in.
    let mut all_used_reads = FxHashSet::default();
    let mut all_used_reads_vec = vec!();
    for reads in all_parts[0].iter(){
        for read in reads.iter(){
            all_used_reads.insert(read);
            all_used_reads_vec.push(read);
        }
    }

    for (i, perm) in perm_vector.iter().enumerate() {
        let part_to_link = &all_parts[i + 1];
        for j in 0..ploidy {
            let set1 = &mut final_part[j];
            let set2 = &part_to_link[perm[current_perm[j]]];

            for read in set2.into_iter() {
//                if !all_used_reads.contains(read){
                    set1.insert(read);
                    all_used_reads.insert(read);
                    all_used_reads_vec.push(read);
//                }
            }
        }
        for j in 0..ploidy {
            current_perm[j] = perm[current_perm[j]];
        }
    }

    final_part
}



fn clone_block_range (hap_block : &HapBlock, start : usize,end : usize) -> HapBlock{
    let mut new_hap_block = HapBlock{blocks : Vec::new()};
    let ploidy = hap_block.blocks.len();
    for _i in 0..ploidy{
        new_hap_block.blocks.push(FxHashMap::default());
    }
    for pos in start..end+1{
        for j in 0..ploidy{
            let hap = &hap_block.blocks[j];
            if hap.contains_key(&pos){
                new_hap_block.blocks[j].insert(pos,hap.get(&pos).unwrap().clone());
            }
        }
    }

    return new_hap_block;
}

fn get_mec_from_part(part : &Vec<FxHashSet<&Frag>>) -> i32{

    let block = utils_frags::hap_block_from_partition(part);
    let mut mec_score = 0;
    for hap in block.blocks.iter(){
        for (_pos,site_map) in hap{
            let mut total_sites = 0;
            let mut sites_max_var = 0;
            for value in site_map.values(){
                if *value > sites_max_var{
                    sites_max_var = *value;
                }
                total_sites += *value;
            }
            mec_score += total_sites - sites_max_var;
        }
    }

    mec_score as i32
}

fn get_mec_positions_hap(new_reads : &FxHashSet<&Frag>, hap_block : &mut HapBlock, pos_sort_vec : &Vec<usize>, l : usize) -> usize {
    
    let hap = &mut hap_block.blocks[l];

    //Add reads to hap_block
    for read in new_reads.iter(){
        for pos in read.positions.iter(){
            let var_at_pos = read.seq_dict.get(pos).unwrap();
            let sites = hap.entry(*pos).or_insert(FxHashMap::default());
            let site_counter = sites.entry(*var_at_pos).or_insert(0);
            *site_counter += 1;
        }
    }

    //Compute new MEC for the new reads
    let mut err = 0;
    for pos in pos_sort_vec[0]..pos_sort_vec[pos_sort_vec.len()-1]{
        if hap.contains_key(&pos){
            let mut total_sites = 0;
            let mut sites_max_var = 0;
            let site_map = hap.get(&pos).unwrap();
            for value in site_map.values(){
                if *value > sites_max_var{
                    sites_max_var = *value;
                }
                total_sites += *value;
            }
            err += total_sites - sites_max_var;
        }
    }

    //Remove reads from hap_block
    for read in new_reads.iter(){
        for pos in read.positions.iter(){
            let var_at_pos = read.seq_dict.get(pos).unwrap();
            let sites = hap.entry(*pos).or_insert(FxHashMap::default());
            let site_counter = sites.entry(*var_at_pos).or_insert(0);
            *site_counter -= 1;
        }
    }

    return err;
}

fn get_iqr(all_scores: &Vec<f64>, factor: f64) -> f64 {
    let mut score_copy = all_scores.clone();
    score_copy.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i_25 = score_copy.len() / 4;
    let i_75 = score_copy.len() * 3 / 4;
    let q25 = score_copy[i_25];
    let q75 = score_copy[i_75];

    let iqr = q75 - q25;
    return q25 - factor * iqr;
}

fn fill_left_block<'a>(left_block: &mut Vec<FxHashSet<&'a Frag>>, reads_interval: Vec<&'a Frag>, epsilon : f64) {
    for frag in reads_interval {
        let mut read_already_used = false;

        for read_set in left_block.iter() {
            if read_set.contains(frag) {
                read_already_used = true;
                break;
            }
        }

        if read_already_used {
            continue;
        }

        //Read hasn't been used : try to fill in.

        let mut max_distances = Vec::new();
        for (i, read_set) in left_block.iter().enumerate() {
            let mut max_dist = 0.0;
            for frag_in_set in read_set.iter() {
                let (same,mec_dist) = utils_frags::distance(frag_in_set, frag);
                let binomial_dist = -1.0 * local_clustering::stable_binom_cdf_p_rev(
                    (same + mec_dist) as usize,
                    mec_dist as usize,
                    2.0 * epsilon * (1.0-epsilon),
                    100.0
                );
                    
                let dist = binomial_dist;
//                let dist = mec_dist as f64;
                if dist > max_dist {
                    max_dist = dist ;
                }
            }
            max_distances.push((i, max_dist));
        }

        max_distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        //        dbg!(&max_distances);
        left_block[max_distances[0].0].insert(frag);
    }
}

pub fn replace_with_filled_blocks<'a>(
    all_scores: &Vec<f64>,
    mut all_parts: Vec<Vec<FxHashSet<&'a Frag>>>,
    factor: f64,
    length_of_block: usize,
    all_frags: &'a Vec<Frag>,
    epsilon : f64
) -> Vec<Vec<FxHashSet<&'a Frag>>> {
    let mut corrected_vec = Vec::new();
    let outlier_score = get_iqr(all_scores, factor);

    //Assume the leftmost block is good. Not a great assumption but
    //otherwise the algorithm would be a bit more painful.
    let first_block = mem::replace(&mut all_parts[0], Vec::new());
    corrected_vec.push(first_block);

    for i in 1..all_parts.len() {
        let score = all_scores[i];

        if score > outlier_score {
            let block_to_push = mem::replace(&mut all_parts[i], Vec::new());
            corrected_vec.push(block_to_push);
            continue;
        }

        println!("Filling in block {}", i);

        //Bad block, fill in from left
        let mut vec_reads_interval: Vec<&Frag> = local_clustering::find_reads_in_interval(
            i * length_of_block,
            (i + 1) * length_of_block,
            all_frags,
        )
        .into_iter()
        .collect();
        vec_reads_interval.sort_by(|a, b| a.first_position.cmp(&b.first_position));
        fill_left_block(corrected_vec.iter_mut().last().unwrap(), vec_reads_interval,epsilon);
    }
    corrected_vec
}

//Modified beam search ... not working right now and incomplete. Maybe
//modify for usage in the future. 
pub fn link_blocks_heur<'a>(all_parts: &Vec<Vec<FxHashSet<&'a Frag>>>, num_sol : usize) -> Vec<FxHashSet<&'a Frag>> {
    let perm_vector: Mutex<Vec<(usize, Vec<Vec<usize>>)>> = Mutex::new(vec![]);
    let ploidy = all_parts[0].len();

    //Get the best permutation for each consecutive pair.
    (1..all_parts.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
            let mut locked_perm = perm_vector.lock().unwrap();
            locked_perm.push((i, get_best_perms(&all_parts[i - 1], &all_parts[i])));
        });

    let mut perm_vector = perm_vector.lock().unwrap().to_vec();
    perm_vector.sort_by(|a, b| a.0.cmp(&b.0));
    let perm_vector: Vec<Vec<Vec<usize>>> = perm_vector.into_iter().map(|x| x.1).collect();

    //Start with the first k different links.
    let first_perms = &perm_vector[0];
    let mut mec_perm_block_part = Vec::new();

    for i in 0..num_sol{
        let perm = &first_perms[i];
        let mut new_part = Vec::new();

        for j in 0..ploidy{
            let set1 = &all_parts[0][j];
            let set2 = &all_parts[1][perm[j]];
            let mut union_set = FxHashSet::default();
            
            for read in set1.iter(){
                union_set.insert(*read);
            }

            for read in set2.iter(){
                union_set.insert(*read);
            }
            new_part.push(union_set);
        }

        let block = utils_frags::hap_block_from_partition(&new_part);
        //TODO Polish stuff
        let (binom_vec,_freq_vec)= local_clustering::get_partition_stats(&new_part,&block);

        let mut errors = 0;
        for tup in binom_vec.iter(){
            errors += tup.1;
        }

        let wrapped_block = RefCell::new(block);

        mec_perm_block_part.push((errors,perm.clone(),wrapped_block,new_part));
    }

    for i in 1..perm_vector.len(){
        let mut new_mec_perm_complete = Vec::new();
        let mut new_mec_perm = Vec::new();
        let perms = &perm_vector[i];
        let mut leftmost_pos = usize::MAX;
        let mut rightmost_pos = usize::MIN;
        //Iterate over each previous candidate solution
        for j in 0..num_sol{
            let (errors,curr_perm,wrapped_block,part) = & mec_perm_block_part[j];
            //Iterate over each permutation to create k^2 different solutions
            for k in 0..num_sol{
                let perm = &perms[k];
                let mut new_errors = *errors;
                let mut vec_new_part = Vec::new();

                //Iterate over each haplotype
                for l in 0..ploidy{
                    
                    let mut list_of_pos_inserts = Vec::new();
                    let set1 = &part[l];
                    let set2 = &all_parts[i+1][perm[curr_perm[l]]];
                    let mut new_reads_set = FxHashSet::default();
                    let mut union_set = FxHashSet::default();

                    for read in set2.iter(){
                        list_of_pos_inserts.push(read.first_position);
                        list_of_pos_inserts.push(read.last_position);

                        if read.first_position < leftmost_pos{
                            leftmost_pos = read.first_position;
                        }

                        if read.last_position > rightmost_pos{
                            rightmost_pos = read.last_position;
                        }

                        let mut contains_read = false;
                        for p in 0..ploidy{
                            if (part[p]).contains(read){
                               contains_read = true 
                            }
                        }

                        if !contains_read{
                            new_reads_set.insert(*read);
                            union_set.insert(*read);
                        }

                    }

                    for read in set1.iter(){
                        union_set.insert(*read);
                    }


                    let emptyset :FxHashSet<&Frag>= FxHashSet::default();
                    list_of_pos_inserts.sort();
                    let mec_before_positions = get_mec_positions_hap(&emptyset, &mut wrapped_block.borrow_mut(), &list_of_pos_inserts, l);
                    let mec_after_positions = get_mec_positions_hap(&new_reads_set, &mut wrapped_block.borrow_mut(), &list_of_pos_inserts, l);
                    vec_new_part.push(union_set);
//
                    new_errors += mec_after_positions - mec_before_positions;
                }
                let mut new_curr_perm = vec![];
                for l in 0..ploidy{
                    new_curr_perm.push(perm[curr_perm[l]]);
                }
                new_mec_perm.push((new_errors,new_curr_perm, wrapped_block ,vec_new_part));
            }
        }

        //Now have a list of k^2 candidate solutions, sort and get the best ones.
        new_mec_perm.sort_by(|a,b| a.0.cmp(&b.0));

        for j in 0..num_sol{
            let (new_errors,new_curr_perm,wrapped_block,new_part) = &mut new_mec_perm[j];
            let swap_part = mem::replace(new_part, Vec::new());
            new_mec_perm_complete.push((*new_errors,new_curr_perm.clone(), RefCell::new(clone_block_range(&wrapped_block.borrow(),leftmost_pos,rightmost_pos)),swap_part));
        }

        mec_perm_block_part = new_mec_perm_complete;
    }

    mec_perm_block_part.into_iter().nth(0).unwrap().3
}

pub fn remove_duplicate_reads(
    part : &mut Vec<FxHashSet<&Frag>>,
    all_frags : &Vec<Frag>,
    block : &HapBlock){

    for frag in all_frags.iter(){
        let mut parts_containing_frag = Vec::new();
        for i in 0..part.len(){
            if part[i].contains(frag){
                parts_containing_frag.push(i);
            }
        }

        if parts_containing_frag.len() > 1{
            let mut mec_hap_vec = Vec::new();
            for i in parts_containing_frag{
                let hap = &block.blocks[i];
                let (_same,diff) = utils_frags::distance_read_haplo(frag,hap);
                mec_hap_vec.push((diff,i));
            }

            mec_hap_vec.sort_by(|a,b| a.0.cmp(&b.0));
            for j in 1..mec_hap_vec.len(){
                //if mec_hap_vec[i].0 > mec_hap_vec[0].0{
                let i = mec_hap_vec[j].1;
                part[i].remove(frag);
                //}
            }

        }
    }

}


pub fn map_reads_against_hap_errors(
    part : &Vec<FxHashSet<&Frag>>,
    block : &HapBlock,
    epsilon : f64){

    let mut list_of_errors = Vec::new();
    let ploidy = part.len();
    for i in 0..ploidy{
        let reads = &part[i];
        let hap = &block.blocks[i];
        for read in reads{
            let length = read.last_position - read.first_position + 1;

            let start_left = read.first_position;
            let end_left = read.first_position + length/2;
            let (same_left,diff_left) = utils_frags::distance_read_haplo_range(read,hap,start_left,end_left);

            let start_right = read.first_position + length/2+1;
            let end_right= read.last_position;
            let (same_right,diff_right) = utils_frags::distance_read_haplo_range(read,hap,start_right,end_right);

            let mut error_est = ((diff_left + diff_right) as f64) / ((same_right + diff_left + same_left + diff_right) as f64);
            if error_est == 0.0{
                error_est = epsilon;
            }

            let left_error_val = local_clustering::stable_binom_cdf_p_rev(same_left + diff_left, diff_left, error_est,1.0);
            let right_error_val = local_clustering::stable_binom_cdf_p_rev(same_right + diff_right, diff_right, error_est,1.0);

            let mut min_error_val = left_error_val;
            if left_error_val > right_error_val{
                min_error_val = right_error_val;
            }

            if !min_error_val.is_finite(){
                continue;
            }

            list_of_errors.push((min_error_val,&read.id,read.first_position,read.last_position,i,diff_left,diff_right));
        }
    }

    list_of_errors.sort_by(|a,b| a.0.partial_cmp(&b.0).unwrap()); 

    let file = File::create("possible_breaks.txt").expect("Can't create file");
    let mut file = LineWriter::new(file);

    for i in 0..list_of_errors.len(){
        write!(file,"{},{},{},{}\n",list_of_errors[i].0,list_of_errors[i].1,list_of_errors[i].3,list_of_errors[i].4).unwrap();
    }
}

pub fn link_blocks_greedy<'a>(all_parts: &Vec<Vec<FxHashSet<&'a Frag>>>,
    all_frags: &'a Vec<Frag>,
    ) -> Vec<FxHashSet<&'a Frag>> {
    //Multithreaded version -- not super useful unless ploidy > 6. Might as well though.
    let mut final_part = all_parts[0].clone();
    let ploidy = final_part.len();
    let mut part_size_distribution_sum = vec!();
    let rangevec: Vec<usize> = (0..ploidy).collect();
    let perms = permute(rangevec);


    let mut all_used_reads = FxHashSet::default();
    for reads in all_parts[0].iter(){
        for read in reads.iter(){
            all_used_reads.insert(read);
        }
    }

    //Find best permutation
    for i in 1..all_parts.len(){
        let part_to_link = &all_parts[i];
        let mut best_total_intersect = 0;
        let mut best_perm = &perms[0];
        for perm in &perms{
            let mut total_intersect = 0;
            for j in 0..ploidy{

                let set1 = &final_part[j];
                let set2 = &part_to_link[perm[j]];

                for read in set2.iter(){
                    if set1.contains(read){
                        total_intersect += 1
                    }
                }
            }

            if total_intersect > best_total_intersect{
                best_total_intersect = total_intersect;
                best_perm = perm;
            }
        }

        for j in 0..ploidy{
            let set1 = &mut final_part[j];
            let set2 = &part_to_link[best_perm[j]];

            for read in set2.iter(){
//                if !all_used_reads.contains(read){
                   set1.insert(read); 
                   all_used_reads.insert(read);
//                }
            }
        }

        //Remove duplicated reads which may occur because of partition joining.
        //TODO need to benchmark which way to remove duplicate reads. 
        //let hap_block_intermediate = utils_frags::hap_block_from_partition(&final_part);
        //remove_duplicate_reads(&mut final_part, &all_frags, &hap_block_intermediate);
        //TEST
        let mut countvec1 = vec!();
        let mut countvec2 = vec!();
        for j in 0..ploidy{
            part_size_distribution_sum.push(vec!());
            countvec1.push(final_part[j].len());
            countvec2.push(part_to_link[best_perm[j]].len());
        }
        dbg!(&countvec1,&countvec2);
        countvec2.sort();
        for j in 0..ploidy{
            part_size_distribution_sum[j].push(countvec2[j]);
        }
    }

    for j in 0..ploidy{
        part_size_distribution_sum[j].sort();
        let dist = &part_size_distribution_sum[j];
        dbg!(dist[(dist.len() as f64 * 0.5) as usize]);
    }
    final_part
}








