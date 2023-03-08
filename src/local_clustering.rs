use crate::types_structs::{Frag, HapBlock, GAP_CHAR};
use crate::types_structs::{GenotypeCount, SnpPosition, Genotype};
use ordered_float::OrderedFloat;
//use rand::rng::Rng;
extern crate time;
use crate::utils_frags;
use fxhash::{FxHashSet};


//Return the set of reads for which every read covers at least one position in the interval
//(inclusive)
pub fn find_reads_in_interval<'a>(
    start: SnpPosition,
    end: SnpPosition,
    //position_to_reads : &FxHashMap<usize,FxHashSet<&Frag>>,
    all_frags: &'a Vec<Frag>,
    max_num_reads: usize,
) -> FxHashSet<&'a Frag> {
    let mut final_set = FxHashSet::default();
    //Original method of doing this : This is slower than just iterating thorugh the entire fragment list. We can speed this up by
    //indexing the fragments as well.

    //    for i in (start..end+1).step_by(10){
    //        let empty_set = &FxHashSet::default();
    //        let set1 = position_to_reads.get(&(i)).unwrap_or(empty_set);
    //        for elt in set1.iter() {
    //            final_set.insert(*elt);
    //        }
    //    }

    //Frags are sorted by first position so we can do this.
    for frag in all_frags.iter() {
        if final_set.len() > max_num_reads {
            break;
        }
        if frag.last_position < start {
            continue;
        }
        if frag.first_position > end {
            break;
        }

        //TODO we use this routine in glopp estimate ploidy, don't want circular mappings.
        if frag.last_position - frag.first_position > 10000 {
            continue;
        }

        //Currently we use 1/3 quantile as the length of the block, i.e.
        //end-start. If a mapping is weird and the fragment
        //spans several regions, we ignore the fragment.
        if false {
            //if frag.last_position - frag.first_position > 60 * (end - start) {
            continue;
        }

        final_set.insert(frag);
    }
    final_set
}

//Use the UPEM optimization procedure to optimize the partition by switching around reads to
//optimize UPEM.
//
//partition : the partition
//epislon : read fragment error rate
//genotype_dict : the known genotypes at positions
//polish : if we polish or not
//max_iters : the maximum number of iterations we do.
//div_factor : a normalizing factor for the binomial test to make the sample size smaller.
//use_mec : we can also use MEC score instead of UPEM is desired
pub fn optimize_clustering<'a>(
    partition: Vec<FxHashSet<&'a Frag>>,
    epsilon: f64,
    max_iters: usize,
) -> (f64, Vec<FxHashSet<&'a Frag>>, HapBlock) {
    let mut not_empty = false;
    for part in partition.iter() {
        if part.len() > 0 {
            not_empty = true;
        }
    }

    if not_empty == false {
        let prev_hap_block = utils_frags::hap_block_from_partition(&partition, true);
        return (0.0, partition, prev_hap_block);
    }

    let mut prev_hap_block = utils_frags::hap_block_from_partition(&partition, true);
    let mut set_of_positions = FxHashSet::default();

    for block in prev_hap_block.blocks.iter() {
        for pos in block.keys() {
            set_of_positions.insert(*pos);
        }
    }

    let binom_vec = get_mec_stats_epsilon(&partition, &prev_hap_block, epsilon, true);
    let mut prev_score = binom_vec.iter().map(|x| x.1).sum();
    prev_score *= -1.;

    let mut best_part = partition;

    //Iterate until an iteration yields a lower UPEM score -- return partition corresponding
    //to the best UPEM score.
    for i in 0..max_iters {
        let new_part = opt_iterate(&best_part, &prev_hap_block, epsilon);
        let new_block = utils_frags::hap_block_from_partition(&new_part, true);
        let new_binom_vec = get_mec_stats_epsilon(&new_part, &new_block, epsilon, true);
        let new_score = new_binom_vec.iter().map(|x| x.1).sum::<f64>() * -1.;
        if new_score > prev_score {
            //            log::trace!("Iter {} successful, new {} prev {}", i, new_score, prev_score);
        }

        if new_score > prev_score {
            prev_score = new_score;
            best_part = new_part;
            prev_hap_block = new_block;
        } else {
            log::trace!(
                "Iter {} unsuccessful, new {} prev {}",
                i,
                new_score,
                prev_score
            );
            return (prev_score, best_part, prev_hap_block);
        }
    }

    return (prev_score, best_part, prev_hap_block);
}

//Get a vector of read frequencies and error rates from a partition and its corresponding
//haplotype block.
pub fn get_partition_stats(
    partition: &Vec<FxHashSet<&Frag>>,
    hap_block: &HapBlock,
) -> (Vec<(usize, usize)>, Vec<usize>) {
    let mut binom_vec = Vec::new();
    let mut freq_vec = Vec::new();
    let ploidy = partition.len();
    for i in 0..ploidy {
        let haplo = &hap_block.blocks[i];
        let mut bases = 0;
        let mut errors = 0;
        for frag in partition[i].iter() {
            let (same, diff) = utils_frags::distance_read_haplo(frag, haplo);
            errors += diff;
            bases += same;
        }
        binom_vec.push((bases, errors));
        freq_vec.push(partition[i].len());
    }

    (binom_vec, freq_vec)
}

pub fn get_mec_stats_epsilon_no_phred(
    read_part: &Vec<FxHashSet<&Frag>>,
    epsilon: f64,
) -> Vec<(f64, f64)> {
    let mut binom_vec = vec![];
    let hap_block_no_phred = utils_frags::hap_block_from_partition(read_part, false);
    for hap in hap_block_no_phred.blocks.iter() {
        let mut errors = 0.;
        let mut bases = 0.;
        for seq_dict in hap.values() {
            let mut allele_counts: Vec<(&Genotype, &GenotypeCount)> = seq_dict.iter().collect();
            if allele_counts.is_empty() {
                continue;
            }
            allele_counts.sort_by(|x, y| x.1.cmp(&y.1));
            let allele_counts: Vec<&GenotypeCount> = allele_counts.iter().map(|x| x.1).collect();
            let cons_bases = **allele_counts.last().unwrap();
            bases += *cons_bases;
            for i in 0..allele_counts.len() - 1 {
                errors += allele_counts[i].into_inner();
            }
            if cons_bases <= OrderedFloat(1.) {
                errors += epsilon;
            }
        }
        binom_vec.push((bases, errors));
    }
    return binom_vec;
}

//Include a pental for single coverage alleles.
pub fn get_mec_stats_epsilon(
    _partition: &Vec<FxHashSet<&Frag>>,
    hap_block: &HapBlock,
    epsilon: f64,
    use_gaps: bool,
) -> Vec<(f64, f64)> {
    let mut binom_vec = vec![];
    for hap in hap_block.blocks.iter() {
        let mut errors = 0.;
        let mut bases = 0.;
        for seq_dict in hap.values() {
            let mut allele_counts: Vec<(&Genotype, &GenotypeCount)> = seq_dict.iter().collect();
            if !use_gaps {
                let mut index_to_remove = None;
                for (i, (allele, _count)) in allele_counts.iter_mut().enumerate() {
                    if **allele == GAP_CHAR {
                        index_to_remove = Some(i);
                    }
                }

                if !index_to_remove.is_none() {
                    allele_counts.remove(index_to_remove.unwrap());
                }
            }
            if allele_counts.is_empty() {
                continue;
            }

            allele_counts.sort_by(|x, y| x.1.cmp(&y.1));
            let allele_counts: Vec<&GenotypeCount> = allele_counts.iter().map(|x| x.1).collect();
            let cons_bases = **allele_counts.last().unwrap();
            bases += *cons_bases;
            for i in 0..allele_counts.len() - 1 {
                errors += allele_counts[i].into_inner();
            }
            if cons_bases <= OrderedFloat(1.) {
                errors += epsilon;
            }
        }
        binom_vec.push((bases, errors));
    }
    binom_vec
}

//Return pem score
pub fn get_pem_score(
    binom_vec: &Vec<(usize, usize)>,
    _freq_vec: &Vec<usize>,
    p: f64,
    div_factor: f64,
) -> f64 {
    let mut score = 0.0;
    for stat in binom_vec.iter() {
        let bincdf = utils_frags::stable_binom_cdf_p_rev(stat.0 + stat.1, stat.1, p, div_factor);
        score += bincdf;
    }
    score
}

//Return mec score
pub fn get_mec_score(
    binom_vec: &Vec<(usize, usize)>,
    _freq_vec: &Vec<usize>,
    _p: f64,
    _div_factor: f64,
) -> f64 {
    let mut score = 0;
    for stat in binom_vec.iter() {
        score += stat.1;
    }
    let score_f64 = score as f64;
    score_f64 * -1.0
}

fn opt_iterate<'a>(
    partition: &Vec<FxHashSet<&'a Frag>>,
    hap_block: &HapBlock,
    epsilon: f64,
) -> Vec<FxHashSet<&'a Frag>> {
    let ploidy = partition.len();
    let mut best_moves = Vec::new();

    for i in 0..ploidy {
        if partition[i].len() <= 1 {
            continue;
        }
        for read in partition[i].iter() {
            let haplo_i = &hap_block.blocks[i];
            let (_bases_good_read, errors_read) =
                utils_frags::distance_read_haplo_epsilon_empty(read, haplo_i, epsilon);
            for j in 0..ploidy {
                if j == i {
                    continue;
                }

                //Test out new move
                let haplo_j = &hap_block.blocks[j];
                let (_read_bases_good_movej, read_errors_movej) =
                    utils_frags::distance_read_haplo_epsilon_empty(read, haplo_j, epsilon);

                let diff_score = errors_read - read_errors_movej;
                if diff_score > 0.0 {
                    best_moves.push((diff_score, (i, read, j)));
                    //dbg!(new_score,old_score,read_bases_good_movej,read_errors_movej);
                }
            }
        }
    }

    let mut moved_reads = FxHashSet::default();
    let mut new_part = partition.clone();
    best_moves.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    //    let num_reads: usize = freq_vec.iter().sum();
    //    let mut number_of_moves = num_reads / 10;
    //    if best_moves.len() / 10 < number_of_moves / 5 {
    //        number_of_moves = best_moves.len() / 5;
    //    }
    let mut number_of_moves = best_moves.len() / 10;
    if number_of_moves == 0 && best_moves.len() > 0 {
        number_of_moves = best_moves.len() / 3 + 1;
    }
    //    log::trace!("Number of moves {}", number_of_moves);

    for (mv_num, mv) in best_moves.iter().enumerate() {
        let (i, read, j) = mv.1;
        if moved_reads.contains(read) {
            continue;
        }
        if new_part[i].len() == 1 {
            continue;
        }
        new_part[j].insert(read);
        new_part[i].remove(read);
        moved_reads.insert(read);
        if mv_num > number_of_moves {
            break;
        }
    }

    new_part
}
