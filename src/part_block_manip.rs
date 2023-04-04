use crate::constants;
use crate::types_structs::{Frag, VcfProfile};
use crate::types_structs::{Options, SnpPosition};
use crate::utils_frags;
use disjoint_sets::UnionFind;
use fxhash::{FxHashMap, FxHashSet};
use rust_lapper::{Interval, Lapper};
use std::fs::File;
use std::io::Write;
use std::mem;
use std::time::Instant;

fn overlap_percent(x1: SnpPosition, x2: SnpPosition, y1: SnpPosition, y2: SnpPosition) -> f64 {
    let intersect = SnpPosition::max(SnpPosition::min(x2 - y1 + 1, y2 - x1 + 1), 0);
    let min_length = x2 - x1 + 1;
    let p = intersect as f64 / min_length as f64;
    if p > 1. {
        return 1.;
    //        println!("{},{},{},{}",x1,x2,y1,y2);
    //        panic!();
    } else {
        return p;
    }
}

//TODO not sure if this is needed. Will implement if needed.
fn separate_broken_haplogroups<'a>(
    all_joined_path_parts: &mut Vec<FxHashSet<&'a Frag>>,
    snp_range_parts_vec: &mut Vec<(SnpPosition, SnpPosition)>,
) {
    let mut all_breaks = vec![];
    for i in 0..snp_range_parts_vec.len() {
        let part = &all_joined_path_parts[i];
        let mut vec_of_frags: Vec<&&Frag> = part.into_iter().collect();
        vec_of_frags.sort_by(|x, y| x.first_position.cmp(&y.first_position));
        let mut current_lastest_pos = 0;
        let mut breaks = vec![];
        for frag in vec_of_frags {
            if current_lastest_pos != 0 && frag.first_position > current_lastest_pos {
                if current_lastest_pos >= snp_range_parts_vec[i].0
                    && current_lastest_pos < snp_range_parts_vec[i].1
                {
                    breaks.push(current_lastest_pos);
                }
            }
            if frag.last_position > current_lastest_pos {
                current_lastest_pos = frag.last_position;
            }
        }
        if !breaks.is_empty() {
            all_breaks.push((i, breaks));
        }
    }

    let mut new_parts = vec![];
    let mut new_ranges = vec![];
    for break_info in all_breaks.iter() {
        let mut spot_index = 0;
        let break_spots = &break_info.1;
        let mut break_start = snp_range_parts_vec[break_info.0].0;
        let part = &all_joined_path_parts[break_info.0];
        let mut vec_of_frags: Vec<&&Frag> = part.iter().collect();
        vec_of_frags.sort_by(|x, y| x.first_position.cmp(&y.first_position));

        let mut end_spot = break_spots[spot_index];
        let mut new_part = FxHashSet::default();

        for frag in vec_of_frags.iter() {
            if frag.last_position <= end_spot {
                new_part.insert(**frag);
            } else {
                let new_snp_range = (break_start, end_spot);
                new_parts.push(new_part);
                new_ranges.push(new_snp_range);

                break_start = end_spot + 1;
                spot_index += 1;
                if spot_index != break_spots.len() {
                    end_spot = break_spots[spot_index];
                } else {
                    end_spot = SnpPosition::MAX;
                }
                new_part = FxHashSet::default();
            }
        }
        let new_snp_range = (break_start, snp_range_parts_vec[break_info.0].1);
        new_parts.push(new_part);
        new_ranges.push(new_snp_range);
    }

    for break_info in all_breaks {
        all_joined_path_parts[break_info.0].clear();
    }
    for i in 0..new_parts.len() {
        all_joined_path_parts.push(mem::take(&mut new_parts[i]));
        snp_range_parts_vec.push(new_ranges[i]);
    }
}
fn merge_overlapping_haplogroups<'a>(
    all_joined_path_parts: &mut Vec<FxHashSet<&'a Frag>>,
    snp_range_parts_vec: &mut Vec<(SnpPosition, SnpPosition)>,
    epsilon: f64,
) {
    let all_parts_block = utils_frags::hap_block_from_partition(&all_joined_path_parts, true);
    let (all_overlaps, _) = find_overlapping_blocks(
        all_joined_path_parts,
        constants::MERGE_CUTOFF,
        snp_range_parts_vec,
    );
    let mut merge_redirect = UnionFind::new(snp_range_parts_vec.len());
    for (index, overlap_vec) in all_overlaps.iter() {
        let mut potential_merges = vec![];
        for inter in overlap_vec.iter() {
            let orig_snp_range = &snp_range_parts_vec[*index];
            let inter_snp_range = &snp_range_parts_vec[inter.val];
            let check_range = (
                SnpPosition::min(orig_snp_range.0, inter_snp_range.0),
                SnpPosition::max(orig_snp_range.1, inter_snp_range.1),
            );
            let (same, diff) = utils_frags::distance_between_haplotypes(
                &all_parts_block.blocks[*index],
                &all_parts_block.blocks[inter.val],
                &check_range,
            );
            if (diff / (same + diff)) < epsilon {
                log::trace!("potential_merge {} {} {} {}", diff, same, index, inter.val);
                potential_merges.push((index, inter.val, check_range.0, check_range.1, same, diff));
            }
        }
        if !potential_merges.is_empty() {
            let best_merge = *potential_merges
                .iter()
                .max_by(|x, y| (x.3 - x.2).cmp(&(y.3 - y.2)))
                .unwrap();
            merge_redirect.union(*best_merge.0, best_merge.1);
        }
    }

    let mut disjoint_set_to_merge = FxHashMap::default();
    for index in 0..snp_range_parts_vec.len() {
        let rep = merge_redirect.find(index);
        let set = disjoint_set_to_merge
            .entry(rep)
            .or_insert(FxHashSet::default());
        set.insert(index);
    }

    for (rep, set) in disjoint_set_to_merge {
        if set.len() <= 1 {
            continue;
        }
        let mut all_range = (SnpPosition::MAX, SnpPosition::MIN);
        for index in set {
            log::trace!("MERGING {} {}", rep, index);
            if all_range.0 > snp_range_parts_vec[index].0 {
                all_range.0 = snp_range_parts_vec[index].0;
            }
            if all_range.1 < snp_range_parts_vec[index].1 {
                all_range.1 = snp_range_parts_vec[index].1;
            }
            if index != rep {
                let old_set = mem::take(&mut all_joined_path_parts[index]);
                let set_parent = &mut all_joined_path_parts[rep];

                for read in old_set {
                    set_parent.insert(read);
                }
            }
        }
        snp_range_parts_vec[rep] = all_range;
    }
}

pub fn process_reads_for_final_parts<'a>(
    mut all_joined_path_parts: Vec<FxHashSet<&'a Frag>>,
    short_frags: &'a Vec<Frag>,
    mut snp_range_parts_vec: Vec<(SnpPosition, SnpPosition)>,
    options: &Options,
    _snp_to_genome_pos: &'a Vec<usize>,
) -> (Vec<FxHashSet<&'a Frag>>, Vec<(SnpPosition, SnpPosition)>) {
    let reassign_short = options.reassign_short;
    let epsilon = options.epsilon;
    //The interval method failed, but may be useful in the future.
    let mut all_parts_block = utils_frags::hap_block_from_partition(&all_joined_path_parts, true);
    let mut read_to_parts_map = FxHashMap::default();
    for (i, set) in all_joined_path_parts.iter().enumerate() {
        for frag in set.iter() {
            let corresponding_partitions = read_to_parts_map
                .entry(*frag)
                .or_insert(FxHashSet::default());
            corresponding_partitions.insert(i);
        }
    }

    for (frag, part_ids) in read_to_parts_map.iter() {
        for id in part_ids.iter() {
            all_joined_path_parts[*id].remove(frag);
            utils_frags::remove_read_from_block(&mut all_parts_block, frag, *id);
        }
    }
    //Only have the long-read fragment be in the best candidate haplotig
    //instead of appearing in multiple
    for (frag, part_ids) in read_to_parts_map {
        let mut diff_part_vec = vec![];
        for id in part_ids.iter() {
            let block_with_id = &all_parts_block.blocks[*id];
            let (same, diff) =
                utils_frags::distance_read_haplo_epsilon_empty(frag, block_with_id, epsilon);
            //This didn't work well
            //            diff_part_vec.push(((diff + 1.) / (same + 1.), id));
            diff_part_vec.push(((diff + 1.), id, (same)));
        }
        let best_part = diff_part_vec
            .iter()
            .min_by(|x, y| x.partial_cmp(&y).unwrap())
            .unwrap()
            .1;
        //        log::trace!("Frag {}, best partitions {:?}", &frag.id, &diff_part_vec);
        all_joined_path_parts[*best_part].insert(frag);
        utils_frags::add_read_to_block(&mut all_parts_block, frag, *best_part);
    }

    if constants::MERGE_SIMILAR_HAPLOGROUPS {
        merge_overlapping_haplogroups(
            &mut all_joined_path_parts,
            &mut snp_range_parts_vec,
            epsilon,
        );
    }
    if constants::SEPARATE_BROKEN_HAPLOGROUPS {
        separate_broken_haplogroups(&mut all_joined_path_parts, &mut snp_range_parts_vec);
    }

    if reassign_short {
        let start_t = Instant::now();
        //Assign short-read fragments to multiple best haplotigs
        for frag in short_frags.iter() {
            let mut best_candidate_block_map = FxHashMap::default();
            for (i, block) in all_parts_block.blocks.iter().enumerate() {
                let mut intersection = false;
                let (a, b) = snp_range_parts_vec[i];
                if a <= frag.first_position && frag.first_position <= b {
                    intersection = true;
                } else if a <= frag.last_position && frag.last_position <= b {
                    intersection = true;
                }
                if intersection {
                    let (same, diff) =
                        utils_frags::distance_read_haplo_epsilon_empty(frag, block, epsilon);
                    let score = ((diff * 10. + 1.) as usize, (same * 10. + 1.) as usize);
                    let candidate_blocks = best_candidate_block_map.entry(score).or_insert(vec![]);
                    candidate_blocks.push(i);
                }
            }
            let best_key = best_candidate_block_map.keys().min_by(|x, y| {
                (x.0 as f64 / x.1 as f64)
                    .partial_cmp(&(y.0 as f64 / y.1 as f64))
                    .unwrap()
            });
            //        dbg!(best_candidate_block_map.keys(), best_key);
            if !best_key.is_none() {
                for i in best_candidate_block_map[&best_key.unwrap()].iter() {
                    all_joined_path_parts[*i].insert(frag);
                }
            }
        }

        log::info!("Time taken for reassign {:?}", Instant::now() - start_t);
    }

    let (sorted_parts, sorted_ranges) = sort_parts(all_joined_path_parts, snp_range_parts_vec);
    return (sorted_parts, sorted_ranges);
}

fn sort_parts<'a>(
    all_joined_path_parts: Vec<FxHashSet<&'a Frag>>,
    snp_range_parts_vec: Vec<(SnpPosition, SnpPosition)>,
) -> (Vec<FxHashSet<&'a Frag>>, Vec<(SnpPosition, SnpPosition)>) {
    let mut zipped: Vec<(FxHashSet<&'a Frag>, (SnpPosition, SnpPosition))> = all_joined_path_parts
        .into_iter()
        .zip(snp_range_parts_vec.into_iter())
        .collect();
    zipped.sort_by(|x, y| x.1.cmp(&y.1));
    let new_range = zipped.iter().map(|x| x.1).collect();
    let new_parts = zipped.into_iter().map(|x| x.0).collect();
    return (new_parts, new_range);
}

pub fn bin_haplogroups<'a>(
    parts: &Vec<FxHashSet<&'a Frag>>,
    snp_endpoints: &Vec<(SnpPosition, SnpPosition)>,
    cov_of_haplogroups: &Vec<Option<f64>>,
    vcf_profile: &VcfProfile,
    contig: &str,
    block_len: usize,
) -> (Vec<(SnpPosition, SnpPosition)>, Vec<FxHashSet<&'a Frag>>) {
    use statrs::distribution::{Discrete, Poisson};

    fn overlap(x1: usize, x2: usize, y1: usize, y2: usize) -> bool {
        if x2 > y1 && x2 < y2 {
            return true;
        } else if y2 > x1 && y2 < x2 {
            return true;
        } else if x1 >= y1 && x2 <= y2 {
            return true;
        } else if x1 <= y1 && x2 >= y2 {
            return true;
        } else {
            return false;
        }
    }

    fn close_enough(x1: usize, x2: usize, y1: usize, y2: usize, block_len: usize) -> bool {
        if (x2 as i64 - y1 as i64).abs() < 2 * block_len as i64
            || (y2 as i64 - x1 as i64).abs() < 2 * block_len as i64
        {
            return true;
        } else {
            return false;
        }
    }

    fn dist(
        x: &Vec<(usize, usize, f64, usize)>,
        y: &Vec<(usize, usize, f64, usize)>,
        block_len: usize,
    ) -> f64 {
        let mut compat_ol = true;
        let mut compat_ce = false;
        for hap1 in x.iter() {
            for hap2 in y.iter() {
                let ol = overlap(hap1.0, hap1.1, hap2.0, hap2.1);
                let ce = close_enough(hap1.0, hap1.1, hap2.0, hap2.1, block_len);
                if ce {
                    compat_ce = true;
                }
                if ol {
                    compat_ol = false;
                    break;
                }
            }
        }

        if !compat_ol || !compat_ce {
            return f64::MAX;
        } else {
            let cov_x = x.iter().map(|x| x.2).sum::<f64>() / x.len() as f64;
            let cov_y = y.iter().map(|x| x.2).sum::<f64>() / y.len() as f64;
            let poix = Poisson::new(cov_x).unwrap();
            let poiy = Poisson::new(cov_y).unwrap();
            let dist = poix.pmf(cov_y as u64) + poiy.pmf(cov_x as u64);
            //            dbg!(x,y,cov_x,cov_y,dist);
            return -1. * (dist / 2.).ln();
        }
    }

    assert!(parts.len() == cov_of_haplogroups.len());
    assert!(parts.len() == snp_endpoints.len());
    let snp_to_gn_pos = &vcf_profile.vcf_snp_pos_to_gn_pos_map[contig];

    let mut clusters = vec![];
    let mut none_clusters = vec![];
    for i in 0..snp_endpoints.len() {
        let left_gn = snp_to_gn_pos[&(snp_endpoints[i].0)] as usize;
        let right_gn = snp_to_gn_pos[&(snp_endpoints[i].1)] as usize;
        let cov = cov_of_haplogroups[i];
        if !cov.is_none() {
            clusters.push(vec![(left_gn, right_gn, cov.unwrap(), i)]);
        } else {
            none_clusters.push(i);
        }
    }

    clusters.sort_by(|x, y| x[0].0.cmp(&y[0].0));
    let cutoff_dist = -(0.01f64.ln());
    loop {
        let mut best_moves = vec![];
        for i in 0..clusters.len() {
            let mut best_moves_i = vec![];
            let h = 100;
            let max_move = usize::min(clusters.len(), i + h);
            let min_move;

            if h > i {
                min_move = 0;
            } else {
                min_move = i - h;
            }

            for j in min_move..max_move {
                if i == j {
                    continue;
                }
                let d = dist(&clusters[i], &clusters[j], block_len);
                if d < cutoff_dist {
                    best_moves_i.push((i, j, d));
                }
            }
            //Only allow very concordant moves
            if best_moves_i.len() == 1 {
                best_moves.extend(best_moves_i);
            }
        }
        if best_moves.is_empty() {
            break;
        } else {
            best_moves.sort_by(|x, y| x.2.partial_cmp(&y.2).unwrap());
            let best_move = best_moves[0];
            let high_ind = usize::max(best_move.0, best_move.1);
            let low_ind = usize::min(best_move.0, best_move.1);
            let removed_clust = clusters.remove(high_ind);
            for item in removed_clust {
                clusters[low_ind].push(item);
            }
        }
        println!("Clusters remaining {}", clusters.len());
    }

    let mut file = File::create("debug_clusters.txt").expect("Can't create file");
    write!(&mut file, "{:?}", &clusters).unwrap();

    let mut new_parts = vec![];
    let mut new_snp_ranges = vec![];
    for cluster in clusters {
        let mut new_snp_range = (SnpPosition::MAX, SnpPosition::MIN);
        let mut new_part = FxHashSet::default();
        for info in cluster {
            let index = info.3;
            let part = &parts[index];
            let range = &snp_endpoints[index];
            for frag in part.iter() {
                new_part.insert(*frag);
            }
            if range.0 < new_snp_range.0 {
                new_snp_range.0 = range.0;
            }
            if range.1 > new_snp_range.1 {
                new_snp_range.1 = range.1;
            }
        }
        new_parts.push(new_part);
        new_snp_ranges.push(new_snp_range);
    }

    for index in none_clusters {
        new_parts.push(parts[index].clone());
        new_snp_ranges.push(snp_endpoints[index].clone());
    }

    return (new_snp_ranges, new_parts);
}

fn find_overlapping_blocks<'a>(
    parts: &Vec<FxHashSet<&'a Frag>>,
    ol_cutoff: f64,
    snp_range_parts_vec: &Vec<(SnpPosition, SnpPosition)>,
) -> (
    FxHashMap<usize, Vec<Interval<SnpPosition, usize>>>,
    FxHashMap<usize, Vec<f64>>,
) {
    let mut interval_vec = vec![];
    type Iv = Interval<SnpPosition, usize>;
    for (i, _parts) in parts.iter().enumerate() {
        //        let mut range = (SnpPosition::MAX, SnpPosition::MIN);
        let range = &snp_range_parts_vec[i];
        //        for read in parts.iter() {
        //            if read.first_position < range.0 {
        //                range.0 = read.first_position;
        //            }
        //            if read.last_position >= range.1 {
        //                range.1 = read.last_position;
        //            }
        //        }
        interval_vec.push(Iv {
            start: range.0,
            stop: range.1,
            val: i,
        });
    }
    let interval_vec_clone = interval_vec.clone();

    let laps = Lapper::new(interval_vec);
    let mut all_overlaps = FxHashMap::default();
    let mut all_overlaps_percentage = FxHashMap::default();
    for (i, range) in interval_vec_clone.iter().enumerate() {
        let overlaps = laps.find(range.start, range.stop);
        let index = i;
        for interval_found in overlaps {
            let index_found = interval_found.val;
            let overlap_p = overlap_percent(
                range.start,
                range.stop,
                interval_found.start,
                interval_found.stop,
            );
            log::trace!(
                "overlap percent {} {} = {}; range {:?}, interval_found {:?}",
                i,
                interval_found.val,
                overlap_p,
                range,
                interval_found,
            );
            if overlap_p > ol_cutoff && index_found != i {
                let vec = all_overlaps.entry(index).or_insert(vec![]);
                vec.push(interval_found.clone());
                let p_vec = all_overlaps_percentage.entry(index).or_insert(vec![]);
                p_vec.push(overlap_p);
            }
        }
    }

    return (all_overlaps, all_overlaps_percentage);
}

pub fn get_hapq<'a>(
    parts: &Vec<FxHashSet<&'a Frag>>,
    snp_to_genome_pos: &'a Vec<usize>,
    snp_range_parts_vec: &Vec<(SnpPosition, SnpPosition)>,
    options: &Options,
) -> (Vec<u8>, Vec<f64>) {
    let mut hapqs = vec![];
    let mut purities = vec![];
    let mut weight = 0.;
    let mut error = 0.;
    let mut total_covs = vec![];
    let mut errs = vec![];
    for (i, part) in parts.iter().enumerate() {
        let (_cov, err, total_err, total_cov) = utils_frags::get_errors_cov_from_frags(
            part,
            snp_range_parts_vec[i].0,
            snp_range_parts_vec[i].1,
        );
        weight += total_cov;
        error += total_err;
        total_covs.push(total_cov);
        errs.push(err);
    }
    let avg_err = error / weight;
    let all_parts_block = utils_frags::hap_block_from_partition(&parts, true);
    let (all_ol, all_overlaps_p) = find_overlapping_blocks(parts, 0.05, snp_range_parts_vec);
    for i in 0..parts.len() {
        log::trace!(
            "hap {} has {} reads, covers range {}-{}",
            i,
            parts[i].len(),
            snp_range_parts_vec[i].0,
            snp_range_parts_vec[i].1
        );
        let mut max_ol = 0.;
        let mut max_penalty = 0.;
        let mut max_ind = 0;
        if all_overlaps_p.contains_key(&i) {
            for (_j, ol) in all_overlaps_p[&i].iter().enumerate() {
                let j = all_ol[&i][_j].val;
                let dist;
                let (same, diff) = utils_frags::distance_between_haplotypes(
                    &all_parts_block.blocks[i],
                    &all_parts_block.blocks[j],
                    &(SnpPosition::MIN, SnpPosition::MAX),
                );
                if (same + diff) == 0. {
                    dist = 1.;
                } else {
                    dist = diff as f64 / (same + diff) as f64;
                }
                if *ol * (1. - dist) > max_penalty {
                    max_ol = *ol;
                    max_penalty = *ol * (1. - dist);
                    max_ind = j;
                }
                log::trace!(
                    "{} overlapping {} with dist, same {} {}, max_ol {}",
                    i,
                    max_ind,
                    same,
                    diff,
                    max_ol
                );
            }
        }
        let mut range = (SnpPosition::MAX, 0);
        for read in parts[i].iter() {
            if read.first_position < range.0 {
                range.0 = read.first_position;
            }
            if read.last_position >= range.1 {
                range.1 = read.last_position;
            }
        }
        let base_range;
        if range.0 > range.1 {
            base_range = 0;
        } else {
            //#base_range = snp_to_genome_pos[(range.1 - 1) as usize] - snp_to_genome_pos[(range.0-1) as usize];
            base_range = snp_to_genome_pos[snp_range_parts_vec[i].1 as usize - 1]
                - snp_to_genome_pos[snp_range_parts_vec[i].0 as usize - 1];
            //base_range = snp_range_parts_vec[i].1 - snp_range_parts_vec[i].0;
        }

        //        let purity_score = utils_frags::stable_binom_cdf_p_rev(usize::max(total_covs[i] as usize,10000), (errs[i] * total_covs[i]) as usize, avg_err, 1.);
        //        let purity_val = f64::min(purity_score * 2.713f64.log(10.), 0.);
        let t1 = constants::HAPQ_CONSTANT * (1. - max_penalty);
        let t2 = f64::min(1., parts[i].len() as f64 / 3.);
        let t3 = f64::max(
            0.0,
            ((base_range as f64 / options.block_length as f64) + 1.).ln(),
        );
        log::trace!("hapq for hap {} = t1 t2 t3 {} {} {}", i, t1, t2, t3);
        let mut hapq = (t1 * t2 * t3) as usize;
        if parts[i].len() == 1 {
            hapq = 0;
        }
        hapqs.push(usize::min(hapq, 60) as u8);
        //        purities.push((-1. * purity_val) as u8);
        purities.push(errs[i] / avg_err)
    }
    return (hapqs, purities);
}
