use crate::types_structs;
use crate::types_structs::{Frag, HapBlock, SearchNode, SnpPosition};
use std::collections::binary_heap::BinaryHeap;
use std::mem;
use std::rc::Rc;
extern crate time;
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};

pub fn beam_search_phasing<'a>(
    clique: Vec<FxHashSet<&'a Frag>>,
    all_reads: &'a Vec<&Frag>,
    epsilon: f64,
    div_factor: f64,
    cutoff_value: f64,
    max_number_solns: usize,
    use_mec: bool,
    _use_ref_bias: bool,
) -> (
    FxHashMap<SnpPosition, FxHashSet<usize>>,
    Vec<FxHashSet<&'a Frag>>,
) {
    if all_reads.len() == 0 {
        return (FxHashMap::default(), vec![]);
    }
    let mut partition = clique.clone();
    let ploidy = clique.len();

    let first_block = utils_frags::hap_block_from_partition(&clique, true);
    let random_frag = &all_reads[0];

    let starting_freq = vec![1; clique.len()];
    let first_node = SearchNode {
        read: &random_frag,
        part: usize::MAX,
        score: 0.0,
        freqs: starting_freq,
        error_vec: vec![(0.0, 0.0); ploidy],
        block_id: 0,
        parent_node: None,
        current_pos: 0,
        broken_blocks: FxHashSet::default(),
    };

    //    let mut search_node_list = vec![(Rc::new(first_node), first_block)];
    let mut search_node_heap = BinaryHeap::new();
    search_node_heap.push((Rc::new(first_node), first_block));

    for i in 0..all_reads.len() {
        let mut max_num_soln_mut = max_number_solns;
        if i < 25 {
            max_num_soln_mut = ploidy * max_number_solns;
        }
        //        let mut search_node_list_next = vec![];
        let mut search_node_heap_next: BinaryHeap<(Rc<SearchNode>, HapBlock)> = BinaryHeap::new();
        let frag = &all_reads[i];
        let mut frag_in_clique = false;
        //If we use the clique construction, we don't
        //want to add multiple copies of the same fragment.
        for j in 0..ploidy {
            if clique[j].contains(frag) {
                frag_in_clique = true;
            }
        }
        if frag_in_clique {
            continue;
        }
        let current_startpos = frag.first_position;
        for (node, block) in search_node_heap.iter() {
            let mut p_value_list = vec![];
//            let mut same_diff_list = vec![];
            for part_index in 0..ploidy {
                let dist;
                let (same, diff) = utils_frags::distance_read_haplo_epsilon_empty(
                    frag,
                    &block.blocks[part_index],
                    epsilon,
                );
                dist = 1.0
                    * utils_frags::stable_binom_cdf_p_rev(
                        (same + diff) as usize,
                        diff as usize,
                        //2.0 * epsilon_ref * (1.0 - epsilon_ref),
                        epsilon,
                        div_factor,
                    );

//                same_diff_list.push((same,diff));
                p_value_list.push(dist);
            }
            let lse = utils_frags::log_sum_exp(&p_value_list);
            //let mut num_prune = 0;
            //let p_value_list_lse : Vec<f64> = p_value_list.iter().map(|x| (x - lse).exp()).collect();
            //dbg!(p_value_list_lse);
            for j in 0..ploidy {
                if p_value_list[j] - lse > cutoff_value {
                //if true{
                    //score is either the PEM or MEC score. I want to play around with using the
                    //iterative sum of p-values as well.
                    let (score, new_error_vec) =
                        read_to_node_value(node, frag, block, j, epsilon, div_factor, use_mec);
                    let new_node_score;
                    new_node_score = -score;

                    let mut new_node = types_structs::build_child_node(
                        frag,
                        j,
                        0,
                        Some(Rc::clone(node)),
                        new_error_vec,
                        new_node_score,
                        current_startpos,
                    );

                    let (broken_blocks_node, new_block) =
                        types_structs::build_truncated_hap_block(block, frag, j, current_startpos);
                    for index in broken_blocks_node {
                        new_node.broken_blocks.insert(index);
                    }
                    let mut project_exists = false;
                    for node in search_node_heap_next.iter() {
                        if node.1 == new_block && node.0.score >= new_node.score {
                            project_exists = true;
                        }
                    }
                    if !project_exists {
                        let toins = (Rc::new(new_node), new_block);
                        search_node_heap_next.push(toins);

                        if search_node_heap_next.len() > max_num_soln_mut {
                            search_node_heap_next.pop();
                        }
                    }
//                    log::trace!("read not pruned {}-{}, {:?}, {:?}, {}", i, j, &p_value_list, &same_diff_list, lse);
                }
                else{
//                    num_prune += 1;
//                    log::trace!("read pruned {}-{}, {:?}, {:?}, {}", i, j, &p_value_list, &same_diff_list, lse);
                }
            }
//            log::trace!("num prune {}, ploidy {}", num_prune, ploidy);
        }

        let _unused = mem::replace(&mut search_node_heap, search_node_heap_next);
    }

    let search_node_heap_to_list = search_node_heap.into_sorted_vec();
    let mut node_pointer = &search_node_heap_to_list[0].0;
    //    log::debug!("Partition count: {:?}", node_pointer.freqs);
    //    log::debug!("Best partition score: {}", node_pointer.score);

    let mut break_positions = FxHashMap::default();
    loop {
        let current_pos = node_pointer.current_pos;
        if node_pointer.broken_blocks.len() > 0 {
            let haps_to_break = break_positions
                .entry(current_pos)
                .or_insert(FxHashSet::default());
            for index in node_pointer.broken_blocks.iter() {
                haps_to_break.insert(*index);
            }
        }
        if node_pointer.parent_node.is_none() {
            break;
        } else {
            partition[node_pointer.part].insert(node_pointer.read);
            //            log::trace!(
            //                "Read {} added to partition {}",
            //                node_pointer.read.id,
            //                node_pointer.part
            //            );
            node_pointer = node_pointer.parent_node.as_ref().unwrap();
        }
    }
    //dbg!(break_positions);
    return (break_positions, partition);
}

fn read_to_node_value(
    node: &SearchNode,
    frag: &Frag,
    block: &HapBlock,
    part_index: usize,
    epsilon: f64,
    _div_factor: f64,
    _use_mec: bool,
) -> (f64, Vec<(f64, f64)>) {
    let ploidy = block.blocks.len();
    //    let (same, diff) = utils_frags::distance_read_haplo(frag, &block.blocks[part_index]);
    let (same, diff) =
        utils_frags::distance_read_haplo_epsilon_empty(frag, &block.blocks[part_index], epsilon);
    let mut new_error_vec = vec![];
    for i in 0..ploidy {
        if i == part_index {
            new_error_vec.push((node.error_vec[i].0 + same, node.error_vec[i].1 + diff));
        } else {
            new_error_vec.push(node.error_vec[i]);
        }
    }
    let mec: f64 = new_error_vec.iter().map(|x| x.1).sum();
    return (
        -1.0 * mec,
        //            local_clustering::get_mec_score(&new_error_vec, &vec![0; 1], epsilon, div_factor),
        new_error_vec,
    );
}

