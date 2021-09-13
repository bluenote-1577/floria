use crate::local_clustering;
use crate::types_structs;
use crate::types_structs::{Frag, HapBlock, SearchNode};
use crate::vcf_polishing;
use rand::prelude::*;
use rand_pcg::Pcg64;
use statrs::distribution::{ChiSquared, Univariate};
use std::collections::BinaryHeap;
use std::mem;
use std::ptr;
use std::rc::Rc;
extern crate time;
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};

pub fn beam_search_phasing<'a>(
    clique: Vec<FxHashSet<&'a Frag>>,
    all_reads: &'a Vec<Frag>,
    epsilon: f64,
    div_factor: f64,
    cutoff_value: f64,
    max_number_solns: usize,
    use_mec: bool,
    use_supp_anchor: bool,
    use_ref_bias: bool,
) -> Vec<FxHashSet<&'a Frag>> {
    let mut partition = clique.clone();
    let ploidy = clique.len();

    let first_block = utils_frags::hap_block_from_partition(&clique);
    let random_frag = &all_reads[0];

    let mut supp_anchors = FxHashMap::default();
    for (i, part) in clique.iter().enumerate() {
        for read in part.iter() {
            if let Some(cont) = &read.supp_aln {
                let map_ent = supp_anchors.entry(cont).or_insert(vec![]);
                map_ent.push(i);
            }
        }
        //Allow frags with supp alignments to go to empty partition too
        if part.len() == 0{
            for (contig,vec) in supp_anchors.iter_mut(){
                vec.push(i);
            }
        }
    }

    let starting_freq = vec![1; clique.len()];
    let first_node = SearchNode {
        read: &random_frag,
        part: usize::MAX,
        score: 0.0,
        freqs: starting_freq,
        error_vec: vec![(0, 0); ploidy],
        block_id: 0,
        parent_node: None,
    };

    let mut search_node_list = vec![(Rc::new(first_node), first_block)];

    for i in 0..all_reads.len() {
        let mut min_score = f64::MIN;
        let mut search_node_list_next = vec![];
        let frag = &all_reads[i];
        let mut frag_in_clique = false;
        for j in 0..ploidy {
            if clique[j].contains(&frag) {
                frag_in_clique = true;
            }
        }
        if frag_in_clique {
            continue;
        }
        let current_startpos = frag.first_position;
        for (node, block) in search_node_list.iter() {
            let mut p_value_list = vec![];
            for part_index in 0..ploidy {
                let mut dist;
                if use_ref_bias {
                    let ((same_ref, diff_ref), (same_alt, diff_alt)) =
                        utils_frags::distance_read_haplo_ref_wild(frag, &block.blocks[part_index]);
                    let epsilon_ref = 0.01;
                    let epsilon_alt = 2.0 * epsilon;
                    let dist_ref = 1.0
                        * utils_frags::stable_binom_cdf_p_rev(
                            (same_ref + diff_ref) as usize,
                            diff_ref as usize,
                            //2.0 * epsilon_ref * (1.0 - epsilon_ref),
                            epsilon_ref,
                            div_factor,
                        );
                    let dist_alt = 1.0
                        * utils_frags::stable_binom_cdf_p_rev(
                            (same_alt + diff_alt) as usize,
                            diff_alt as usize,
                            //2.0 * epsilon_alt* (1.0 - epsilon_alt),
                            epsilon_alt,
                            div_factor,
                        );
                    dist = dist_alt + dist_ref;
                }
                else{
                    let (same,diff) = utils_frags::distance_read_haplo(frag, &block.blocks[part_index]);
                    dist = 1.0
                        * utils_frags::stable_binom_cdf_p_rev(
                            (same + diff) as usize,
                            diff as usize,
                            //2.0 * epsilon_ref * (1.0 - epsilon_ref),
                            epsilon,
                            div_factor,
                        );
                }
            

                if use_supp_anchor {
                    if let Some(cont) = &frag.supp_aln {
                        if supp_anchors.contains_key(cont) {
                            if !supp_anchors.get(cont).unwrap().contains(&part_index) {
                                dist = f64::MIN;
                            }
                        }
                    }
                }
                p_value_list.push(dist);
            }
            let lse = utils_frags::log_sum_exp(&p_value_list);
            //let p_value_list_lse : Vec<f64> = p_value_list.iter().map(|x| (x - lse).exp()).collect();
            //dbg!(p_value_list_lse);
            for j in 0..ploidy {
                if p_value_list[j] - lse > cutoff_value {
                    let (pem_score, new_error_vec) =
                        read_to_node_value(node, frag, block, j, epsilon, div_factor, use_mec);
                    let new_node_score;

                    if !use_mec {
                        //new_node_score = node.score + p_value_list[j];
                        new_node_score = pem_score;
                    } else {
                        new_node_score = pem_score;
                    }

                    let new_node = types_structs::build_child_node(
                        frag,
                        j,
                        0,
                        Some(Rc::clone(node)),
                        new_error_vec,
                        new_node_score,
                    );

                    if search_node_list_next.len() >= max_number_solns {
                        if pem_score < min_score {
                            continue;
                        } else {
                            let new_block = types_structs::build_truncated_hap_block(
                                block,
                                frag,
                                j,
                                current_startpos,
                            );
                            let toins = (Rc::new(new_node), new_block);

                            match search_node_list_next.binary_search_by(
                                |x: &(Rc<SearchNode>, HapBlock)| {
                                    x.0.score
                                        .partial_cmp(&(*(toins.0)).score)
                                        .expect("Couldn't compare")
                                },
                            ) {
                                Ok(pos) => search_node_list_next.insert(pos, toins),
                                Err(pos) => search_node_list_next.insert(pos, toins),
                            }
                            search_node_list_next.pop();
                        }
                    } else {
                        let new_block = types_structs::build_truncated_hap_block(
                            block,
                            frag,
                            j,
                            current_startpos,
                        );

                        search_node_list_next.push((Rc::new(new_node), new_block));
                    }

                    //first time get max number of solutions, sort
                    if search_node_list_next.len() == max_number_solns {
                        search_node_list_next
                            .sort_by(|a, b| b.0.score.partial_cmp(&a.0.score).unwrap());
                    }

                    min_score = search_node_list_next.iter().last().unwrap().0.score;
                }
            }
        }

        search_node_list_next.sort_by(|a, b| b.0.score.partial_cmp(&a.0.score).unwrap());
        mem::replace(&mut search_node_list, search_node_list_next);
        if search_node_list.len() > max_number_solns {
            search_node_list.drain(max_number_solns..);
        }

        if i % 100 == 0 {
            log::trace!(
                "Number of solutions before trimming: {}, {}, {}",
                search_node_list.len(),
                i,
                all_reads.len()
            );

            let test_pointer = &search_node_list[0].0;
            log::trace!("Partition best count:{}, {:?}", i, test_pointer.freqs);
            let test_pointer = &search_node_list.iter().last().unwrap().0;
            log::trace!("Partition worst count:{}, {:?}", i, test_pointer.freqs);
        }
    }

    let mut node_pointer = &search_node_list[0].0;
    log::debug!("Partition count: {:?}", node_pointer.freqs);
    log::debug!("Best partition score: {}", node_pointer.score);

    loop {
        if node_pointer.parent_node.is_none() {
            break;
        } else {
            partition[node_pointer.part].insert(node_pointer.read);
            //            log::trace!(
            //                "Read {} added to partition {}",
            //                node_pointer.read.id,
            //                node_pointer.part
            //            );
            node_pointer = &node_pointer.parent_node.as_ref().unwrap();
        }
    }
    partition
}

fn read_to_node_value(
    node: &SearchNode,
    frag: &Frag,
    block: &HapBlock,
    part_index: usize,
    epsilon: f64,
    div_factor: f64,
    use_mec: bool,
) -> (f64, Vec<(usize, usize)>) {
    let ploidy = block.blocks.len();
    let (same, diff) = utils_frags::distance_read_haplo(frag, &block.blocks[part_index]);
    let mut new_error_vec = vec![];
    for i in 0..ploidy {
        if i == part_index {
            new_error_vec.push((node.error_vec[i].0 + same, node.error_vec[i].1 + diff));
        } else {
            new_error_vec.push(node.error_vec[i]);
        }
    }
    if use_mec {
        return (
            local_clustering::get_mec_score(&new_error_vec, &vec![0; 1], epsilon, div_factor),
            new_error_vec,
        );
    } else {
        return (
            local_clustering::get_pem_score(&new_error_vec, &vec![0; 1], epsilon, div_factor),
            new_error_vec,
        );
    }
}

pub fn get_read_graph<'a>(
    vec_all_reads: &Vec<&'a Frag>,
    epsilon: f64,
    use_binomial_dist: bool,
) -> (Vec<(f64, i32, i32)>, Vec<Vec<(f64, i32)>>) {
    let mut vec_all_edges = Vec::new();
    let mut adj_list_edges = Vec::new();
    for _i in 0..vec_all_reads.len() {
        adj_list_edges.push(Vec::new());
    }

    //Get local read-read graph, a.k.a the distance matrix. In the future, we can speed this up by precomputing a
    //
    //global read-read graph or precomputing the distances between reads.
    //
    for (i, r1) in vec_all_reads.iter().enumerate() {
        for j in i + 1..vec_all_reads.len() {
            let r2 = &vec_all_reads[j];
            if !utils_frags::check_overlap(r1, r2) {
                continue;
            }

            let dist;
            let (same, mec_dist) = utils_frags::distance(r1, r2);

            //BINOMIAL DIST
            if use_binomial_dist {
                dist = -1.0
                    * utils_frags::stable_binom_cdf_p_rev(
                        (same + mec_dist) as usize,
                        mec_dist as usize,
                        2.0 * epsilon * (1.0 - epsilon),
                        100.0,
                    );
            } else {
                dist = mec_dist as f64;
            }

            let i_type = i as i32;
            let j_type = j as i32;
            vec_all_edges.push((dist, i_type, j_type));
            let edge_list1 = &mut adj_list_edges[i];
            edge_list1.push((dist, j_type));
            let edge_list2 = &mut adj_list_edges[j];
            edge_list2.push((dist, i_type));
        }
    }

    return (vec_all_edges, adj_list_edges);
}

pub fn get_initial_clique<'a>(
    all_reads: &'a Vec<Frag>,
    ploidy: usize,
    epsilon: f64,
) -> Vec<FxHashSet<&'a Frag>> {
    let starting_clique_reads = local_clustering::find_reads_in_interval(1, 20, all_reads);
    let use_binomial_dist = true;
    let mut partition = Vec::new();
    //let vec_all_reads: Vec<_> = all_reads.iter().collect();
    let vec_all_reads: Vec<_> = starting_clique_reads.into_iter().collect();
    let (mut vec_all_edges, _adj_list_edges) =
        get_read_graph(&vec_all_reads, epsilon, use_binomial_dist);

    if vec_all_edges.len() == 0 {
        let mut clusters: Vec<FxHashSet<&Frag>> = Vec::new();
        for _i in 0..ploidy {
            clusters.push(FxHashSet::default());
        }

        return clusters;
    }

    //Finding max clique
    //println!("{:?}",vec_all_edges);
    vec_all_edges.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let best_edge = vec_all_edges.last().unwrap();

    let mut used_vertices = FxHashSet::default();
    used_vertices.insert(best_edge.1);
    used_vertices.insert(best_edge.2);

    //Greedily find a max clique once the first two vertices are found. Greedily do this by adding
    //edges which maximizes the minimum of the distance between all other vertices in the putative clique.
    for _i in 0..ploidy - 2 {
        //min_dist_map contains the minimum distance from each vertex in the putative clique to
        //the other every other vertex. I.e. the key is a different vertex and the value is the
        //minimum distance from the clique to the key.

        let mut vertices_meeting_clique_map = FxHashMap::default();
        let mut min_dist_map = FxHashMap::default();

        for edge in vec_all_edges.iter() {
            if used_vertices.contains(&edge.1) && !used_vertices.contains(&edge.2) {
                let met_cliques_set = vertices_meeting_clique_map
                    .entry(edge.2)
                    .or_insert(FxHashSet::default());
                met_cliques_set.insert(edge.1);

                if min_dist_map.contains_key(&edge.2) {
                    if *min_dist_map.get(&edge.2).unwrap() < edge.0 {
                        continue;
                    }
                }
                min_dist_map.insert(edge.2, edge.0);
            } else if used_vertices.contains(&edge.2) && !used_vertices.contains(&edge.1) {
                let met_cliques_set = vertices_meeting_clique_map
                    .entry(edge.1)
                    .or_insert(FxHashSet::default());
                met_cliques_set.insert(edge.2);

                if min_dist_map.contains_key(&edge.1) {
                    if *min_dist_map.get(&edge.1).unwrap() < edge.0 {
                        continue;
                    }
                }
                min_dist_map.insert(edge.1, edge.0);
            }
        }

        let mut sorted_dict_to_vec: Vec<_> = min_dist_map.into_iter().collect();
        sorted_dict_to_vec.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        if sorted_dict_to_vec.len() == 0 {
            continue;
        }

        for vertex in sorted_dict_to_vec.iter().rev() {
            if *vertices_meeting_clique_map.get(&vertex.0).unwrap() == used_vertices {
                let best_vertex = vertex;
                used_vertices.insert(best_vertex.0);
                break;
            }
        }
    }

    let mut clusters = Vec::new();
    for vertex in used_vertices.iter() {
        let mut cluster = FxHashSet::default();
        cluster.insert(*vertex);
        clusters.push(cluster);
    }
    if clusters.len() < ploidy {
        loop {
            clusters.push(FxHashSet::default());
            if clusters.len() == ploidy {
                break;
            }
        }
    }

    for cluster in clusters.iter() {
        let mut frag_set = FxHashSet::default();
        for vertex in cluster.iter() {
            let vertex_usize = *vertex as usize;
            frag_set.insert(vec_all_reads[vertex_usize]);
        }
        partition.push(frag_set)
    }

    for frag_set in partition.iter() {
        for frag in frag_set.iter() {
            log::debug!(
                "Clique: {},{},{}",
                frag.id,
                frag.first_position,
                frag.last_position
            );
        }
    }

    partition
}

pub fn get_initial_from_anchor<'a>(
    all_reads: &'a Vec<Frag>,
    ploidy: usize,
    epsilon: f64,
    contig_anchors: &Vec<String>,
) -> Vec<FxHashSet<&'a Frag>> {
    //get anchors
    let num_anchors = contig_anchors.len();
    let mut used_anchors = vec![];
    if num_anchors > ploidy {
        panic!("Number of anchor contigs greater than ploidy. Exiting");
    }
    let mut starting_clique_reads = FxHashMap::default();
    let mut starting_clique_reads_vec = vec![];
    for frag in all_reads.iter() {
        if let Some(contig) = &frag.supp_aln {
            if contig_anchors.contains(contig) {
                let anchor_to_frags = starting_clique_reads.entry(contig).or_insert(vec![]);
                anchor_to_frags.push(frag);
                starting_clique_reads_vec.push(frag);
            }
        }
    }

    let use_binomial_dist = true;
    let vec_all_reads = starting_clique_reads_vec;
    let mut partition = Vec::new();
    let (mut vec_all_edges, _adj_list_edges) =
        get_read_graph(&vec_all_reads, epsilon, use_binomial_dist);

    //Return empty partitions if no edges are found
    if vec_all_edges.len() == 0 {
        let mut clusters: Vec<FxHashSet<&Frag>> = Vec::new();
        for _i in 0..ploidy {
            clusters.push(FxHashSet::default());
        }

        log::debug!("Empty initial partitioning!");
        return clusters;
    }

    //Finding max clique
    //println!("{:?}",vec_all_edges);
    vec_all_edges.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let mut best_edge = vec_all_edges.first().unwrap();
    for edge in vec_all_edges.iter() {
        let cont1 = vec_all_reads[edge.1 as usize].supp_aln.as_ref().unwrap();
        let cont2 = vec_all_reads[edge.2 as usize].supp_aln.as_ref().unwrap();
        if cont1 != cont2 {
            best_edge = edge;
            used_anchors.push(cont1);
            used_anchors.push(cont2);
            log::debug!("Clique between anchors {}-{}", cont1, cont2);
            break;
        }
    }

    let mut used_vertices = FxHashSet::default();
    used_vertices.insert(best_edge.1);
    used_vertices.insert(best_edge.2);

    //Greedily find a max clique once the first two vertices are found. Greedily do this by adding
    //edges which maximizes the minimum of the distance between all other vertices in the putative clique.
    for _i in 0..ploidy - 2 {
        //min_dist_map contains the minimum distance from each vertex in the putative clique to
        //the other every other vertex. I.e. the key is a different vertex and the value is the
        //minimum distance from the clique to the key.

        let mut vertices_meeting_clique_map = FxHashMap::default();
        let mut min_dist_map = FxHashMap::default();

        for edge in vec_all_edges.iter() {
            if used_vertices.contains(&edge.1) && !used_vertices.contains(&edge.2) {
                let met_cliques_set = vertices_meeting_clique_map
                    .entry(edge.2)
                    .or_insert(FxHashSet::default());
                met_cliques_set.insert(edge.1);

                if min_dist_map.contains_key(&edge.2) {
                    if *min_dist_map.get(&edge.2).unwrap() < edge.0 {
                        continue;
                    }
                }
                min_dist_map.insert(edge.2, edge.0);
            } else if used_vertices.contains(&edge.2) && !used_vertices.contains(&edge.1) {
                let met_cliques_set = vertices_meeting_clique_map
                    .entry(edge.1)
                    .or_insert(FxHashSet::default());
                met_cliques_set.insert(edge.2);

                if min_dist_map.contains_key(&edge.1) {
                    if *min_dist_map.get(&edge.1).unwrap() < edge.0 {
                        continue;
                    }
                }
                min_dist_map.insert(edge.1, edge.0);
            }
        }

        let mut sorted_dict_to_vec: Vec<_> = min_dist_map.into_iter().collect();
        sorted_dict_to_vec.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        if sorted_dict_to_vec.len() == 0 {
            continue;
        }

        for vertex in sorted_dict_to_vec.iter().rev() {
            if *vertices_meeting_clique_map.get(&vertex.0).unwrap() == used_vertices {
                let best_vertex = vertex;
                if used_anchors.len() != num_anchors {
                    let frag = &vec_all_reads[best_vertex.1 as usize];
                    let cont = frag.supp_aln.as_ref().unwrap();
                    if !used_anchors.contains(&cont) {
                        used_vertices.insert(best_vertex.0);
                        used_anchors.push(&cont);
                        log::debug!("Clique for anchors {}", cont);
                        break;
                    }
                } else {
                    //When using GLOPP, don't need a full clique output. Can just return partial
                    //clique. 
                    //used_vertices.insert(best_vertex.0);
                    break;
                }
            }
        }
    }

    let mut clusters = Vec::new();
    for vertex in used_vertices.iter() {
        let mut cluster = FxHashSet::default();
        cluster.insert(*vertex);
        clusters.push(cluster);
    }

    if clusters.len() < ploidy {
        loop {
            clusters.push(FxHashSet::default());
            if clusters.len() == ploidy {
                break;
            }
        }
    }

    for cluster in clusters.iter() {
        let mut frag_set = FxHashSet::default();
        for vertex in cluster.iter() {
            let vertex_usize = *vertex as usize;
            frag_set.insert(vec_all_reads[vertex_usize]);
        }
        partition.push(frag_set)
    }

    for frag_set in partition.iter() {
        for frag in frag_set.iter() {
            log::debug!(
                "Clique: {},{},{}",
                frag.id,
                frag.first_position,
                frag.last_position
            );
        }
    }

    partition
}