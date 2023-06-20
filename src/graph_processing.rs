use crate::constants;
use log::Level::{Debug, Trace};
use log::{log_enabled};
use crate::file_writer;
use crate::global_clustering;
use crate::local_clustering;
use crate::part_block_manip;
use crate::types_structs::{Frag, GnPosition, HapNode, SnpPosition, TraceBackNode, VcfProfile, Options, FlowUpVec};
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use std::sync::Mutex;
//use osqp::{CscMatrix, Problem, Settings};
use petgraph::algo;
use petgraph::dot::Dot;
use petgraph::prelude::*;
use std::fs::File;
use std::io::Write;
use std::mem;


fn update_hap_graph(hap_graph: &mut Vec<Vec<HapNode>>) {
    //    let pseudo_count = 10.;
    let mut out_edges_block_hap = vec![];
    for i in 0..hap_graph.len() - 1 {
        let mut out_edges_block = vec![];
        let hap_block1 = &hap_graph[i];
        let hap_block2 = &hap_graph[i + 1];
        for hap_node1 in hap_block1 {
            let mut out_weights = vec![0.0; hap_block2.len()];
            //We have to make sure that ambiguous reads do not
            //get accounted for. For short reads, this is important.
            for read in hap_node1.frag_set.iter() {
                let mut read_to_hap_sim = vec![];
                let mut hap_id_in = usize::MAX;
                for (l, hap_node2) in hap_block2.iter().enumerate() {
                    if hap_node2.frag_set.contains(read) {
                        hap_id_in = l;
                    }
                    let (_same, diff) = utils_frags::distance_read_haplo(read, &hap_node2.hap_map);
                    read_to_hap_sim.push((diff, l));
                }
                read_to_hap_sim.sort();
                if read_to_hap_sim.len() > 1 {
                    if read_to_hap_sim[0].0 != read_to_hap_sim[1].0 {
                        if hap_id_in != usize::MAX {
                            out_weights[hap_id_in] += 1.;
                        }
                    } else {
                        //                        dbg!(i,&read.id, &read.first_position, &read.last_position);
                    }
                } else {
                    if hap_id_in != usize::MAX {
                        out_weights[hap_id_in] += 1.;
                    }
                }
            }
            //            let sum: f64 = out_weights.iter().sum();
            //            let _normalized_out_weights: Vec<f64> = out_weights.iter().map(|x| x / sum).collect();
            let mut out_edges_hap = vec![];
            for l in 0..hap_block2.len() {
                if out_weights[l] >= constants::MIN_SHARED_READS_UNAMBIG {
                    //                    out_edges_hap.push((l, normalized_out_weights[l]));
                    out_edges_hap.push((l, out_weights[l]));
                }
            }
            //            let new_normalization: f64 = out_edges_hap.iter().map(|x| x.1).sum();
            let new_normalization: f64 = 1.;
            let new_normal_out_edges_hap: Vec<_> = out_edges_hap
                .iter()
                .map(|x| (x.0, x.1 / new_normalization))
                .collect();
            out_edges_block.push(new_normal_out_edges_hap);
        }
        out_edges_block_hap.push(out_edges_block);
    }

    let num_blocks = hap_graph.len();
    for (i, hap_block) in hap_graph.iter_mut().enumerate() {
        if i != num_blocks - 1 {
            let out_edges_block = &out_edges_block_hap[i];
            for (j, hap_node) in hap_block.iter_mut().enumerate() {
                let out_edges_hap = &out_edges_block[j];
                for (k, prob) in out_edges_hap {
                    hap_node.out_edges.push((*k, *prob));
                    log::trace!("BLOCK {}: {}-{} weight {}", i, j, k, prob);
                }
            }
        }
        if i != 0 {
            let prev_out_edges_block = &out_edges_block_hap[i - 1];
            for (j, out_edges_hap) in prev_out_edges_block.iter().enumerate() {
                for (k, prob) in out_edges_hap.iter() {
                    let hap_node = &mut hap_block[*k];
                    hap_node.in_edges.push((j, *prob));
                }
            }
        }
    }
}


fn get_local_hap_blocks<'a>(
    all_frags: &'a Vec<Frag>,
    snp_to_genome_pos: &'a Vec<GnPosition>,
    floria_out_dir: &str,
    j: usize,
    snp_range_vec: &Vec<(SnpPosition, SnpPosition)>,
    options: &Options
) -> Option<Vec<Vec<HapNode<'a>>>> {
    let max_ploidy = options.max_ploidy;
    let epsilon = options.epsilon;
    let max_number_solns = options.max_number_solns;
    let ploidy_start = 1;
    let ploidy_end = max_ploidy + 1;
    let num_ploidies = ploidy_end - ploidy_start;
    let mut mec_vector = vec![0.; num_ploidies];
    let mut parts_vector = vec![];
    let mut expected_errors_ref = vec![];
    let mut endpoints_vector = vec![];
    let reads = local_clustering::find_reads_in_interval(
        snp_range_vec[j].0,
        snp_range_vec[j].1,
        all_frags,
        usize::MAX,
    );

    let mut best_ploidy = ploidy_start;
    if reads.is_empty() {
        return None;
    }
    for ploidy in ploidy_start..ploidy_end {
        best_ploidy = ploidy;
        let mut num_alleles = 0.0;
        let mut vec_reads_own = vec![];
        for read in reads.iter() {
            vec_reads_own.push(*read);
        }
        vec_reads_own.sort();
        let (break_pos, part) = global_clustering::beam_search_phasing(
            vec![FxHashSet::default(); ploidy],
            &vec_reads_own,
            epsilon,
            constants::DIV_FACTOR,
            //            f64::MIN,
            constants::PROB_CUTOFF.ln(),
            max_number_solns,
            true,
            false,
        );

        //            let optimized_part = part;
        let (_new_score, optimized_part, _block) =
            local_clustering::optimize_clustering(part, epsilon, constants::NUM_ITER_OPTIMIZE);

        let binom_vec = local_clustering::get_mec_stats_epsilon_no_phred(&optimized_part, epsilon);
//        let binom_vec = local_clustering::get_mec_stats_epsilon_yes_phred(&optimized_part, epsilon);
        for (good, bad) in binom_vec {
            mec_vector[ploidy - ploidy_start] += bad;
            num_alleles += good;
            num_alleles += bad;
        }

        let split_part_merge;
        let split_part_endpoints;
        if constants::WEIRD_SPLIT{
            let split_part =
                utils_frags::split_part_using_breaks(&break_pos, &optimized_part, &all_frags);
            let endpoints;
            if j != 0 {
                endpoints = (snp_range_vec[j].0, snp_range_vec[j].1);
            } else {
                endpoints = (snp_range_vec[j].0, snp_range_vec[j].1);
            }
            let merge_result = merge_split_parts(split_part, break_pos, endpoints);
            split_part_merge = merge_result.0;
            split_part_endpoints = merge_result.1;
        }
        else{
            split_part_merge = vec![optimized_part];
            split_part_endpoints = vec![snp_range_vec[j]];
        }

        let mut ind_parts = vec![];
        for part in split_part_merge {
            let mut ind_part = vec![];
            for set in part {
                let index_set: FxHashSet<usize> = set.iter().map(|x| x.counter_id).collect();
                ind_part.push(index_set);
            }
            ind_parts.push(ind_part);
        }
        parts_vector.push(ind_parts);
        endpoints_vector.push(split_part_endpoints);

        expected_errors_ref.push(num_alleles as f64 * epsilon);

        if ploidy > ploidy_start {
            //            let mec_threshold = 1.0 / (1.0 - error_rate) / (1.0 + 1.0 / (ploidy + 1) as f64);
            //            log::trace!(
            //                "MEC vector {:?}, error_thresh {:?}",
            //                &mec_vector,
            //                expected_errors_ref
            //            );
            let mec_threshold;
            if options.ploidy_sensitivity == 1{
                mec_threshold =
                1.0 / (1.0 - epsilon) / (1.0 + 1.0 / ((ploidy as f64).powf(0.50) + 1.00) as f64);

            }
            else if options.ploidy_sensitivity == 2{
//                mec_threshold =
//                1.0 / (1.0 - epsilon) / (1.0 + 1.0 / ((ploidy as f64).powf(0.75) + 1.32) as f64);
                mec_threshold =
                1.0 / (1.0 - epsilon) / (1.0 + 1.0 / ((ploidy as f64).powf(1.00) + 1./3.) as f64);

            }
            else{
                mec_threshold =
                1.0 / (1.0 - epsilon) / (1.0 + 1.0 / ((ploidy as f64).powf(1.00) + 1.00) as f64);

            }
            log::trace!(
                "Expected MEC ratio {}, observed MEC ratio {}",
                mec_threshold,
                mec_vector[ploidy - ploidy_start] as f64
                    / mec_vector[ploidy - ploidy_start - 1] as f64
            );
            if (mec_vector[ploidy - ploidy_start] as f64
                / mec_vector[ploidy - ploidy_start - 1] as f64)
                < mec_threshold
            {
                //do nothing
            } 
            else{
                if options.stopping_heuristic{
                    log::trace!("MEC decrease thereshold, returning ploidy {}.", ploidy - 1);
                    best_ploidy -= 1;
                    break;
                }
            }
            if mec_vector[ploidy - ploidy_start] < expected_errors_ref[ploidy - ploidy_start] {
                log::trace!("MEC error threshold, returning ploidy {}.", ploidy);
                break;
            }
        } else {
            if mec_vector[ploidy - ploidy_start] < expected_errors_ref[ploidy - ploidy_start] {
                log::trace!("MEC error threshold, returning ploidy {}.", ploidy);
                break;
            }
        }
    }

    if best_ploidy == max_ploidy{
        log::debug!("Max ploidy {} reached at SNPs {:?} . Consider increasing the maximum ploidy (-p option)",max_ploidy, &snp_range_vec[j]);
    }

    log::trace!("DIFF\t{}\t{}", mec_vector[0], mec_vector[1]);

    log::trace!(
        "MEC vector {:?}, error_thresh {:?}, SNPs interval  {} {}",
        &mec_vector,
        expected_errors_ref,
        snp_range_vec[j].0,
        snp_range_vec[j].1,
    );

    let best_parts = mem::take(&mut parts_vector[best_ploidy - ploidy_start]);
    let best_endpoints = mem::take(&mut endpoints_vector[best_ploidy - ploidy_start]);
    let local_part_dir = format!("{}/local_parts/", floria_out_dir);
    let mut hap_node_blocks = vec![];
    //    if best_ploidy == 1{
    //        return None;
    //    }

    for (l, best_part) in best_parts.iter().enumerate() {
        let mut hap_node_block = vec![];
        let mut frag_best_part = vec![];

        for ind_part in best_part.iter() {
            let frag_set: FxHashSet<&Frag> = ind_part.iter().map(|x| &all_frags[*x]).collect();
            frag_best_part.push(frag_set.clone());
            let mut hap_node = HapNode::new(frag_set, best_endpoints[l]);
            hap_node.row = hap_node_block.len();
            hap_node_block.push(hap_node);
        }
        hap_node_blocks.push(hap_node_block);

        if log_enabled!(Debug) || log_enabled!(Trace){
            file_writer::write_all_parts_file(
                &frag_best_part,
                "",
                &vec![],
                &local_part_dir,
                &format!("{}-{}-{}-{}", j, l, snp_range_vec[j].0, best_ploidy),
                &snp_to_genome_pos,
                &vec![],
                &vec![],
            );
        }
    }

    return Some(hap_node_blocks);
}

fn process_chunks(mut chunks: Vec<(usize, Vec<Vec<HapNode>>)>) -> Vec<Vec<HapNode>> {
    chunks.sort_by(|x, y| x.0.cmp(&y.0));
    let mut return_blocks = vec![];
    for (_i, chunk) in chunks {
        for block in chunk {
            return_blocks.push(block);
        }
    }
    let mut id_counter = 0;
    for (i, block) in return_blocks.iter_mut().enumerate() {
        for node in block.iter_mut() {
            node.column = i;
            node.id = id_counter;
            id_counter += 1;
        }
    }
    return return_blocks;
}

pub fn generate_hap_graph<'a>(
    all_frags: &'a Vec<Frag>,
    snp_to_genome_pos: &'a Vec<usize>,
    floria_out_dir: String,
    options: &Options,
) -> Vec<Vec<HapNode<'a>>> {
    let block_length = options.block_length;
    let minimal_density = options.snp_density;

    let iter_vec: Vec<(SnpPosition, SnpPosition)> = utils_frags::get_range_with_lengths(
        snp_to_genome_pos,
        block_length,
        block_length / 3,
        minimal_density,
    );

    let interval_vec = iter_vec.to_vec();
    log::trace!("SNP Endpoints {:?}", &interval_vec);

    let block_chunks: Mutex<Vec<_>> = Mutex::new(vec![]);
    (0..interval_vec.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|j| {
            let block_chunk = get_local_hap_blocks(
                all_frags,
                snp_to_genome_pos,
                &floria_out_dir,
                j,
                &interval_vec,
                options,
            );
            //If the ploidy is 1, we return nothing.
            if !block_chunk.is_none() {
                let mut locked = block_chunks.lock().unwrap();
                locked.push((j, block_chunk.unwrap()));
            }
        });

    let block_chunks = block_chunks.into_inner().unwrap();
    let mut hap_node_blocks = process_chunks(block_chunks);
    if hap_node_blocks.is_empty() {
        return vec![];
    } else {
        update_hap_graph(&mut hap_node_blocks);
    }
    hap_node_blocks
}

fn merge_split_parts(
    mut split_part: Vec<Vec<FxHashSet<&Frag>>>,
    break_pos: FxHashMap<SnpPosition, FxHashSet<usize>>,
    original_snp_endpoints: (SnpPosition, SnpPosition),
) -> (Vec<Vec<FxHashSet<&Frag>>>, Vec<(SnpPosition, SnpPosition)>) {
    let mut breaks_with_min_sorted = vec![];
    for (key, value) in break_pos.iter() {
        if value.len() > 1 {
            breaks_with_min_sorted.push(key);
        }
    }
    breaks_with_min_sorted.sort();
    assert!(breaks_with_min_sorted.len() == split_part.len() - 1);

    let mut snp_breakpoints = vec![];
    let mut split_part_merge = vec![];
    let mut tomerge = vec![];
    let mut left_endpoint_merged = original_snp_endpoints.0;
    for k in 0..split_part.len() - 1 {
        let total_cov_1: usize = split_part[k].iter().map(|x| x.len()).sum();
        let total_cov_2: usize = split_part[k + 1].iter().map(|x| x.len()).sum();
        let cov_rat = total_cov_1 as f64 / (total_cov_2 + total_cov_1) as f64;
        //TODO Merging stuff
        if *breaks_with_min_sorted[k] <= left_endpoint_merged {
            tomerge.push(k);
            //            snp_breakpoints.push((left_endpoint_merged, *breaks_with_min_sorted[k]));
            //            left_endpoint_merged = *breaks_with_min_sorted[k];
            continue;
        }
        if cov_rat > 0.95 || cov_rat < 0.05 {
            tomerge.push(k);
        } else {
            snp_breakpoints.push((left_endpoint_merged, *breaks_with_min_sorted[k]));
            left_endpoint_merged = *breaks_with_min_sorted[k];
        }
    }
    snp_breakpoints.push((left_endpoint_merged, original_snp_endpoints.1));

    let mut new_part_empty = true;
    let mut new_part = vec![FxHashSet::default(); split_part[0].len()];
    for k in 0..split_part.len() {
        if tomerge.contains(&k) {
            for (q, hap) in split_part[k].iter().enumerate() {
                for read in hap {
                    new_part[q].insert(*read);
                }
            }
            for (q, hap) in split_part[k + 1].iter().enumerate() {
                for read in hap {
                    new_part[q].insert(*read);
                }
            }
            new_part_empty = false;
        } else {
            if !new_part_empty {
                split_part_merge.push(new_part);
                new_part = vec![FxHashSet::default(); split_part[0].len()];
                new_part_empty = true;
            } else {
                split_part_merge.push(mem::take(&mut split_part[k]));
            }
        }
    }
    if split_part_merge.len() > 1 {
        //        dbg!(&break_pos, split_part_merge.len());
        //        dbg!(&snp_breakpoints);
        for part in split_part_merge.iter() {
            log::trace!("---------------PART---------------");
            log::trace!(
                "BREAKPOS {:?}, ORIG {:?}",
                &break_pos,
                &original_snp_endpoints
            );
            for (s, hap) in part.iter().enumerate() {
                log::trace!("---------------HAP {} ---------------", s);
                log::trace!("{}", hap.len());
                //                        for read in hap{
                //                            println!("{}-{}",read.first_position, read.last_position);
                //                        }
            }
        }
    }
    //TODO
    //    let snp_breakpoints =
    //        vec![(original_snp_endpoints.0, original_snp_endpoints.1); snp_breakpoints.len()];
    return (split_part_merge, snp_breakpoints);
}

pub fn get_disjoint_paths_rewrite<'a>(
    hap_graph: &'a mut Vec<Vec<HapNode>>,
    flow_update_vec: FlowUpVec,
    floria_out_dir: String,
    vcf_profile: &VcfProfile,
    contig: &str,
    options: &Options,
) -> (Vec<FxHashSet<&'a Frag>>, Vec<(SnpPosition, SnpPosition)>) {
    let block_len = options.block_length;
    let do_binning = options.do_binning;
    let mut hap_petgraph = StableGraph::<(usize, usize), f64>::new();
    //Update the graph to include flows.
    for (n1_inf, n2_inf, flow) in flow_update_vec {
        //if flow < constants::FLOW_CUTOFF_MULT * options.epsilon {
        if flow < constants::MIN_SHARED_READS_UNAMBIG{
            continue;
        }
        hap_graph[n1_inf.0][n1_inf.1]
            .out_flows
            .push((n2_inf.1, flow));
    }
    //Populate the traceback vector and the container for the NodeIndex
    let mut node_index_container = vec![];
    for block in hap_graph.iter() {
        let mut node_index_block = vec![];
        for node in block.iter() {
            node_index_block.push(hap_petgraph.add_node((node.column, node.row)));
        }
        node_index_container.push(node_index_block);
    }

    //Add edges using NodeIndices
    for block in hap_graph.iter() {
        for node in block.iter() {
            let node1 = &node_index_container[node.column][node.row];
            for (next_row, flow) in node.out_flows.iter() {
                let node2 = &node_index_container[node.column + 1][*next_row];
                hap_petgraph.add_edge(*node1, *node2, *flow);
            }
        }
    }

    let num_starting_nodes = hap_petgraph.node_count();
    let mut trace_back_vec = vec![
        TraceBackNode {
            score: 0.,
            prev_ind: None,
            is_sink: false,
            is_source: false
        };
        num_starting_nodes
    ];
    for node_index in hap_petgraph.node_indices() {
        let in_neigh: Vec<_> = hap_petgraph
            .neighbors_directed(node_index, Direction::Incoming)
            .into_iter()
            .collect();
        let out_neigh: Vec<_> = hap_petgraph
            .neighbors_directed(node_index, Direction::Outgoing)
            .into_iter()
            .collect();
        let is_source = in_neigh.is_empty();
        let is_sink = out_neigh.is_empty();
        let score;
        if is_source {
            score = f64::MAX;
        } else {
            score = 0.;
        }
        trace_back_vec[node_index.index()] = TraceBackNode {
            score: score,
            prev_ind: None,
            is_sink: is_sink,
            is_source: is_source,
        };
    }

    if log_enabled!(Debug) || log_enabled!(Trace){
        let mut pet_graph_file =
            File::create(format!("{}/pet_graph.dot", floria_out_dir)).expect("Can't create file");
        write!(pet_graph_file, "{:?}", Dot::new(&hap_petgraph)).unwrap();
    }
    let mut iter_count = 0;
    let mut all_joined_path_parts = vec![];
    let mut cov_of_haplogroups = vec![];
    let mut path_parts_snp_endspoints = vec![];
    let mut best_paths = vec![];
    let mut best_pathscolrow = vec![];

    while hap_petgraph.node_count() > 0 {
        if iter_count != 0 {
            trace_back_vec = vec![
                TraceBackNode {
                    score: 0.,
                    prev_ind: None,
                    is_sink: false,
                    is_source: false
                };
                num_starting_nodes
            ];
            for node_index in hap_petgraph.node_indices() {
                let in_neigh: Vec<_> = hap_petgraph
                    .neighbors_directed(node_index, Direction::Incoming)
                    .into_iter()
                    .collect();
                let out_neigh: Vec<_> = hap_petgraph
                    .neighbors_directed(node_index, Direction::Outgoing)
                    .into_iter()
                    .collect();
                let is_source = in_neigh.is_empty();
                let is_sink = out_neigh.is_empty();
                let score;
                if is_source {
                    score = f64::MAX;
                } else {
                    score = 0.;
                }
                trace_back_vec[node_index.index()] = TraceBackNode {
                    score: score,
                    prev_ind: None,
                    is_sink: is_sink,
                    is_source: is_source,
                };
            }
        }
        //Top sort and find the maximal path
        let top_order = algo::toposort(&hap_petgraph, None).unwrap();
        let mut flow_cut_edges = vec![];
        for node_index in top_order {
            let out_edges = hap_petgraph.edges(node_index);
            for edge in out_edges {
                let source = edge.source();
                let target = edge.target();
                let flow = edge.weight();
                //The second condition means that the new flow has to be at least a third
                //of the smallest flow for the previous path: large dropoff indicates
                //that the main strain for the next node is diff. than previous node.
                if f64::min(trace_back_vec[source.index()].score, *flow)
                    > trace_back_vec[target.index()].score
                {
                    if *flow < trace_back_vec[source.index()].score * 0.33
                        && !trace_back_vec[source.index()].is_source
                    {
                        //Also cut off the edge from main strain to low cov strain if
                        //the in-node has only one high cov in-edge
                        //
                        //              O
                        //              |  100
                        //              O
                        //          90 / \  10
                        //            O   O
                        //
                        //            Cut off 10-edge
                        let in_neigh_source: Vec<_> = hap_petgraph
                            .neighbors_directed(source, Direction::Incoming)
                            .into_iter()
                            .collect();
                        let in_neigh_target: Vec<_> = hap_petgraph
                            .neighbors_directed(target, Direction::Incoming)
                            .into_iter()
                            .collect();

                        if in_neigh_source.len() == 1 {
                            flow_cut_edges.push(edge.id());
                            //                                flow_graph_versatile.remove_edge(edge.id());
                        }
                        if in_neigh_target.len() == 1 {
                            trace_back_vec[target.index()].score = f64::MAX;
                            trace_back_vec[target.index()].is_source = true;
                        }
                    } else {
                        trace_back_vec[target.index()].score =
                            f64::min(trace_back_vec[source.index()].score, *flow);
                        trace_back_vec[target.index()].prev_ind = Some(source.index());
                    }
                }
            }
        }

        for edge in flow_cut_edges {
            hap_petgraph.remove_edge(edge);
        }

        //Get best paths
        let mut index_of_best_end_node = None;
        let mut best_score = f64::MIN;
        let mut best_path = vec![];
        let mut best_path_colrow = vec![];
        for (i, trace_back_node) in trace_back_vec.iter().enumerate() {
            if trace_back_node.score > best_score && trace_back_node.is_sink {
                index_of_best_end_node = Some(i);
                best_score = trace_back_node.score;
            }
        }
        if let None = index_of_best_end_node {
            dbg!(&hap_petgraph);
            dbg!(&trace_back_vec.iter().enumerate());
            panic!("Shouldn't get here");
        }
        log::trace!(
            "Iter {}, Node count {}",
            iter_count,
            hap_petgraph.node_count()
        );
        let mut joined_path_part = FxHashSet::default();
        let mut snp_endpoints = (SnpPosition::MAX, SnpPosition::MIN);
        let mut haplogroup_flows = vec![];
        while !index_of_best_end_node.is_none() {
            let node_index = NodeIndex::new(index_of_best_end_node.unwrap());
            let out_edges = hap_petgraph.edges(node_index);
            for edge in out_edges {
                let flow = edge.weight();
                haplogroup_flows.push(*flow);
            }
            if let None = hap_petgraph.node_weight(node_index) {
                dbg!(&hap_petgraph, node_index, &trace_back_vec.len());
                panic!();
            }
            let (col, row) = hap_petgraph.node_weight(node_index).unwrap();
            let hap_graph_node = &hap_graph[*col][*row];

            if hap_graph_node.snp_endpoints.0 < snp_endpoints.0 {
                snp_endpoints.0 = hap_graph_node.snp_endpoints.0;
            }
            if hap_graph_node.snp_endpoints.1 > snp_endpoints.1 {
                snp_endpoints.1 = hap_graph_node.snp_endpoints.1;
            }

            for frag in hap_graph_node.frag_set.iter() {
                joined_path_part.insert(*frag);
            }

            //            dbg!(index_of_best_end_node,trace_back_vec[index_of_best_end_node.0][index_of_best_end_node.1]);
            best_path.push(index_of_best_end_node);
            best_path_colrow.push((col.clone(), row.clone()));
            log::trace!("{:?}, {:?}", &index_of_best_end_node, snp_endpoints);
            index_of_best_end_node = trace_back_vec[index_of_best_end_node.unwrap()].prev_ind;
        }
        let haplogroup_flow: Option<f64>;
        if haplogroup_flows.is_empty() {
            haplogroup_flow = None;
        } else {
            haplogroup_flow =
                Some(haplogroup_flows.iter().sum::<f64>() / haplogroup_flows.len() as f64);
        }

        for index in best_path.iter() {
            let node_index = NodeIndex::new(index.unwrap());
            hap_petgraph.remove_node(node_index);
        }

        all_joined_path_parts.push(joined_path_part);
        path_parts_snp_endspoints.push(snp_endpoints);

        iter_count += 1;

        best_paths.push(best_path);
        best_pathscolrow.push(best_path_colrow);
        cov_of_haplogroups.push(haplogroup_flow);
    }

    log::debug!("Number of haplogroups/disjoint paths: {}", best_paths.len());
//    let glopp_out_dir_copy = glopp_out_dir.clone();
//    let mut path_debug_file =
//        File::create(format!("{}/debug_paths.txt", glopp_out_dir_copy)).expect("Can't create file");
//    for (i, path) in best_pathscolrow.iter().enumerate() {
//        writeln!(path_debug_file, "{}", i).unwrap();
//        writeln!(path_debug_file, "{:?}", path).unwrap();
//        //        writeln!(path_debug_file, "{:?}", path_parts_snps_endpoints_copy[i]).unwrap();
//        writeln!(path_debug_file, "{:?}", cov_of_haplogroups[i]).unwrap();
//    }

    //Put read into best haplotig.
    if do_binning {
        let (binned_path_parts_snp_endspoints, binned_all_joined_path_parts) =
            part_block_manip::bin_haplogroups(
                &all_joined_path_parts,
                &path_parts_snp_endspoints,
                &cov_of_haplogroups,
                &vcf_profile,
                contig,
                block_len,
            );
        all_joined_path_parts = binned_all_joined_path_parts;
        path_parts_snp_endspoints = binned_path_parts_snp_endspoints;
    }

    return (all_joined_path_parts, path_parts_snp_endspoints);
}
