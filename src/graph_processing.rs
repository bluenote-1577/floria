use crate::file_reader;
use crate::global_clustering;
use crate::local_clustering;
use crate::types_structs::{Frag, HapNode, TraceBackNode, VcfProfile};
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use highs::{RowProblem, Sense};
use rayon::prelude::*;
use std::sync::Mutex;
use std::time::Instant;
//use osqp::{CscMatrix, Problem, Settings};
use petgraph::algo;
use petgraph::dot::Dot;
use petgraph::prelude::*;
use std::fs::File;
use std::io::Write;
use std::mem;

//Not using osqp anymore.

//pub fn entry_list_to_csc<'a>(
//    mut entry_list: Vec<(usize, usize, f64)>,
//    ncols: usize,
//    nrows: usize,
//) -> CscMatrix<'a> {
//    //Traverse each column first
//
//    let mut indptr = vec![];
//    let mut indices = vec![];
//    let mut data = vec![];
//    entry_list.sort_by(|x, y| x.0.cmp(&y.0));
//    let mut current_col = usize::MAX;
//    for entry in entry_list {
//        let col_id = entry.0;
//        if col_id != current_col {
//            if current_col == usize::MAX {
//                current_col = 0;
//                indptr.push(0);
//            }
//            for _i in current_col..col_id {
//                indptr.push(data.len());
//            }
//            current_col = col_id;
//        }
//        let row_id = entry.1;
//        let val = entry.2;
//        indices.push(row_id);
//        data.push(val);
//    }
//    indptr.push(data.len());
//
//    CscMatrix {
//        nrows: nrows,
//        ncols: ncols,
//        indptr: indptr.into(),
//        indices: indices.into(),
//        data: data.into(),
//    }
//}

type FlowUpVec = Vec<((usize, usize), (usize, usize), f64)>;

pub fn update_hap_graph(hap_graph: &mut Vec<Vec<HapNode>>) {
    //    let pseudo_count = 10.;
    let cutoff_val = 3.0;
    let mut out_edges_block_hap = vec![];
    for i in 0..hap_graph.len() - 1 {
        let mut out_edges_block = vec![];
        let hap_block1 = &hap_graph[i];
        let hap_block2 = &hap_graph[i + 1];
        for hap_node1 in hap_block1 {
            let mut out_weights = vec![0.0; hap_block2.len()];
            //We have to make sure that ambiguous reads do not
            //get accounted for.
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
                if out_weights[l] > cutoff_val {
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

pub fn solve_lp_graph(hap_graph: &Vec<Vec<HapNode>>, glopp_out_dir: String) -> FlowUpVec {
    let flow_cutoff = 3.0;
    let mut ae = vec![];

    //LP values
    let mut pb = RowProblem::default();
    let mut t = vec![];
    let mut x = vec![];

    //QP variables
    // Vector looks like [x,t]
    let mut num_constraints = 0;
    let mut a_entry = vec![];
    let mut p_entry = vec![];
    let mut l = vec![];
    let mut u = vec![];

    let mut hap_graph_vec = vec![];
    for hap_block in hap_graph.iter() {
        for hap_node in hap_block.iter() {
            hap_graph_vec.push(hap_node);
        }
    }

    let mut edge_to_nodes = vec![];
    let mut nodes_to_edges = FxHashMap::default();

    for hap_node in hap_graph_vec.iter() {
        let id1 = hap_node.id;
        for edge in hap_node.out_edges.iter() {
            let id2 = hap_graph[hap_node.column + 1][edge.0].id;
            edge_to_nodes.push((id1, id2));
            nodes_to_edges.insert((id1, id2), edge_to_nodes.len() - 1);
            //            ae.push(edge.1 * hap_node.cov());
            ae.push(edge.1);
        }
    }

    for _i in 0..edge_to_nodes.len() {
        x.push(pb.add_column(0., 0..));
    }
    for _i in 0..edge_to_nodes.len() {
        t.push(pb.add_column(1., 0..));
    }

    for i in edge_to_nodes.len()..2 * edge_to_nodes.len() {
        p_entry.push((i, i, 1. / ae[i - edge_to_nodes.len()]));
    }

    for (column_ind, hap_block) in hap_graph.iter().enumerate() {
        if column_ind == 0 || column_ind == hap_graph.len() - 1 {
            continue;
        }
        for hap_node in hap_block.iter() {
            if !hap_node.in_edges.is_empty() && !hap_node.out_edges.is_empty() {
                let node_id = hap_node.id;
                let mut in_edge_ids = vec![];
                let mut out_edge_ids = vec![];
                for in_edge in hap_node.in_edges.iter() {
                    let node_id2 = hap_graph[hap_node.column - 1][in_edge.0].id;
                    let in_edge_id = nodes_to_edges.get(&(node_id2, node_id)).unwrap();
                    in_edge_ids.push(in_edge_id);
                }
                for out_edge in hap_node.out_edges.iter() {
                    let node_id2 = hap_graph[hap_node.column + 1][out_edge.0].id;
                    let out_edge_id = nodes_to_edges.get(&(node_id, node_id2)).unwrap();
                    out_edge_ids.push(out_edge_id);
                }
                let mut constraint_row = vec![];

                for in_edge_id in in_edge_ids.iter() {
                    constraint_row.push((x[**in_edge_id], 1.));
                }
                for out_edge_id in out_edge_ids.iter() {
                    constraint_row.push((x[**out_edge_id], -1.));
                }

                //                dbg!(&constraint_row, &hap_node.in_edges);
                pb.add_row(..0, &constraint_row);
                pb.add_row(0.., &constraint_row);

                //QP variables
                for in_edge_id in in_edge_ids {
                    a_entry.push((*in_edge_id, num_constraints, 1.));
                }
                for out_edge_id in out_edge_ids {
                    a_entry.push((*out_edge_id, num_constraints, -1.));
                }
                l.push(0.);
                u.push(0.);
                num_constraints += 1;
            }
        }
    }

    for i in 0..t.len() {
        pb.add_row(-1.0 * ae[i].., &[(t[i], 1.), (x[i], -1.)]);
        pb.add_row(1.0 * ae[i].., &[(t[i], 1.), (x[i], 1.)]);
        pb.add_row(0.0.., &[(x[i], 1.)]);
    }

    //QP variables
    for i in 0..t.len() {
        a_entry.push((i, num_constraints, 1.));
        a_entry.push((i + t.len(), num_constraints, -1.));
        l.push(ae[i]);
        u.push(ae[i]);
        num_constraints += 1;

        a_entry.push((i, num_constraints, 1.));
        l.push(0.);
        u.push(f64::MAX);
        num_constraints += 1;
    }

    //QP SOLVING
    //    let p = entry_list_to_csc(p_entry, t.len() * 2, t.len() * 2);
    //    let a = entry_list_to_csc(a_entry, t.len() * 2, num_constraints);
    //    let mut q = vec![0.; t.len() * 2];
    //    let lambda = 0.;
    //    for i in 0..x.len() {
    //        q[i] = lambda;
    //    }
    //
    //    // Disable verbose output
    //    let settings = Settings::default().verbose(false);
    //
    //    // Create an OSQP problem
    //    let mut prob = Problem::new(p, &q, a, &l, &u, &settings).expect("failed to setup problem");
    //
    //    // Solve problem
    //    let qp_solution = prob.solve();
    //    println!("{:?}", qp_solution.x().expect("failed to solve problem"));

    let solved = pb.optimise(Sense::Minimise).solve();
    let solution = solved.get_solution();

    let mut file = File::create(format!("{}/graph.csv", glopp_out_dir)).expect("Can't create file");
    for i in 0..edge_to_nodes.len() {
        if solution.columns()[i] < flow_cutoff {
            continue;
        }
        writeln!(
            file,
            "{},{}-{},{}-{},{}",
            &solution.columns()[i],
            hap_graph_vec[edge_to_nodes[i].0].column,
            edge_to_nodes[i].0,
            hap_graph_vec[edge_to_nodes[i].1].column,
            edge_to_nodes[i].1,
            ae[i]
        )
        .unwrap();
    }
    drop(file);

    //    let mut file = File::create(format!("{}/qp_graph.csv", glopp_out_dir)).expect("Can't create file");
    //    for i in 0..edge_to_nodes.len() {
    //        if qp_solution.x().unwrap()[i] < flow_cutoff {
    //            continue;
    //        }
    //        writeln!(
    //            file,
    //            "{},{}-{},{}-{},{}",
    //            &qp_solution.x().unwrap()[i],
    //            hap_graph_vec[edge_to_nodes[i].0].column,
    //            edge_to_nodes[i].0,
    //            hap_graph_vec[edge_to_nodes[i].1].column,
    //            edge_to_nodes[i].1,
    //            ae[i]
    //        )
    //        .unwrap();
    //    }
    //    drop(file);

    let mut flow_update_vec = vec![];
    for i in 0..edge_to_nodes.len() {
        let (node1_id, node2_id) = edge_to_nodes[i];
        let node1 = hap_graph_vec[node1_id];
        let node2 = hap_graph_vec[node2_id];
        let flow = solution.columns()[i];
        flow_update_vec.push(((node1.column, node1.row), (node2.column, node2.row), flow));
    }

    println!("Linear program finished.");
    return flow_update_vec;
}

fn get_local_hap_blocks<'a>(
    _num_blocks: usize,
    _num_iters: usize,
    all_frags: &'a Vec<Frag>,
    epsilon: f64,
    snp_to_genome_pos: &'a Vec<usize>,
    max_number_solns: usize,
    _block_length: usize,
    glopp_out_dir: &str,
    j: usize,
    random_vec: &Vec<(usize, usize)>,
) -> Vec<Vec<HapNode<'a>>> {
    let ploidy_start = 1;
    let ploidy_end = 6;
    let error_rate = epsilon;
    let num_ploidies = ploidy_end - ploidy_start;
    let mut mec_vector = vec![0.; num_ploidies];
    let mut parts_vector = vec![];
    let mut expected_errors_ref = vec![];
    let mut endpoints_vector = vec![];
    // NOTE THE 1 INDEXING!
    let reads = local_clustering::find_reads_in_interval(
        random_vec[j].0 + 1,
        random_vec[j].1 + 1,
        all_frags,
        usize::MAX,
    );
    let mut best_ploidy = ploidy_start;
    if reads.is_empty() {
        return vec![];
    }
    for ploidy in ploidy_start..ploidy_end {
        best_ploidy = ploidy;
        let mut num_alleles = 0.0;
        let mut vec_reads_own = vec![];
        for read in reads.iter() {
            vec_reads_own.push(*read);
        }
        vec_reads_own.sort_by(|a, b| a.first_position.cmp(&b.first_position));
        let (break_pos, part) = global_clustering::beam_search_phasing(
            vec![FxHashSet::default(); ploidy],
            &vec_reads_own,
            epsilon,
            0.05,
            //            f64::MIN,
            0.0001_f64.ln(),
            max_number_solns,
            true,
            false,
        );

        //            let optimized_part = part;
        let (_new_score, optimized_part, block) =
            local_clustering::optimize_clustering(part, epsilon, 20);

        let split_part =
            utils_frags::split_part_using_breaks(&break_pos, &optimized_part, &all_frags);
        let endpoints;
        if j != 0 {
            endpoints = (random_vec[j - 1].1 + 1, random_vec[j].1 + 1);
        } else {
            endpoints = (random_vec[j].0 + 1, random_vec[j].1 + 1);
        }
        let (split_part_merge, split_part_endpoints) =
            merge_split_parts(split_part, break_pos, endpoints);
        //            let block = utils_frags::hap_block_from_partition(&optimized_part);
        //            let (binom_vec, _freq_vec) = local_clustering::get_partition_stats(&part, &block);
        let binom_vec = local_clustering::get_mec_stats_epsilon(&optimized_part, &block, epsilon, false);
        for (good, bad) in binom_vec {
            mec_vector[ploidy - ploidy_start] += bad;
            num_alleles += good;
            num_alleles += bad;
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

        expected_errors_ref.push(num_alleles as f64 * error_rate);

        if ploidy > ploidy_start {
            //                let mec_threshold = 1.0 / (1.0 - error_rate) / (1.0 + 1.0 / (ploidy + 1) as f64);
            //            log::trace!(
            //                "MEC vector {:?}, error_thresh {:?}",
            //                &mec_vector,
            //                expected_errors_ref
            //            );
            let mec_threshold =
                1.0 / (1.0 - error_rate) / (1.0 + 1.0 / ((ploidy as f64).powf(0.75) + 1.32) as f64);
            //            log::trace!(
            //                "Expected MEC ratio {}, observed MEC ratio {}",
            //                mec_threshold,
            //                mec_vector[ploidy - ploidy_start] as f64
            //                    / mec_vector[ploidy - ploidy_start - 1] as f64
            //            );
            if (mec_vector[ploidy - ploidy_start] as f64
                / mec_vector[ploidy - ploidy_start - 1] as f64)
                < mec_threshold
            {
            } else {
                log::trace!("MEC decrease thereshold, returning ploidy {}.", ploidy - 1);
                best_ploidy -= 1;
                break;
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

    log::trace!("DIFF\t{}\t{}", mec_vector[0], mec_vector[1]);

    log::trace!(
        "MEC vector {:?}, error_thresh {:?}, SNPs interval  {} {}",
        &mec_vector,
        expected_errors_ref,
        random_vec[j].0 + 1,
        random_vec[j].1 + 1,
    );

    let best_parts = mem::take(&mut parts_vector[best_ploidy - ploidy_start]);
    let best_endpoints = mem::take(&mut endpoints_vector[best_ploidy - ploidy_start]);
    let local_part_dir = format!("{}/local_parts/", glopp_out_dir);
    let mut hap_node_blocks = vec![];

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

        file_reader::write_output_partition_to_file(
            &frag_best_part,
            &vec![],
            local_part_dir.clone(),
            &format!("{}-{}-{}-{}", j, l, random_vec[j].0, best_ploidy),
            &snp_to_genome_pos,
            false
        );
    }

    return hap_node_blocks;
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
    num_blocks: usize,
    num_iters: usize,
    all_frags: &'a Vec<Frag>,
    epsilon: f64,
    snp_to_genome_pos: &'a Vec<usize>,
    max_number_solns: usize,
    block_length: usize,
    glopp_out_dir: String,
    minimal_density: f64
) -> Vec<Vec<HapNode<'a>>> {
    let using_bam;
    //Using frags instead of bam
    if snp_to_genome_pos.len() == 0 {
        using_bam = false;
    } else {
        using_bam = true;
    }

    let mut iter_vec: Vec<(usize, usize)> = vec![];
    if using_bam == false {
        let temp_iter_vec: Vec<usize> = (0..num_blocks).step_by(num_blocks / num_iters).collect();
        for i in 0..temp_iter_vec.len() - 1 {
            iter_vec.push((temp_iter_vec[i], temp_iter_vec[i + 1]));
        }
    } else {
        //        iter_vec = (0..num_blocks).step_by(num_blocks / num_iters).collect();
        iter_vec =
            utils_frags::get_range_with_lengths(snp_to_genome_pos, block_length, block_length / 3, minimal_density);
    }

    let interval_vec = iter_vec[0..iter_vec.len()].to_vec();
    log::trace!("SNP Endpoints {:?}", &interval_vec);

    let block_chunks: Mutex<Vec<_>> = Mutex::new(vec![]);
    (0..interval_vec.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|j| {
            //    for j in 0..random_vec.len() {
            //            if j % 1 == 0 {
            //                println!(
            //                    "Iteration {}/{}, SNP coords {} ",
            //                    j,
            //                    random_vec.len(),
            //                    random_vec[j].1
            //                );
            //            }

            let block_chunk = get_local_hap_blocks(
                num_blocks,
                num_iters,
                all_frags,
                epsilon,
                snp_to_genome_pos,
                max_number_solns,
                block_length,
                &glopp_out_dir,
                j,
                &interval_vec,
            );

            let mut locked = block_chunks.lock().unwrap();
            locked.push((j, block_chunk));
        });

    let block_chunks = block_chunks.into_inner().unwrap();
    let mut hap_node_blocks = process_chunks(block_chunks);
    println!("Phasing done");
    update_hap_graph(&mut hap_node_blocks);
    hap_node_blocks
}

fn merge_split_parts(
    mut split_part: Vec<Vec<FxHashSet<&Frag>>>,
    break_pos: FxHashMap<usize, FxHashSet<usize>>,
    original_snp_endpoints: (usize, usize),
) -> (Vec<Vec<FxHashSet<&Frag>>>, Vec<(usize, usize)>) {
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

pub fn get_disjoint_paths_rewrite(
    hap_graph: &mut Vec<Vec<HapNode>>,
    flow_update_vec: FlowUpVec,
    epsilon: f64,
    glopp_out_dir: String,
    snp_to_genome_pos: &Vec<usize>,
    short_frags: &Vec<Frag>,
    reassign_short: bool,
    vcf_profile: &VcfProfile,
    contig: &str,
    block_len: usize,
    do_binning: bool,
    extend_read_clipping: bool
) {
    let flow_cutoff = 3.0;
    let mut hap_petgraph = StableGraph::<(usize, usize), f64>::new();
    //Update the graph to include flows.
    for (n1_inf, n2_inf, flow) in flow_update_vec {
        if flow < flow_cutoff {
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

    let mut pet_graph_file =
        File::create(format!("{}/pet_graph.dot", glopp_out_dir)).expect("Can't create file");
    write!(pet_graph_file, "{:?}", Dot::new(&hap_petgraph)).unwrap();
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
        let mut snp_endpoints = (usize::MAX, usize::MIN);
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

    println!("Number of haplogroups/disjoint paths: {}", best_paths.len());

    //Put read into best haplotig.
    if do_binning{
        let (binned_path_parts_snp_endspoints, binned_all_joined_path_parts) = bin_haplogroups(
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

    process_reads_for_final_parts(
        &mut all_joined_path_parts,
        epsilon,
        short_frags,
        &path_parts_snp_endspoints,
        reassign_short,
    );

    
    let glopp_out_dir_copy = glopp_out_dir.clone();
    let path_parts_snps_endpoints_copy = path_parts_snp_endspoints.clone();
    file_reader::write_output_partition_to_file(
        &all_joined_path_parts,
        &path_parts_snp_endspoints,
        glopp_out_dir,
        &format!("all"),
        &snp_to_genome_pos,
        extend_read_clipping
    );

    let mut path_debug_file =
        File::create(format!("{}/debug_paths.txt", glopp_out_dir_copy)).expect("Can't create file");
    for (i, path) in best_pathscolrow.iter().enumerate() {
        writeln!(path_debug_file, "{}", i).unwrap();
        writeln!(path_debug_file, "{:?}", path).unwrap();
//        writeln!(path_debug_file, "{:?}", path_parts_snps_endpoints_copy[i]).unwrap();
        writeln!(path_debug_file, "{:?}", cov_of_haplogroups[i]).unwrap();
    }
}

fn process_reads_for_final_parts<'a>(
    all_joined_path_parts: &mut Vec<FxHashSet<&'a Frag>>,
    epsilon: f64,
    short_frags: &'a Vec<Frag>,
    snp_range_parts_vec: &Vec<(usize, usize)>,
    reassign_short: bool,
) {
    let mut all_parts_block = utils_frags::hap_block_from_partition(&all_joined_path_parts);
    let mut read_to_parts_map = FxHashMap::default();
    for (i, set) in all_joined_path_parts.iter().enumerate() {
        for frag in set.iter() {
            let corresponding_partitions = read_to_parts_map
                .entry(*frag)
                .or_insert(FxHashSet::default());
            corresponding_partitions.insert(i);
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
            diff_part_vec.push(((diff + 1.) / (same + 1.), id));
        }
        let best_part = diff_part_vec
            .iter()
            .min_by(|x, y| x.partial_cmp(&y).unwrap())
            .unwrap()
            .1;
        for id in part_ids.iter() {
            if *id != *best_part {
                all_joined_path_parts[*id].remove(frag);
                utils_frags::remove_read_from_block(&mut all_parts_block, frag, *id);
            }
        }
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

        println!("Time taken for reassign {:?}", Instant::now() - start_t);
    }
}

fn bin_haplogroups<'a>(
    parts: &Vec<FxHashSet<&'a Frag>>,
    snp_endpoints: &Vec<(usize, usize)>,
    cov_of_haplogroups: &Vec<Option<f64>>,
    vcf_profile: &VcfProfile,
    contig: &str,
    block_len: usize,
) -> (Vec<(usize, usize)>, Vec<FxHashSet<&'a Frag>>) {
    use statrs::distribution::{Discrete, Poisson};
    use statrs::statistics::Distribution;

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
    let gn_to_snp_pos = &vcf_profile.vcf_pos_to_snp_counter_map[contig];

    let mut clusters = vec![];
    let mut none_clusters = vec![];
    for i in 0..snp_endpoints.len() {
        let left_gn = snp_to_gn_pos[&(snp_endpoints[i].0 as i64)] as usize;
        let right_gn = snp_to_gn_pos[&(snp_endpoints[i].1 as i64)] as usize;
        let cov = cov_of_haplogroups[i];
        if !cov.is_none() {
            clusters.push(vec![(left_gn, right_gn, cov.unwrap(), i)]);
        }
        else{
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

            if h > i{
                min_move = 0;
            }

            else{
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
            if best_moves_i.len() == 1{
                best_moves.extend(best_moves_i);
            }
        }
        if best_moves.is_empty(){
            break;
        } else {
            best_moves.sort_by(|x,y| x.2.partial_cmp(&y.2).unwrap());
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
    for cluster in clusters{
        let mut new_snp_range = (usize::MAX,usize::MIN);
        let mut new_part = FxHashSet::default();
        for info in cluster{
            let index = info.3;
            let part = &parts[index];
            let range = &snp_endpoints[index];
            for frag in part.iter(){
                new_part.insert(*frag);
            }
            if range.0 < new_snp_range.0{
                new_snp_range.0 = range.0;
            }
            if range.1 > new_snp_range.1{
                new_snp_range.1 = range.1;
            }
        }
        new_parts.push(new_part);
        new_snp_ranges.push(new_snp_range);
    }

    for index in none_clusters{
        new_parts.push(parts[index].clone());
        new_snp_ranges.push(snp_endpoints[index].clone());
    }
    
    return ( new_snp_ranges, new_parts);
}
