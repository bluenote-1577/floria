use crate::file_reader;
use crate::global_clustering;
use crate::local_clustering;
use crate::types_structs::{Frag, HapNode};
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use highs::{RowProblem, Sense};
use osqp::{CscMatrix, Problem, Settings};
use petgraph::algo;
use petgraph::dot::{Config, Dot};
use petgraph::prelude::*;
use std::fs::File;
use std::io::Write;
use std::mem;

pub fn entry_list_to_csc<'a>(
    mut entry_list: Vec<(usize, usize, f64)>,
    ncols: usize,
    nrows: usize,
) -> CscMatrix<'a> {
    //Traverse each column first

    let mut indptr = vec![];
    let mut indices = vec![];
    let mut data = vec![];
    entry_list.sort_by(|x, y| x.0.cmp(&y.0));
    let mut current_col = usize::MAX;
    for entry in entry_list {
        let col_id = entry.0;
        if col_id != current_col {
            if current_col == usize::MAX {
                current_col = 0;
                indptr.push(0);
            }
            for i in current_col..col_id {
                indptr.push(data.len());
            }
            current_col = col_id;
        }
        let row_id = entry.1;
        let val = entry.2;
        indices.push(row_id);
        data.push(val);
    }
    indptr.push(data.len());

    CscMatrix {
        nrows: nrows,
        ncols: ncols,
        indptr: indptr.into(),
        indices: indices.into(),
        data: data.into(),
    }
}

type Flow_up_vec = Vec<((usize, usize), (usize, usize), f64)>;

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
            let sum: f64 = out_weights.iter().sum();
            let normalized_out_weights: Vec<f64> = out_weights.iter().map(|x| x / sum).collect();
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
                    println!("BLOCK {}: {}-{} weight {}", i, j, k, prob);
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

pub fn solve_lp_graph(hap_graph: &Vec<Vec<HapNode>>) -> Flow_up_vec {
    let flow_cutoff = 3.0;
    let mut ae = vec![];

    //LP values
    let mut pb = RowProblem::default();
    let mut t = vec![];
    let mut x = vec![];

    //QP variables
    // Vector looks like [x,t]
    let mut num_constraints = 0;
    let mut A_entry = vec![];
    let mut P_entry = vec![];
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
        P_entry.push((i, i, 1. / ae[i - edge_to_nodes.len()]));
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
                    A_entry.push((*in_edge_id, num_constraints, 1.));
                }
                for out_edge_id in out_edge_ids {
                    A_entry.push((*out_edge_id, num_constraints, -1.));
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
        A_entry.push((i, num_constraints, 1.));
        A_entry.push((i + t.len(), num_constraints, -1.));
        l.push(ae[i]);
        u.push(ae[i]);
        num_constraints += 1;

        A_entry.push((i, num_constraints, 1.));
        l.push(0.);
        u.push(f64::MAX);
        num_constraints += 1;
    }

    //QP SOLVING
    let P = entry_list_to_csc(P_entry, t.len() * 2, t.len() * 2);
    let A = entry_list_to_csc(A_entry, t.len() * 2, num_constraints);
    let mut q = vec![0.; t.len() * 2];
    let lambda = 0.;
    for i in 0..x.len() {
        q[i] = lambda;
    }

    // Disable verbose output
    let settings = Settings::default().verbose(false);

    // Create an OSQP problem
    let mut prob = Problem::new(P, &q, A, &l, &u, &settings).expect("failed to setup problem");

    // Solve problem
    let qp_solution = prob.solve();
    println!("{:?}", qp_solution.x().expect("failed to solve problem"));

    let solved = pb.optimise(Sense::Minimise).solve();
    let solution = solved.get_solution();

    let mut file = File::create("graph.csv").expect("Can't create file");
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

    let mut file = File::create("qp_graph.csv").expect("Can't create file");
    for i in 0..edge_to_nodes.len() {
        if qp_solution.x().unwrap()[i] < flow_cutoff {
            continue;
        }
        writeln!(
            file,
            "{},{}-{},{}-{},{}",
            &qp_solution.x().unwrap()[i],
            hap_graph_vec[edge_to_nodes[i].0].column,
            edge_to_nodes[i].0,
            hap_graph_vec[edge_to_nodes[i].1].column,
            edge_to_nodes[i].1,
            ae[i]
        )
        .unwrap();
    }
    drop(file);

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

pub fn get_disjoint_paths(hap_graph: &mut Vec<Vec<HapNode>>, flow_update_vec: Flow_up_vec) {
    //Max flow = true doesn't work right now
    let max_flow = false;
    let flow_cutoff = 3.0;
    let mut flow_graph_versatile = StableGraph::<(usize, usize), f64>::new();
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
    let mut trace_back_vec = vec![];
    for block in hap_graph.iter() {
        let mut trace_back_block = vec![];
        let mut node_index_block = vec![];
        for node in block.iter() {
            node_index_block.push(flow_graph_versatile.add_node((node.column, node.row)));
            let is_sink = node.out_flows.is_empty();
            let is_source = node.in_edges.is_empty();
            if max_flow {
                trace_back_block.push((0., (usize::MAX, usize::MAX), is_sink, is_source));
            } else {
                if is_source {
                    trace_back_block.push((
                        f64::MAX - 1.,
                        (usize::MAX, usize::MAX),
                        is_sink,
                        is_source,
                    ));
                } else {
                    trace_back_block.push((0., (usize::MAX, usize::MAX), is_sink, is_source));
                }
            }
        }
        node_index_container.push(node_index_block);
        trace_back_vec.push(trace_back_block);
    }

    //Add edges using NodeIndices
    for block in hap_graph.iter() {
        for node in block.iter() {
            let node1 = &node_index_container[node.column][node.row];
            for (next_row, flow) in node.out_flows.iter() {
                let node2 = &node_index_container[node.column + 1][*next_row];
                flow_graph_versatile.add_edge(*node1, *node2, *flow);
            }
        }
    }

    let mut pet_graph_file = File::create("pet_graph.dot").expect("Can't create file");
    write!(pet_graph_file, "{:?}", Dot::new(&flow_graph_versatile)).unwrap();
    let mut iter_count = 0;
    let mut all_joined_path_parts = vec![];

    while flow_graph_versatile.node_count() > 0 {
        //Top sort and find the maximal path
        let top_order = algo::toposort(&flow_graph_versatile, None).unwrap();
        for node_index in top_order {
            let out_edges = flow_graph_versatile.edges(node_index);
            for edge in out_edges {
                let source = flow_graph_versatile.node_weight(edge.source()).unwrap();
                let target = flow_graph_versatile.node_weight(edge.target()).unwrap();
                let flow = edge.weight();
                if max_flow {
                    if trace_back_vec[source.0][source.1].0 + flow
                        > trace_back_vec[target.0][target.1].0
                    {
                        trace_back_vec[target.0][target.1].0 =
                            trace_back_vec[source.0][source.1].0 + flow;
                        trace_back_vec[target.0][target.1].1 .0 = source.0;
                        trace_back_vec[target.0][target.1].1 .1 = source.1;
                    }
                } else {
                    if f64::min(trace_back_vec[source.0][source.1].0, *flow)
                        > trace_back_vec[target.0][target.1].0
                    {
                        trace_back_vec[target.0][target.1].0 =
                            f64::min(trace_back_vec[source.0][source.1].0, *flow);
                        trace_back_vec[target.0][target.1].1 .0 = source.0;
                        trace_back_vec[target.0][target.1].1 .1 = source.1;
                    }
                }
            }
        }

        //Get best paths
        let mut index_of_best_end_node = (usize::MAX, usize::MAX);
        let mut best_score = f64::MIN;
        let mut best_path = vec![];
        for (i, block) in trace_back_vec.iter().enumerate() {
            for (j, node) in block.iter().enumerate() {
                if node.0 > best_score && node.2 {
                    index_of_best_end_node = (i, j);
                    best_score = node.0;
                }
            }
        }
        if index_of_best_end_node == (usize::MAX, usize::MAX) {
            dbg!(&flow_graph_versatile);
            panic!("Shouldn't get here");
        }
        println!("Iter {}, Node count {}", iter_count, flow_graph_versatile.node_count());
        let mut joined_path_part = FxHashSet::default();
        while index_of_best_end_node != (usize::MAX, usize::MAX) {
            for frag in hap_graph[index_of_best_end_node.0][index_of_best_end_node.1]
                .frag_set
                .iter()
            {
                joined_path_part.insert(*frag);
            }
//            dbg!(index_of_best_end_node,trace_back_vec[index_of_best_end_node.0][index_of_best_end_node.1]);
            best_path.push(index_of_best_end_node);
            println!("{:?}", &index_of_best_end_node);
            index_of_best_end_node =
                trace_back_vec[index_of_best_end_node.0][index_of_best_end_node.1].1;
        }
        all_joined_path_parts.push(joined_path_part);

        for (col, row) in best_path {
                        //Reset the trace back vectors for nodes adjacent to the path
            let node_index = node_index_container[col][row];
            let out_neigh =
                flow_graph_versatile.neighbors_directed(node_index, Direction::Outgoing).into_iter().collect::<Vec<_>>();
            let in_neigh =
                flow_graph_versatile.neighbors_directed(node_index, Direction::Incoming).into_iter().collect::<Vec<_>>();
            //Remove the node
            flow_graph_versatile
                .remove_node(node_index)
                .unwrap();

            for neighbor in out_neigh {
                let target = flow_graph_versatile.node_weight(neighbor).unwrap();
                let is_source = flow_graph_versatile
                    .neighbors_directed(
                        node_index_container[target.0][target.1],
                        Direction::Incoming,
                    )
                    .into_iter()
                    .collect::<Vec<_>>()
                    .is_empty();
                let init_weight = if is_source { f64::MAX-1. } else { 0. };
                trace_back_vec[target.0][target.1] = (
                    init_weight,
                    (usize::MAX, usize::MAX),
                    //Removing in-nodes doesn't change sink-ness
                    trace_back_vec[target.0][target.1].2,
                    //May become a source now
                    is_source,
                    );
            }

            for neighbor in in_neigh {
                let target = flow_graph_versatile.node_weight(neighbor).unwrap();
                let is_source = trace_back_vec[target.0][target.1].3;
                let init_weight = if is_source { f64::MAX-1. } else { 0. };
                trace_back_vec[target.0][target.1] = (
                    init_weight,
                    (usize::MAX, usize::MAX),
                    //May become a sink now
                    flow_graph_versatile
                        .neighbors_directed(
                            node_index_container[target.0][target.1],
                            Direction::Outgoing,
                        )
                        .into_iter()
                        .collect::<Vec<_>>()
                        .is_empty(),
                    //Doesn't change source-ness
                    trace_back_vec[target.0][target.1].3,
                );
            }

            trace_back_vec[col][row].0 = f64::MIN;
            trace_back_vec[col][row].1 = (usize::MAX, usize::MAX);
        }
        iter_count += 1;
    }
    file_reader::write_output_partition_to_file(
        &all_joined_path_parts,
        "./parts_joined/",
        &format!("{}", iter_count),
        &FxHashMap::default(),
    );
}

pub fn get_best_paths(hap_graph: &mut Vec<Vec<HapNode>>, flow_update_vec: Flow_up_vec) {
    println!("Getting best paths");
    let sum_weight_max = true;
    let min_flow_max = false;
    let flow_cutoff = 3.0;
    for (n1_inf, n2_inf, flow) in flow_update_vec {
        if flow < flow_cutoff {
            continue;
        }
        hap_graph[n1_inf.0][n1_inf.1]
            .out_flows
            .push((n2_inf.1, flow));
    }
    let max_num_paths = 8;
    //TODO
    let run_min_flow_cutoff = 2;

    let mut cfile = File::create("color_info.csv").expect("Can't create file");

    let mut path_counter = 0;
    let mut already_seen_edges = FxHashMap::default();
    loop {
        let mut min_flow = f64::MAX;
        let mut run_min_flow = f64::MAX;
        let mut node_path_colrow = vec![];
        {
            let mut node_path = vec![];
            let mut trace_back_vec = vec![];

            for i in 0..hap_graph.len() {
                trace_back_vec.push(vec![(f64::MIN, usize::MAX); hap_graph[i].len()]);
            }

            if min_flow_max {
                for tup in trace_back_vec[0].iter_mut() {
                    tup.0 = f64::MAX;
                }
            } else {
                for tup in trace_back_vec[0].iter_mut() {
                    tup.0 = 0.;
                }
            }

            for i in 0..hap_graph.len() - 1 {
                for (j, node) in hap_graph[i].iter().enumerate() {
                    for (next_row_id, flow) in node.out_flows.iter() {
                        //Negative flows means the edge has been removed.
                        let mut mod_flow = *flow;
                        if *flow < flow_cutoff {
                            continue;
                        }
                        if already_seen_edges
                            .contains_key(&(node.id, hap_graph[i + 1][*next_row_id].id))
                        {
                            let mult = already_seen_edges
                                .get(&(node.id, hap_graph[i + 1][*next_row_id].id))
                                .unwrap();
                            //                           mod_flow = *flow/(1.5_f64.powf(*mult));
                            mod_flow = *flow;
                        }
                        if sum_weight_max {
                            if trace_back_vec[i][node.row].0 + mod_flow
                                > trace_back_vec[i + 1][*next_row_id].0
                            {
                                trace_back_vec[i + 1][*next_row_id].0 =
                                    trace_back_vec[i][node.row].0 + mod_flow;
                                trace_back_vec[i + 1][*next_row_id].1 = j;
                            }
                        } else {
                            if f64::min(mod_flow, trace_back_vec[i][node.row].0)
                                > trace_back_vec[i + 1][*next_row_id].0
                            {
                                trace_back_vec[i + 1][*next_row_id].0 =
                                    f64::min(mod_flow, trace_back_vec[i][node.row].0);
                                trace_back_vec[i + 1][*next_row_id].1 = j;
                            }
                        }
                    }
                }
            }

            let mut best_i = usize::MAX;
            let mut best_score = f64::MIN;
            let last_col = trace_back_vec.last().unwrap();
            for i in 0..last_col.len() {
                if last_col[i].0 > best_score {
                    best_i = i;
                    best_score = last_col[i].0;
                }
            }

            if best_score == f64::MIN {
                break;
            }

            let mut run_length = 0;
            for i in (0..hap_graph.len()).rev() {
                let best_node = &hap_graph[i][best_i];
                if best_node.out_flows.len() == 1 {
                    run_length += 1;
                    if run_length > run_min_flow_cutoff {
                        if best_node.out_flows[0].1 < run_min_flow {
                            run_min_flow = best_node.out_flows[0].1;
                        }
                    }
                } else {
                    run_length = 0;
                }
                let score = trace_back_vec[i][best_i].0;
                node_path.push(best_node);
                node_path_colrow.push((best_node.column, best_node.row));
                best_i = trace_back_vec[i][best_i].1;

                //TODO do min_flow properly, new avg flow termination conditoin
                if i != 0 {
                    if sum_weight_max {
                        let flow = score - trace_back_vec[i - 1][best_i].0;
                        if flow < min_flow {
                            min_flow = flow;
                        }
                    } else {
                        min_flow = best_score;
                    }
                }
            }

            node_path.reverse();
            node_path_colrow.reverse();
            for (l, node) in node_path.iter().enumerate() {
                if l == node_path.len() - 1 {
                    write!(cfile, "{}-{}\n", node.column, node.id).unwrap();
                } else {
                    write!(cfile, "{}-{},", node.column, node.id).unwrap();
                }
            }
            println!(
                "Highest weighted flow path is {}, with min flow {}, run min flow {}, avg coverage is {}",
                best_score, min_flow, run_min_flow, best_score/(node_path.len()-1) as f64
            );
        }
        //Update edges
        for i in 0..node_path_colrow.len() - 1 {
            let (col1, row1) = node_path_colrow[i];
            let (col2, row2) = node_path_colrow[i + 1];
            let id1 = hap_graph[col1][row1].id;
            let id2 = hap_graph[col2][row2].id;
            let edge_mult_count = already_seen_edges.entry((id1, id2)).or_insert(0.);
            *edge_mult_count += 1.;
        }
        for i in 0..node_path_colrow.len() - 1 {
            let (col, row) = node_path_colrow[i];
            let hap_node_mut = &mut hap_graph[col][row];
            for out_flow_edge in hap_node_mut.out_flows.iter_mut() {
                if out_flow_edge.0 == node_path_colrow[i + 1].1 {
                    if run_min_flow > 5000000.0 {
                        out_flow_edge.1 -= min_flow * 1.0;
                    } else {
                        out_flow_edge.1 -= run_min_flow * 1.0;
                    }
                    if out_flow_edge.1 < flow_cutoff {
                        out_flow_edge.1 = f64::MIN;
                    }
                }
            }
        }

        path_counter += 1;
        if path_counter > max_num_paths {
            break;
        }
    }
}

pub fn generate_hap_graph<'a>(
    num_blocks: usize,
    num_iters: usize,
    all_frags: &'a Vec<Frag>,
    epsilon: f64,
    snp_to_genome_pos: &'a Vec<usize>,
    max_number_solns: usize,
    block_length: usize,
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
            utils_frags::get_range_with_lengths(snp_to_genome_pos, block_length, block_length / 3);
    }

    let random_vec = iter_vec[0..iter_vec.len()].to_vec();
    dbg!(&random_vec);

    let ploidy_start = 1;
    let ploidy_end = 6;
    let error_rate = epsilon;
    let num_ploidies = ploidy_end - ploidy_start;
    let mut hap_node_blocks = vec![];
    let mut hap_node_counter = 0;
    let mut column_counter = 0;
    for j in 0..random_vec.len() {
        let mut mec_vector = vec![0.; num_ploidies];
        let mut parts_vector = vec![];
        let mut expected_errors_ref = vec![];
        // NOTE THE 1 INDEXING!
        let reads = local_clustering::find_reads_in_interval(
            random_vec[j].0 + 1,
            random_vec[j].1 + 1,
            all_frags,
            usize::MAX,
        );
        let mut best_ploidy = ploidy_start;
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
                0.0001_f64.ln(),
                max_number_solns,
                true,
                false,
                false,
            );

            //            let optimized_part = part;
            let (_new_score, optimized_part, block) =
                local_clustering::optimize_clustering(part, epsilon, 10);

            let split_part =
                utils_frags::split_part_using_breaks(&break_pos, &optimized_part, &all_frags);
            let split_part_merge = merge_split_parts(split_part, break_pos);
            //            let block = utils_frags::hap_block_from_partition(&optimized_part);
            //            let (binom_vec, _freq_vec) = local_clustering::get_partition_stats(&part, &block);
            let binom_vec =
                local_clustering::get_mec_stats_epsilon(&optimized_part, &block, epsilon);
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

            expected_errors_ref.push(num_alleles as f64 * error_rate);

            if ploidy > ploidy_start {
                //                let mec_threshold = 1.0 / (1.0 - error_rate) / (1.0 + 1.0 / (ploidy + 1) as f64);
                log::debug!(
                    "MEC vector {:?}, error_thresh {:?}",
                    &mec_vector,
                    expected_errors_ref
                );
                let mec_threshold = 1.0
                    / (1.0 - error_rate)
                    / (1.0 + 1.0 / ((ploidy as f64).powf(0.75) + 1.32) as f64);
                log::debug!(
                    "Expected MEC ratio {}, observed MEC ratio {}",
                    mec_threshold,
                    mec_vector[ploidy - ploidy_start] as f64
                        / mec_vector[ploidy - ploidy_start - 1] as f64
                );
                if (mec_vector[ploidy - ploidy_start] as f64
                    / mec_vector[ploidy - ploidy_start - 1] as f64)
                    < mec_threshold
                {
                } else {
                    println!("MEC decrease thereshold, returning ploidy {}.", ploidy - 1);
                    best_ploidy -= 1;
                    break;
                }
                if mec_vector[ploidy - ploidy_start] < expected_errors_ref[ploidy - ploidy_start] {
                    println!("MEC error threshold, returning ploidy {}.", ploidy);
                    break;
                }
            } else {
                if mec_vector[ploidy - ploidy_start] < expected_errors_ref[ploidy - ploidy_start] {
                    println!("MEC error threshold, returning ploidy {}.", ploidy);
                    break;
                }
            }
        }

        let best_parts = mem::take(&mut parts_vector[best_ploidy - ploidy_start]);
        for best_part in best_parts.iter() {
            let mut hap_node_block = vec![];
            let mut frag_best_part = vec![];

            for ind_part in best_part.iter() {
                let frag_set: FxHashSet<&Frag> = ind_part.iter().map(|x| &all_frags[*x]).collect();
                frag_best_part.push(frag_set.clone());
                let mut hap_node = HapNode::new(frag_set);
                hap_node.column = column_counter;
                hap_node.row = hap_node_block.len();
                hap_node.id = hap_node_counter;
                hap_node_counter += 1;
                hap_node_block.push(hap_node);
            }
            column_counter += 1;
            hap_node_blocks.push(hap_node_block);
            file_reader::write_output_partition_to_file(
                &frag_best_part,
                "./parts_graph/",
                &format!("{}-{}-{}", column_counter - 1, random_vec[j].0, best_ploidy),
                &FxHashMap::default(),
            );
            println!("Coords: {}", random_vec[j].0)
        }
    }
    update_hap_graph(&mut hap_node_blocks);
    hap_node_blocks
}

fn merge_split_parts(
    mut split_part: Vec<Vec<FxHashSet<&Frag>>>,
    break_pos: FxHashMap<usize, FxHashSet<usize>>,
) -> Vec<Vec<FxHashSet<&Frag>>> {
    let mut split_part_merge = vec![];
    let mut tomerge = vec![];
    for k in 0..split_part.len() - 1 {
        let total_cov_1: usize = split_part[k].iter().map(|x| x.len()).sum();
        let total_cov_2: usize = split_part[k + 1].iter().map(|x| x.len()).sum();
        let cov_rat = total_cov_1 as f64 / (total_cov_2 + total_cov_1) as f64;
        //TODO Merging stuff
        if cov_rat > 0.75 || cov_rat < 0.25 {
            //        if true {
            tomerge.push(k);
        }
    }

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
        dbg!(&break_pos, split_part_merge.len());
        for part in split_part_merge.iter() {
            println!("---------------PART---------------");
            for (s, hap) in part.iter().enumerate() {
                println!("---------------HAP {} ---------------", s);
                println!("{}", hap.len());
                //                        for read in hap{
                //                            println!("{}-{}",read.first_position, read.last_position);
                //                        }
            }
        }
    }
    return split_part_merge;
}
