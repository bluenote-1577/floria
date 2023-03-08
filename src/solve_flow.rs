use crate::types_structs::{HapNode, FlowUpVec};
//use std::fs::File;
//use crate::constants;
//use std::io::Write;
use fxhash::{FxHashMap};

#[cfg(feature = "highs")]
pub fn solve_lp_graph(hap_graph: &Vec<Vec<HapNode>>) -> FlowUpVec {
    use highs::{RowProblem, Sense};
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

//    let mut file = File::create(format!("{}/graph.csv", glopp_out_dir)).expect("Can't create file");
//    for i in 0..edge_to_nodes.len() {
//        writeln!(
//            file,
//            "{},{}-{},{},{}-{},{},{}",
//            &solution.columns()[i],
//            hap_graph_vec[edge_to_nodes[i].0].column,
//            hap_graph_vec[edge_to_nodes[i].0].row,
//            edge_to_nodes[i].0,
//            hap_graph_vec[edge_to_nodes[i].1].column,
//            hap_graph_vec[edge_to_nodes[i].1].row,
//            edge_to_nodes[i].1,
//            ae[i]
//        )
//        .unwrap();
//    }
//    drop(file);

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

    log::info!("Linear program finished.");
    return flow_update_vec;
}

#[cfg(not(feature = "highs"))]
pub fn solve_lp_graph(hap_graph: &Vec<Vec<HapNode>>) -> FlowUpVec{
    use minilp::{Problem, OptimizationDirection, ComparisonOp};

    // Maximize an objective function x + 2 * y of two variables x >= 0 and 0 <= y <= 3
    let mut problem = Problem::new(OptimizationDirection::Minimize);
    let mut ae = vec![];

    //LP values
    let mut t = vec![];
    let mut x = vec![];
    
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
//        x.push(pb.add_column(0., 0..));
        x.push(problem.add_var(0., (0., f64::INFINITY)));

    }
    for _i in 0..edge_to_nodes.len() {
        t.push(problem.add_var(1., (0., f64::INFINITY)));
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
//                pb.add_row(..0, &constraint_row);
//                pb.add_row(0.., &constraint_row);
                problem.add_constraint(&constraint_row, ComparisonOp::Eq, 0.0);

            }
        }
    }

    for i in 0..t.len() {
        problem.add_constraint(&[(t[i], 1.), (x[i], -1.)], ComparisonOp::Ge, -1.0 * ae[i]);
        problem.add_constraint(&[(t[i], 1.), (x[i], 1.)], ComparisonOp::Ge, 1.0 * ae[i]);
    }
    let solution = problem.solve().unwrap();
    let mut flow_update_vec = vec![];
    for (i,x_var) in x.iter().enumerate(){
        let (node1_id, node2_id) = edge_to_nodes[i];
        let node1 = hap_graph_vec[node1_id];
        let node2 = hap_graph_vec[node2_id];
        let flow = solution[*x_var];
        flow_update_vec.push(((node1.column, node1.row), (node2.column, node2.row), flow));
    }

    return flow_update_vec;

}
