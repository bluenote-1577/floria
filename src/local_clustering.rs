use crate::types_structs::{Frag, HapBlock};
use rand::prelude::*;
use rand_core::SeedableRng;
//use rand::rng::Rng;
use rand_pcg::Pcg64;
use statrs::distribution::{ChiSquared};
use statrs::distribution::ContinuousCDF;
extern crate time;
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};

#[cfg(not(debug_assertions))]
macro_rules! debug {
    ($x:expr) => {
        std::convert::identity($x)
    };
}

//use std::time::Instant;

//Return the set of reads for which every read covers at least one position in the interval
//(inclusive)
pub fn find_reads_in_interval<'a>(
    start: usize,
    end: usize,
    //position_to_reads : &FxHashMap<usize,FxHashSet<&Frag>>,
    all_frags: &'a Vec<Frag>,
    max_num_reads: usize
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
        if final_set.len() > max_num_reads{
            break;
        }
        if frag.last_position < start {
            continue;
        }
        if frag.first_position > end {
            break;
        }

        //TODO we use this routine in glopp estimate ploidy, don't want circular mappings.
        if frag.last_position - frag.first_position > 10000{
            continue;
        }

        //Currently we use 1/3 quantile as the length of the block, i.e.
        //end-start. If a mapping is weird and the fragment
        //spans several regions, we ignore the fragment.
        if false{
        //if frag.last_position - frag.first_position > 60 * (end - start) {
            continue;
        }

        final_set.insert(frag);
    }
    final_set
}

//Return a partition from a set of reads using our local clustering method.
pub fn generate_hap_block<'a>(
    start: usize,
    end: usize,
    //position_to_reads : &'a FxHashMap<usize,FxHashSet<&Frag>>,
    ploidy: usize,
    all_frags: &'a Vec<Frag>,
    epsilon: f64,
) -> Vec<FxHashSet<&'a Frag>> {
    //debug!(start);
    //debug!(end);
    let all_reads = find_reads_in_interval(start, end, all_frags, 100);
    let partition = cluster_reads(&all_reads, ploidy, epsilon);
    partition
}

//Compute distance between two reads from the precomputed hash map which can be useful for speeding
//up. We don't end up using this
//method right now because there is some weirdness with hashmaps being slow, I think I fixed this
//so maybe we'll reimplement in the future.
pub fn dist_from_graph(
    r1: &Frag,
    r2: &Frag,
    all_distances: &FxHashMap<&Frag, FxHashMap<&Frag, i32>>,
) -> Option<i32> {
    if !utils_frags::check_overlap(r1, r2) {
        None
    } else {
        let map1 = all_distances.get(r1).unwrap();
        let dist;

        if map1.contains_key(r2) {
            dist = map1.get(r2).unwrap();
        } else {
            dist = all_distances.get(r2).unwrap().get(r1).unwrap();
        }
        Some(*dist)
    }
}

//Local clustering method for a set of reads -> partition.
//Given a set of reads covering an interval, we greedily find a max k-clique and break this
//k-clique up iinto k different clusters. iteratively add reads to the best k-clique where the max
//of the intracluster distances is minimized.
//Importantly, the order in which we itertively add reads is sorted by the minimum of the maximum
//overlap of the read within the clusters.
pub fn cluster_reads<'a>(
    all_reads: &FxHashSet<&'a Frag>,
    ploidy: usize,
    epsilon: f64,
) -> Vec<FxHashSet<&'a Frag>> {
    let use_binomial_dist = true;

    //Get the largest distance edge
    let mut vec_all_edges = Vec::new();
    let vec_all_reads: Vec<_> = all_reads.iter().collect();
    //    dbg!(vec_all_reads.len());
    //    for read in vec_all_reads.iter(){
    //        dbg!(read.last_position - read.first_position,&read.id,&read.seq_dict.len(),&read.positions.len());
    //    }
    let mut adj_list_edges = Vec::new();
    for _i in 0..vec_all_reads.len() {
        adj_list_edges.push(Vec::new());
    }

    //Get local read-read graph, a.k.a the distance matrix. In the future, we can speed this up by precomputing a
    //
    //global read-read graph or precomputing the distances between reads.
    //
    //debug!("Computing read-graph");
    //debug!(vec_all_reads.len());
    for (i, r1) in vec_all_reads.iter().enumerate() {
        for j in i + 1..vec_all_reads.len() {
            //
            let r2 = &vec_all_reads[j];
            //Originally tried to use a precomputed read distance, this was much slower, however.
            //            let dist_try = dist_from_graph(*r1,*r2,&all_distances);
            //            let dist = match dist_try{
            //                Some(x) => x,
            //                None => continue,
            //            };
            //
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

    //debug!("Done computing read-graph");
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

    //DEBUGGING/TESTING
//        {
//            println!("CLIQUE VERTICES");
//            for vertex in used_vertices.iter() {
//                let id = &vec_all_reads[(*vertex) as usize].id;
//                println!("{}", id);
//                }
//                let mut split = id.split("/");
//                let cluster = split.next();
//                id_set.insert(cluster);
//            }
    
    //        if id_set.len() < ploidy {
    //            println!("clique partition not effective");
    //            for vertex in used_vertices.iter() {
    //                let frag1 = &vec_all_reads[(*vertex) as usize];
    //                for vertex in used_vertices.iter() {
    //                    let frag2 = &vec_all_reads[(*vertex) as usize];
    //                    println!(
    //                        "{},{},{},{}",
    //                        &frag1.id,
    //                        &frag2.id,
    //                        utils_frags::distance(frag1, frag2).1,
    //                        utils_frags::distance(frag1, frag2).0
    //                    );
    //                }
    //            }
    //        }
    //    }

    //Populate the clusters -- experimenting with diffferent methods...
    //    populate_clusters1(&mut clusters, &mut used_vertices, &vec_all_reads, &adj_list_edges, ploidy);
    //    populate_clusters2_binom(&mut clusters, &mut used_vertices, &vec_all_reads, &adj_list_edges, ploidy,10);
    populate_clusters3(
        &mut clusters,
        &mut used_vertices,
        &vec_all_reads,
        &adj_list_edges,
        ploidy,
        10,
    );

    //Turn the vertex indices into actual fragments -- could probably come up with a more elegant
    //solution using some sort of map function...
    let mut partition = Vec::new();
    for cluster in clusters.iter() {
        let mut frag_set = FxHashSet::default();
        for vertex in cluster.iter() {
            let vertex_usize = *vertex as usize;
            frag_set.insert(*vec_all_reads[vertex_usize]);
        }
        partition.push(frag_set)
    }
    partition
}

fn populate_clusters3(
    clusters: &mut Vec<FxHashSet<i32>>,
    used_vertices: &mut FxHashSet<i32>,
    vec_all_reads: &Vec<&&Frag>,
    adj_list_edges: &Vec<Vec<(f64, i32)>>,
    ploidy: usize,
    iters: usize,
) {
    //Sort reads(vertices) by the minimum of the maximum overlaps between clusters
    for iteration in 0..iters {
        let mut max_num_iters = vec_all_reads.len() / iters;
        if iteration == iters - 1 {
            max_num_iters = usize::MAX;
        }
        //First sort vertices by the minimum overlap for a read with all vertices in the clusters.
        //We do a loop here because once

        //Sort reads(vertices) by the minimum of the maximum overlaps between clusters
        let mut sorted_vec_overlap_reads = Vec::new();

        if used_vertices.len() == vec_all_reads.len() {
            break;
        }

        for (i, read) in vec_all_reads.iter().enumerate() {
            let index = i as i32;
            if used_vertices.contains(&index) {
                continue;
            }

            let edges = &adj_list_edges[i];
            let mut read_overlaps_between_clusters = Vec::new();
            for _i in 0..ploidy {
                read_overlaps_between_clusters.push(0);
            }

            for edge in edges.iter() {
                for (j, cluster) in clusters.iter().enumerate() {
                    if cluster.contains(&edge.1) {
                        let edge_index = edge.1 as usize;
                        let read2 = vec_all_reads[edge_index];
                        let t: Vec<_> = read.positions.intersection(&read2.positions).collect();
                        let overlap_len = t.len();
                        if overlap_len > read_overlaps_between_clusters[j] {
                            read_overlaps_between_clusters[j] = overlap_len;
                        }
                    }
                }
            }
            sorted_vec_overlap_reads.push((
                i,
                *read_overlaps_between_clusters.iter().min().unwrap() as usize,
            ));
        }

        //Obtained sorted vertices
        sorted_vec_overlap_reads.sort_by(|a, b| b.1.cmp(&a.1));

        //        if sorted_vec_overlap_reads.len() > 0{
        //            dbg!(&sorted_vec_overlap_reads[0], vec_all_reads[sorted_vec_overlap_reads[0].0 as usize].seq_dict.len(),iteration);
        //        }

        //Now greedily add vertices to the partition where the maximum distance to the cluster is
        //minimized.
        for (j, (vertex, _dist)) in sorted_vec_overlap_reads.iter().enumerate() {
            if j > max_num_iters {
                break;
            }
            let edges = &adj_list_edges[*vertex];
            let mut max_dist = Vec::new();
            for _i in 0..ploidy {
                max_dist.push(-1.0);
            }
            for edge in edges.iter() {
                for (j, cluster) in clusters.iter().enumerate() {
                    if cluster.contains(&edge.1) {
                        let dist = edge.0;
                        if max_dist[j] < dist {
                            max_dist[j] = dist;
                        }
                    }
                }
            }
            //            println!("{:?},{:?} max_dist, overlap",max_dist,overlap);

            //Find index of minimum distance cluster
            let mut min_index = 0;
            let mut min_score = f64::MAX;
            for i in 0..ploidy {
                if max_dist[i] < min_score {
                    min_index = i;
                    min_score = max_dist[i];
                }
            }

            let vertex_i32 = *vertex as i32;
            used_vertices.insert(vertex_i32);
            clusters[min_index].insert(vertex_i32);
        }
    }
}

fn _populate_clusters2_binom(
    clusters: &mut Vec<FxHashSet<i32>>,
    used_vertices: &mut FxHashSet<i32>,
    vec_all_reads: &Vec<&&Frag>,
    adj_list_edges: &Vec<Vec<(f64, i32)>>,
    ploidy: usize,
    iters: usize,
) {
    //Sort reads(vertices) by the minimum of the maximum overlaps between clusters
    for iteration in 0..iters {
        let mut sorted_vec_overlap_reads = Vec::new();

        let mut max_num_iters = vec_all_reads.len() / iters;
        if iteration == iters - 1 {
            max_num_iters = usize::MAX;
        }

        for (i, _read) in vec_all_reads.iter().enumerate() {
            let index = i as i32;
            if used_vertices.contains(&index) {
                continue;
            }

            let edges = &adj_list_edges[i];
            let mut dist_between_clusters = Vec::new();
            for _i in 0..ploidy {
                dist_between_clusters.push(-1.0);
            }

            for edge in edges.iter() {
                for (j, cluster) in clusters.iter().enumerate() {
                    if cluster.contains(&edge.1) {
                        let dist = edge.0;
                        if dist > dist_between_clusters[j] {
                            dist_between_clusters[j] = dist;
                        }
                    }
                }
            }

            let mut sorted_distances: Vec<_> = dist_between_clusters.iter().collect();
            sorted_distances.sort_by(|a, b| a.partial_cmp(&b).unwrap());
            //            dbg!(&sorted_distances);
            //            let mut sort_val = sorted_distances[0] - (sorted_distances[1] + sorted_distances[2] + sorted_distances[3]);
            let mut sort_val = sorted_distances[0] - sorted_distances[1];
            if *sorted_distances[0] == -1.0 {
                sort_val = f64::MAX;
            }
            sorted_vec_overlap_reads.push((i, sort_val));
        }

        //Obtained sorted vertices
        sorted_vec_overlap_reads.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        //        if sorted_vec_overlap_reads.len() > 0{
        //            dbg!(&sorted_vec_overlap_reads[0], vec_all_reads[sorted_vec_overlap_reads[0].0 as usize].seq_dict.len(),iteration);
        //        }

        //Now greedily add vertices to the partition where the maximum distance to the cluster is
        //minimized.
        for (j, (vertex, _dist)) in sorted_vec_overlap_reads.iter().enumerate() {
            if j > max_num_iters {
                break;
            }
            let edges = &adj_list_edges[*vertex];
            let mut max_dist = Vec::new();
            for _i in 0..ploidy {
                max_dist.push(-1.0);
            }
            for edge in edges.iter() {
                for (j, cluster) in clusters.iter().enumerate() {
                    if cluster.contains(&edge.1) {
                        let dist = edge.0;
                        if max_dist[j] < dist {
                            max_dist[j] = dist;
                        }
                    }
                }
            }
            //            println!("{:?},{:?} max_dist, overlap",max_dist,overlap);

            //Find index of minimum distance cluster
            let mut min_index = 0;
            let mut min_score = f64::MAX;
            for i in 0..ploidy {
                if max_dist[i] < min_score {
                    min_index = i;
                    min_score = max_dist[i];
                }
            }

            let vertex_i32 = *vertex as i32;
            used_vertices.insert(vertex_i32);
            clusters[min_index].insert(vertex_i32);
        }
    }
}

fn _populate_clusters1(
    clusters: &mut Vec<FxHashSet<i32>>,
    used_vertices: &mut FxHashSet<i32>,
    vec_all_reads: &Vec<&&Frag>,
    adj_list_edges: &Vec<Vec<(f64, i32)>>,
    ploidy: usize,
) {
    //A read must overlap with the intial clusters at least this much in order for it to be processed.
    //Once all good reads are processed, we sort the vertices again based on overlap and update the
    //read's overlap with the new clusters.
    let min_overlap = 50;
    let mut prev_used = 0;
    //If some reads just don't overlap super well, we need to relax the condition on the
    //min_overlap.
    let mut relax = false;
    loop {
        //First sort vertices by the minimum overlap for a read with all vertices in the clusters.
        //We do a loop here because once
        if prev_used == used_vertices.len() {
            relax = true;
        }

        prev_used = used_vertices.len();

        //Sort reads(vertices) by the minimum of the maximum overlaps between clusters
        let mut sorted_vec_overlap_reads = Vec::new();

        if used_vertices.len() == vec_all_reads.len() {
            break;
        }

        for (i, read) in vec_all_reads.iter().enumerate() {
            let index = i as i32;
            if used_vertices.contains(&index) {
                continue;
            }

            let edges = &adj_list_edges[i];
            let mut read_overlaps_between_clusters = Vec::new();
            for _i in 0..ploidy {
                read_overlaps_between_clusters.push(0);
            }

            for edge in edges.iter() {
                for (j, cluster) in clusters.iter().enumerate() {
                    if cluster.contains(&edge.1) {
                        let edge_index = edge.1 as usize;
                        let read2 = vec_all_reads[edge_index];
                        let t: Vec<_> = read.positions.intersection(&read2.positions).collect();
                        let overlap_len = t.len();
                        if overlap_len > read_overlaps_between_clusters[j] {
                            read_overlaps_between_clusters[j] = overlap_len;
                        }
                    }
                }
            }
            sorted_vec_overlap_reads.push((
                i,
                *read_overlaps_between_clusters.iter().min().unwrap() as usize,
            ));
        }

        //Obtained sorted vertices
        sorted_vec_overlap_reads.sort_by(|a, b| b.1.cmp(&a.1));

        //Now greedily add vertices to the partition where the maximum distance to the cluster is
        //minimized.
        for (vertex, overlap) in sorted_vec_overlap_reads.iter() {
            if !relax && *overlap < min_overlap {
                break;
            }
            let edges = &adj_list_edges[*vertex];
            let mut max_dist = Vec::new();
            for _i in 0..ploidy {
                max_dist.push(-1.0);
            }
            for edge in edges.iter() {
                for (j, cluster) in clusters.iter().enumerate() {
                    if cluster.contains(&edge.1) {
                        let dist = edge.0;
                        if max_dist[j] < dist {
                            max_dist[j] = dist;
                        }
                    }
                }
            }
            //            println!("{:?},{:?} max_dist, overlap",max_dist,overlap);

            //Find index of minimum distance cluster
            let mut min_index = 0;
            let mut min_score = f64::MAX;
            for i in 0..ploidy {
                if max_dist[i] < min_score {
                    min_index = i;
                    min_score = max_dist[i];
                }
            }

            let vertex_i32 = *vertex as i32;
            used_vertices.insert(vertex_i32);
            clusters[min_index].insert(vertex_i32);
        }
    }
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
        let prev_hap_block = utils_frags::hap_block_from_partition(&partition);
        return (0.0, partition, prev_hap_block);
    }

    let mut prev_hap_block = utils_frags::hap_block_from_partition(&partition);
    let mut set_of_positions = FxHashSet::default();

    for block in prev_hap_block.blocks.iter() {
        for pos in block.keys() {
            set_of_positions.insert(*pos);
        }
    }

    let binom_vec = get_mec_stats_epsilon(&partition, &prev_hap_block, epsilon);
    let mut prev_score = binom_vec.iter().map(|x| x.1).sum();
    prev_score *= -1.;

    let mut best_part = partition;

    //Iterate until an iteration yields a lower UPEM score -- return partition corresponding
    //to the best UPEM score.
    for i in 0..max_iters {
        let new_part = opt_iterate(&best_part, &prev_hap_block, epsilon);
        let new_block = utils_frags::hap_block_from_partition(&new_part);
        let new_binom_vec = get_mec_stats_epsilon(&new_part, &new_block,epsilon);
        let new_score = new_binom_vec.iter().map(|x| x.1).sum::<f64>() * -1.;
        if new_score > prev_score{
            log::trace!("Iter {} successful, new {} prev {}", i, new_score, prev_score);
        }

        if new_score > prev_score {
            prev_score = new_score;
            best_part = new_part;
            prev_hap_block = new_block;
        } else {
            log::trace!("Iter {} unsuccessful, new {} prev {}", i, new_score, prev_score);
            return (prev_score, best_part, prev_hap_block);
        }
    }

    return (prev_score, best_part, prev_hap_block);
}

//Get the chiq-square log p value from a vector of frequencies.
fn _chi_square_p(freqs: &Vec<usize>) -> f64 {
    let dof = (freqs.len() - 1) as f64;

    let mean: usize = freqs.iter().sum();
    let mean = mean as f64;
    let mean = mean / (freqs.len() as f64);
    let mut chi_stat = 0.0;
    for freq in freqs.iter() {
        chi_stat += ((*freq as f64) - mean).powf(2.0);
    }
    chi_stat /= mean;
    //We have to handle the case where all frequencies are the same separately or rust crashes.
    if chi_stat <= 0.00 {
        return 0.000;
    }
    let rv = ChiSquared::new(dof).unwrap();
    let rv_res = 1.0 - rv.cdf(chi_stat);
    return rv_res.ln();
}

//Can also use a normal approximation. This formula is taken from wikipedia.
pub fn log_erfc(x: f64) -> f64 {
    let p = 0.47047;
    let a1 = 0.3480242;
    let a2 = -0.0958798;
    let a3 = 0.7478556;
    let t = 1.0 / (1.0 + p * x);

    let polynomial_term = a1 * t + a2 * t * t + a3 * t * t * t;
    if polynomial_term < 0.0 {
        return 0.0;
    }
    let log_erfc = -(x * x) + polynomial_term.ln();

    return log_erfc;
}

pub fn norm_approx(n: usize, k: usize, p: f64, _div_factor: f64) -> f64 {
    let div_factor = 1.0;
    let samp_size = (n as f64) / div_factor;
    let mu = samp_size * p;
    let sigma = (mu * (1.0 - p)).sqrt();
    let z = (k as f64 + 0.5 - mu) / (sigma);
    return log_erfc(z);
}

//Get a vector of read frequencies and error rates from a partition and its corresponding
//haplotype block.
pub fn get_partition_stats(
    partition: &Vec<FxHashSet<&Frag>>,
    hap_block: &HapBlock,
) -> (Vec<(usize, usize)>, Vec<usize>) {
    let debug_chunks = false;
    let mut chunk_read_errors = vec![];
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

            //Collecting info about statistics here... TODO remove
            if debug_chunks{
                utils_frags::chunk_vec_update(frag, haplo, &mut chunk_read_errors);
            }

        }
        binom_vec.push((bases, errors));
        freq_vec.push(partition[i].len());
    }
    if debug_chunks{
        dbg!(chunk_read_errors);
    }

    (binom_vec, freq_vec)
}

//Include a pental for single coverage alleles. 
pub fn get_mec_stats_epsilon(
    _partition: &Vec<FxHashSet<&Frag>>,
    hap_block: &HapBlock,
    epsilon : f64
) -> Vec<(f64, f64)> {
    let mut binom_vec = vec![];
    for hap in hap_block.blocks.iter(){
        let mut errors = 0.;
        let mut bases = 0.;
        for seq_dict in hap.values(){
            let mut allele_counts : Vec<&usize>= seq_dict.values().collect();
            allele_counts.sort();
            let cons_bases = **allele_counts.last().unwrap();
            bases += cons_bases as f64;
            for i in 0..allele_counts.len()-1{
                errors += *allele_counts[i] as f64;
            }
            //TODO Test if this is worthwhile. 
            if cons_bases == 1{
                errors += epsilon;
            }
        }
        binom_vec.push((bases,errors));
    }
    binom_vec
}

//Get a vector of read frequencies and error rates from a partition and its corresponding
//haplotype block.
pub fn get_partition_stats_ref_wild(
    partition: &Vec<FxHashSet<&Frag>>,
    hap_block: &HapBlock,
) -> (Vec<((usize, usize),(usize,usize))>, Vec<usize>) {
    let mut binom_vec = Vec::new();
    let mut freq_vec = Vec::new();
    let ploidy = partition.len();
    for i in 0..ploidy {
        let haplo = &hap_block.blocks[i];
        let mut bases_ref = 0;
        let mut bases_wild= 0;
        let mut errors_ref = 0;
        let mut errors_wild = 0;
        for frag in partition[i].iter() {
            let ((same_ref, diff_ref),(same_wild,diff_wild)) = utils_frags::distance_read_haplo_ref_wild(frag, haplo);
            errors_ref += diff_ref;
            errors_wild+= diff_wild;
            bases_ref += same_ref;
            bases_wild+= same_wild;
        }
        binom_vec.push(((bases_ref, errors_ref),(bases_wild,errors_wild)));
        freq_vec.push(partition[i].len());
    }

    (binom_vec, freq_vec)
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

//Return upem score
fn _get_upem_score(
    binom_vec: &Vec<(usize, usize)>,
    freq_vec: &Vec<usize>,
    p: f64,
    div_factor: f64,
) -> f64 {
    let mut score = 0.0;
    for stat in binom_vec.iter() {
        let bincdf = utils_frags::stable_binom_cdf_p_rev(stat.0 + stat.1, stat.1, p, div_factor);
        score += bincdf;
    }
    score += _chi_square_p(freq_vec);
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
            let (_bases_good_read, errors_read) = utils_frags::distance_read_haplo_epsilon_empty(read, haplo_i, epsilon);
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
    if number_of_moves == 0 && best_moves.len() > 0{
        number_of_moves = best_moves.len()/3 + 1;
    }
    log::trace!("Number of moves {}", number_of_moves);

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

pub fn estimate_epsilon(
    num_iters: usize,
    num_tries: usize,
    ploidy: usize,
    all_frags: &Vec<Frag>,
    block_len: usize,
    initial_epsilon: f64,
) -> f64 {
    let mut rng = Pcg64::seed_from_u64(1);
    let mut random_vec = Vec::new();

    let mut epsilons = Vec::new();

    for _ in 0..num_tries {
        random_vec.push(rng.gen_range(0..num_iters));
    }

    let mut part_stats = vec![Vec::new();ploidy];
    for i in random_vec.into_iter() {
        let part = generate_hap_block(
            i * block_len,
            (i + 1) * block_len,
            ploidy,
            all_frags,
            initial_epsilon,
        );
        let mut part_lens: Vec<usize> = part.iter().map(|x| x.len()).collect();
        part_lens.sort();
        for j in 0..ploidy{
            part_stats[j].push(part_lens[j]);
        }
        let block = utils_frags::hap_block_from_partition(&part);
        let (binom_vec, _freq_vec) = get_partition_stats(&part, &block);
        for (good, bad) in binom_vec {
            if good + bad == 0 {
                break;
            }
            let epsilon = (bad as f64) / ((good + bad) as f64);
            epsilons.push(epsilon);
        }
    }

    let percentile_index = epsilons.len() / 10;
    let mut median_part_stats = vec![];
    for j in 0..ploidy{
        part_stats[j].sort();
        median_part_stats.push(part_stats[j][num_tries/2]);
    }

    log::debug!("Average count partitions from local clustering : {:?}",median_part_stats);
    epsilons.sort_by(|a, b| a.partial_cmp(&b).unwrap());
    epsilons[percentile_index]
}

pub fn estimate_ploidy_flopp(
    num_iters: usize,
    num_tries: usize,
    all_frags: &Vec<Frag>,
    initial_epsilon: f64,
) -> usize {
    let mut rng = Pcg64::seed_from_u64(1);
    let mut random_vec = Vec::new();

    for _ in 0..num_tries {
        random_vec.push(rng.gen_range(0..num_iters));
    }
    //println!("{:?}",random_vec);

    let block_len = 1;
    let ploidy_start = 2;
    let ploidy_end = 6;
    let error_rate = initial_epsilon;
    let num_ploidies = ploidy_end - ploidy_start;
    let mut mec_vector = vec![0;num_ploidies];
    let mut expected_errors_ref = vec![];
    for ploidy in ploidy_start..ploidy_end{
        let mut num_alleles = 0;
        for i in random_vec.iter(){
            let part = generate_hap_block(
                i * block_len,
                (i + 1) * block_len,
                ploidy,
                all_frags,
                initial_epsilon,
            );
            let mut part_lens: Vec<usize> = part.iter().map(|x| x.len()).collect();
            part_lens.sort();
            let block = utils_frags::hap_block_from_partition(&part);
            let (binom_vec, _freq_vec) = get_partition_stats(&part, &block);
            for (good,bad) in binom_vec {
                mec_vector[ploidy-2] += bad;
                num_alleles +=good;
                num_alleles +=bad;
            }
        }
        expected_errors_ref.push(num_alleles as f64 * error_rate);
    }

    log::debug!("{:?}, {:?}",mec_vector, &expected_errors_ref);
    for i in 0..num_ploidies-1{
        let mec_threshold = 1.0 / (1.0 - error_rate) / (1.0 + 1.0/(i+3) as f64 );
        if mec_vector[i] < expected_errors_ref[i] as usize{
            println!("MEC error threshold, returning ploidy.");
            return i + 2;
        }
        log::debug!("Expected MEC ratio {}, observed MEC ratio {}",mec_threshold, mec_vector[i+1] as f64 / mec_vector[i] as f64);
        if (mec_vector[i+1] as f64 / mec_vector[i] as f64) < mec_threshold{
        }
        else{
            println!("MEC decrease thereshold, returning ploidy.");
            return i + 2;
        }
    }

    
    return ploidy_end-1;
}
