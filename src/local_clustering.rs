use crate::types_structs::{Frag, HapBlock};
use crate::vcf_polishing;
use statrs::distribution::{Binomial, ChiSquared, Univariate};
extern crate time;
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use std::time::Instant;

//Return the set of reads for which every read covers at least one position in the interval
//(inclusive)
pub fn find_reads_in_interval<'a>(
    start: usize,
    end: usize,
    //position_to_reads : &FxHashMap<usize,FxHashSet<&Frag>>,
    all_frags: &'a Vec<Frag>,
) -> FxHashSet<&'a Frag> {
    let mut final_set = FxHashSet::default();
    //This is slower than just iterating thorugh the entire fragment list. We can speed this up by
    //indexing the fragments as well.

    //    for i in (start..end+1).step_by(10){
    //        let empty_set = &FxHashSet::default();
    //        let set1 = position_to_reads.get(&(i)).unwrap_or(empty_set);
    //        for elt in set1.iter() {
    //            final_set.insert(*elt);
    //        }
    //    }

    for frag in all_frags.iter() {
        if frag.last_position < start {
            continue;
        }
        if frag.first_position > end {
            break;
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
) -> Vec<FxHashSet<&'a Frag>> {
    let start_t = Instant::now();
    let all_reads = find_reads_in_interval(start, end, all_frags);
    println!(
        "Time taken get overlaps {:?}",
        start_t.elapsed().as_millis()
    );
    let partition = cluster_reads(&all_reads, ploidy);
    partition
}

//Compute distance between two reads from the precomputed hash map
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

//Local clustering method for a set of reads -> partition
pub fn cluster_reads<'a>(
    all_reads: &FxHashSet<&'a Frag>,
    ploidy: usize,
) -> Vec<FxHashSet<&'a Frag>> {
    //Get the largest distance edge
    let mut vec_all_edges = Vec::new();
    let vec_all_reads: Vec<_> = all_reads.iter().collect();
    let mut adj_list_edges = Vec::new();
    for _i in 0..vec_all_reads.len() {
        adj_list_edges.push(Vec::new());
    }
    println!("Computing edges for local cluster");
    let start_t = Instant::now();
    for (i, r1) in vec_all_reads.iter().enumerate() {
        for j in i + 1..vec_all_reads.len() {
            //
            let r2 = &vec_all_reads[j];
            //            let dist_try = dist_from_graph(*r1,*r2,&all_distances);
            //            let dist = match dist_try{
            //                Some(x) => x,
            //                None => continue,
            //            };
            //
            let mut dist;
            if !utils_frags::check_overlap(r1, r2) {
                continue;
            } else {
                dist = utils_frags::distance(r1, r2);
            }

            let i_type = i as i32;
            let j_type = j as i32;
            vec_all_edges.push(vec![dist, i_type, j_type]);
            let edge_list1 = &mut adj_list_edges[i];
            edge_list1.push(vec![dist, j_type]);
            let edge_list2 = &mut adj_list_edges[j];
            edge_list2.push(vec![dist, i_type]);
        }
    }

    println!(
        "Time edges local cluster {:?}",
        start_t.elapsed().as_millis()
    );
    println!("Finding max clique");
    let start_t = Instant::now();
    vec_all_edges.sort_by(|a, b| a[0].cmp(&b[0]));
    let best_edge = vec_all_edges.last().unwrap();
    //    println!("{:?}",vec_all_edges);

    let mut used_vertices = FxHashSet::default();
    used_vertices.insert(best_edge[1]);
    used_vertices.insert(best_edge[2]);

    //Greedily find a max clique once the first two vertices are found
    for _i in 0..ploidy - 2 {
        let mut min_dist_map = FxHashMap::default();
        for edge in vec_all_edges.iter() {
            if used_vertices.contains(&edge[1]) && !used_vertices.contains(&edge[2]) {
                if min_dist_map.contains_key(&edge[2]) {
                    if *min_dist_map.get(&edge[2]).unwrap() < edge[0] {
                        continue;
                    }
                }
                min_dist_map.insert(edge[2], edge[0]);
            } else if used_vertices.contains(&edge[2]) && !used_vertices.contains(&edge[1]) {
                if min_dist_map.contains_key(&edge[1]) {
                    if *min_dist_map.get(&edge[1]).unwrap() < edge[0] {
                        continue;
                    }
                }
                min_dist_map.insert(edge[1], edge[0]);
            }
        }
        let mut sorted_dict_to_vec: Vec<_> = min_dist_map.into_iter().collect();
        sorted_dict_to_vec.sort_by(|a, b| a.1.cmp(&b.1));
        //        println!("{:?}",sorted_dict_to_vec);
        //        println!("{:?}",vec_all_edges);
        let best_vertex = sorted_dict_to_vec.last().unwrap();
        used_vertices.insert(best_vertex.0);
    }

    let mut clusters = Vec::new();
    for vertex in used_vertices.iter() {
        let mut cluster = FxHashSet::default();
        cluster.insert(*vertex);
        clusters.push(cluster);
    }

    println!("Greedy partitioning...");
    //Once seed vertices for each cluster is found, greedily add edges to each cluster based on
    //minimizing the max dist over clusters

    let mut relax = false;
    let min_overlap = 2;
    let mut prev_used = 0;
    loop {
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
                    if cluster.contains(&edge[1]) {
                        let edge_index = edge[1] as usize;
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
        //Obtain sorted vertices
        sorted_vec_overlap_reads.sort_by(|a, b| b.1.cmp(&a.1));
        //        println!("{:?} sorted_vec_overlap",sorted_vec_overlap_reads);

        //Now greedily add vertices to the partition where the maximum distance to the cluster is
        //minimized.
        for (vertex, overlap) in sorted_vec_overlap_reads.iter() {
            if !relax && *overlap < min_overlap {
                break;
            }
            let edges = &adj_list_edges[*vertex];
            let mut max_dist = Vec::new();
            for _i in 0..ploidy {
                max_dist.push(i32::MIN);
            }
            for edge in edges.iter() {
                for (j, cluster) in clusters.iter().enumerate() {
                    if cluster.contains(&edge[1]) {
                        let dist = edge[0];
                        if max_dist[j] < dist {
                            max_dist[j] = dist;
                        }
                    }
                }
            }
            //            println!("{:?},{:?} max_dist, overlap",max_dist,overlap);

            //Find index of minimum distance cluster
            let mut min_index = 0;
            let mut min_score = i32::MAX;
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
    println!(
        "Time greedly local clustering {:?}",
        Instant::now() - start_t
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

pub fn optimize_clustering<'a>(
    partition: Vec<FxHashSet<&'a Frag>>,
    epsilon: f64,
    genotype_dict: &FxHashMap<usize, FxHashMap<usize, usize>>,
    polish: bool,
    max_iters: usize,
) -> (f64, Vec<FxHashSet<&'a Frag>>,HapBlock) {

    let mut prev_hap_block = utils_frags::hap_block_from_partition(&partition);
    let mut set_of_positions = FxHashSet::default();

    for block in prev_hap_block.blocks.iter(){
        for pos in block.keys(){
            set_of_positions.insert(*pos);
        }
    }

    let position_vec: Vec<usize> = set_of_positions.into_iter().collect();

    if polish {
        prev_hap_block =
            vcf_polishing::polish_using_vcf(genotype_dict, &prev_hap_block, &position_vec);
    }

    let (binom_vec, freq_vec) = get_partition_stats(&partition, &prev_hap_block);
    //dbg!(&binom_vec,&freq_vec);
    let mut prev_score = get_upem_score(&binom_vec, &freq_vec, epsilon);
    let mut best_part = partition;


    dbg!(prev_score,freq_vec,binom_vec);

    for _i in 0..max_iters {
        let new_part = opt_iterate(&best_part, &prev_hap_block,epsilon);
        let mut new_block = utils_frags::hap_block_from_partition(&new_part);
        if polish{
            new_block = vcf_polishing::polish_using_vcf(genotype_dict,&new_block,&position_vec);
        }
        let (new_binom_vec,new_freq_vec) = get_partition_stats(&new_part,&new_block);
        let new_score = get_upem_score(&new_binom_vec,&new_freq_vec,epsilon);
        dbg!(new_score,new_freq_vec,new_binom_vec);
        if new_score > prev_score {
            prev_score = new_score;
            best_part = new_part;
            prev_hap_block = new_block;
        } else {
            return (prev_score, best_part,prev_hap_block);
        }
    }

    return (prev_score, best_part,prev_hap_block);
}

//Get the chiq-square log p value from a vector of frequencies.
fn chi_square_p(freqs: &Vec<usize>) -> f64 {
    let dof = (freqs.len() - 1) as f64;

    let mean: usize = freqs.iter().sum();
    let mean = mean as f64;
    let mean = mean / (freqs.len() as f64);
    let mut chi_stat = 0.0;
    for freq in freqs.iter() {
        chi_stat += ((*freq as f64) - mean).powf(2.0);
    }
    chi_stat /= mean;
    let rv = ChiSquared::new(dof).unwrap();
    let rv_res = 1.0 - rv.cdf(chi_stat);
    return rv_res.ln();
}

//Get the log p-value for a 1-sided binomial test.
fn binom_cdf_p(n: usize, k: usize, p: f64) -> f64 {
    let rv = Binomial::new(p, n as u64).unwrap();
    let rv_result = rv.cdf(k as f64);
    //dbg!((1.0-rv_result).ln());
    return (1.0 - rv_result).ln();
}

fn stable_binom_cdf_p_rev(n : usize, k : usize, p : f64) -> f64{
    let n64 = n as f64;
    let k64 = k as f64;
    let a = (n64-k64)/n64;
    let p = 1.0-p;
    let rel_ent = a * (a/p).ln() + (1.0 - a) * ((1.0-a)/(1.0-p)).ln();
    return -1.0*n64/25.0*rel_ent;
}

//Get a vector of read frequencies and error rates from a partition and its corresponding
//haplotype block.
fn get_partition_stats(
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

//Return upem score
fn get_upem_score(binom_vec: &Vec<(usize, usize)>, freq_vec: &Vec<usize>, p: f64) -> f64 {
    let mut score = 0.0;
    for stat in binom_vec.iter() {
        let bincdf = stable_binom_cdf_p_rev(stat.0 + stat.1, stat.1, p);
        score += bincdf;
    }
    //dbg!(chi_square_p(freq_vec));
    score += chi_square_p(freq_vec);
    score
}

fn opt_iterate<'a>(
    partition: &Vec<FxHashSet<&'a Frag>>,
    hap_block: &HapBlock,
    epsilon: f64,
) -> Vec<FxHashSet<&'a Frag>> {
    let ploidy = partition.len();
    let (binom_vec, freq_vec) = get_partition_stats(partition, hap_block);
    let mut freq_vec = freq_vec;
    let mut binom_p_vec = Vec::new();
    let chi_square_val = chi_square_p(&freq_vec);

    for bases_errors in binom_vec.iter() {
        let bases = bases_errors.0;
        let errors = bases_errors.1;
        let binom_logp_val = stable_binom_cdf_p_rev(bases + errors, errors, epsilon);
        binom_p_vec.push(binom_logp_val);
    }

    let mut best_moves = Vec::new();

    for i in 0..ploidy {
        if partition[i].len() <= 1{
            continue;
        }
        for read in partition[i].iter() {
            let haplo_i = &hap_block.blocks[i];
            let (bases_good_read, errors_read) = utils_frags::distance_read_haplo(read, haplo_i);
            let bases_good_after = binom_vec[i].0 - bases_good_read;
            let errors_after = binom_vec[i].1 - errors_read;
            let new_binom_val_i =
                stable_binom_cdf_p_rev(bases_good_after + errors_after, errors_after, epsilon);
            for j in 0..ploidy {
                if j == i {
                    continue;
                }

                //Test out new move
                let haplo_j = &hap_block.blocks[j];
                let (read_bases_good_movej, read_errors_movej) =
                    utils_frags::distance_read_haplo(read, haplo_j);

                let bases_good_after_movej = binom_vec[j].0 + read_bases_good_movej;
                let errors_after_movej = binom_vec[j].1 + read_errors_movej;
                let new_binom_val_j = stable_binom_cdf_p_rev(
                    bases_good_after_movej + errors_after_movej,
                    errors_after_movej,
                    epsilon,
                );

                freq_vec[j] += 1;
                freq_vec[i] -= 1;
                let new_chi_square_val = chi_square_p(&freq_vec);

                let new_score = new_binom_val_i + new_binom_val_j + new_chi_square_val;
                let old_score = binom_p_vec[i] + binom_p_vec[j] + chi_square_val;
                if new_score - old_score > 0.0 {
                    best_moves.push((new_score-old_score,(i,read,j)));
//                    dbg!(new_score,old_score,read_bases_good_movej,read_errors_movej);
                }

                freq_vec[i] += 1;
                freq_vec[j] -= 1;

            }
        }
    }

    let mut moved_reads = FxHashSet::default();
    let mut new_part = partition.clone();
    best_moves.sort_by(|a,b| b.0.partial_cmp(&a.0).unwrap());
    let num_reads : usize = freq_vec.iter().sum();
    let mut number_of_moves = num_reads/10;
    if best_moves.len()/10 < number_of_moves/5{
        number_of_moves = best_moves.len()/5;
    }
    dbg!(number_of_moves);

    for (mv_num,mv) in best_moves.iter().enumerate(){
        let (i,read,j) = mv.1;
        if moved_reads.contains(read){
            continue;
        }
        if freq_vec[i] == 1{
            continue;
        }
        new_part[j].insert(read);
        new_part[i].remove(read);
        freq_vec[j] += 1;
        freq_vec[i] -= 1;
        moved_reads.insert(read);
        if mv_num > number_of_moves{
            break;
        }
    }

    new_part
}
