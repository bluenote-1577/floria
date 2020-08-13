use crate::types_structs::Frag;
extern crate time;
use time::PreciseTime;
use crate::utils_frags;
use fnv::{FnvHashMap,FnvHashSet};

//Return the set of reads for which every read covers at least one position in the interval
//(inclusive)
pub fn find_reads_in_interval<'a>(start : usize,
                              end : usize,
                              position_to_reads : &'a FnvHashMap<usize,FnvHashSet<&Frag>>) 
    -> FnvHashSet<&'a Frag>{

    let mut final_set =  FnvHashSet::default();
    for i in start..end+1{
        let empty_set = &FnvHashSet::default();
        let set1 = position_to_reads.get(&(i)).unwrap_or(empty_set);
        for elt in set1.iter() {
            final_set.insert(*elt);
        }
    }

    final_set
}

//Return a partition from a set of reads using our local clustering method. 
pub fn generate_hap_block<'a>(start : usize, 
                          end : usize, 
                          position_to_reads : &'a FnvHashMap<usize,FnvHashSet<&Frag>>,
                          all_distances : &FnvHashMap<&Frag,FnvHashMap<&Frag,i32>>,
                          ploidy : usize) -> Vec<FnvHashSet<&'a Frag>>{

    let start_t = PreciseTime::now();
    let all_reads = find_reads_in_interval(start,end,position_to_reads);
    println!("Time taken get overlaps {:?}",start_t.to(PreciseTime::now()));
    let partition = cluster_reads(&all_reads,&all_distances,ploidy);
    partition
    
}

//Compute distance between two reads from the precomputed hash map 
pub fn dist_from_graph(r1 : &Frag,
                       r2 : &Frag,
                       all_distances : &FnvHashMap<&Frag,FnvHashMap<&Frag,i32>>) -> Option<i32>{
    if !utils_frags::check_overlap(r1,r2){
        None
    }
    else{
        let map1 = all_distances.get(r1).unwrap();
        let dist;

        if map1.contains_key(r2){
            dist = map1.get(r2).unwrap();
        }

        else{
            dist = all_distances.get(r2).unwrap().get(r1).unwrap();
        }
        Some(*dist)
    }
}

//Local clustering method for a set of reads -> partition
pub fn cluster_reads<'a>(all_reads : &FnvHashSet<&'a Frag>, 
                      all_distances : &FnvHashMap<&Frag,FnvHashMap<&Frag,i32>>,
                      ploidy : usize) -> Vec<FnvHashSet<&'a Frag>>{

    //Get the largest distance edge
    let mut vec_all_edges = Vec::new();
    let mut adj_list_edges = FnvHashMap::default();
    let vec_all_reads : Vec<_> = all_reads.iter().collect();
    println!("Computing edges for local cluster");
    for (i,r1) in vec_all_reads.iter().enumerate(){
        for j in i+1..vec_all_reads.len(){

            let r2 = &vec_all_reads[j];
//            let dist_try = dist_from_graph(*r1,*r2,&all_distances);
//            let dist = match dist_try{
//                Some(x) => x,
//                None => continue,
//            };

            let mut dist;
            if !utils_frags::check_overlap(r1,r2){
                continue;
            }
            else{
                dist = utils_frags::distance(r1,r2);
            }

            let i_type = i as i32;
            let j_type = j as i32;
            vec_all_edges.push(vec![dist,i_type,j_type]);
            let edge_list1 = adj_list_edges.entry(i_type).or_insert(Vec::new());
            edge_list1.push(vec![dist,j_type]);
            let edge_list2 = adj_list_edges.entry(j_type).or_insert(Vec::new());
            edge_list2.push(vec![dist,i_type]);

        }
    }

    println!("Finding max clique");
    vec_all_edges.sort_by(|a,b| a[0].cmp(&b[0]));
    let best_edge = vec_all_edges.last().unwrap();
//    println!("{:?}",vec_all_edges);

    let mut used_vertices = FnvHashSet::default();
    used_vertices.insert(best_edge[1]);
    used_vertices.insert(best_edge[2]);

    //Greedily find a max clique once the first two vertices are found 
    for _i in 0..ploidy-2{
        let mut min_dist_map = FnvHashMap::default();
        for edge in vec_all_edges.iter(){
            if used_vertices.contains(&edge[1]) &&
                !used_vertices.contains(&edge[2]){
                if min_dist_map.contains_key(&edge[2]){
                    if *min_dist_map.get(&edge[2]).unwrap() <
                        edge[0]{
                        continue
                    }
                }
                min_dist_map.insert(edge[2],edge[0]);
            }
            else if used_vertices.contains(&edge[2])&&
                !used_vertices.contains(&edge[1]){

                if min_dist_map.contains_key(&edge[1]){
                    if *min_dist_map.get(&edge[1]).unwrap() < edge[0]{
                        continue
                    }
                }
                min_dist_map.insert(edge[1],edge[0]);
            }
        }
        let mut sorted_dict_to_vec : Vec<_> = min_dist_map.into_iter().collect();
        sorted_dict_to_vec.sort_by(|a,b| a.1.cmp(&b.1));
//        println!("{:?}",sorted_dict_to_vec);
//        println!("{:?}",vec_all_edges);
        let best_vertex = sorted_dict_to_vec.last().unwrap();
        used_vertices.insert(best_vertex.0);
    }

    let mut clusters = Vec::new();
    for vertex in used_vertices.iter(){
        let mut cluster = FnvHashSet::default();
        cluster.insert(*vertex);
        clusters.push(cluster);
    }


    println!("Greedy partitioning...");
    //Once seed vertices for each cluster is found, greedily add edges to each cluster based on
    //minimizing the max dist over clusters
    
    let mut relax = false;
    let min_overlap= 2;
    let mut prev_used = 0;
    loop{

        if prev_used == used_vertices.len(){
            relax = true;
        }

        prev_used = used_vertices.len();

        //Sort reads(vertices) by the minimum of the maximum overlaps between clusters
        let mut sorted_vec_overlap_reads = Vec::new();

        if used_vertices.len() == vec_all_reads.len(){
            break;
        }

        for (i,read) in vec_all_reads.iter().enumerate(){

            let index = i as i32;
            if used_vertices.contains(&index){
                continue
            }

            let edges = adj_list_edges.get(&index).unwrap();
            let mut read_overlaps_between_clusters = Vec::new();
            for _i in 0..ploidy{
                read_overlaps_between_clusters.push(0);
            }

            for edge in edges.iter(){
                for (j,cluster) in clusters.iter().enumerate(){
                    if cluster.contains(&edge[1]){
                        let edge_index = edge[1] as usize;
                        let read2 = vec_all_reads[edge_index];
                        let t : Vec<_> = read.positions.intersection(&read2.positions).collect();
                        let overlap_len = t.len();
                        if overlap_len > read_overlaps_between_clusters[j]{
                            read_overlaps_between_clusters[j] = overlap_len;
                        }
                    }
                }
            }
            sorted_vec_overlap_reads.push((i,*read_overlaps_between_clusters.iter().min().unwrap() as usize));
        }
        //Obtain sorted vertices 
        sorted_vec_overlap_reads.sort_by(|a,b| b.1.cmp(&a.1));
//        println!("{:?} sorted_vec_overlap",sorted_vec_overlap_reads);

        //Now greedily add vertices to the partition where the maximum distance to the cluster is
        //minimized. 
        for (vertex,overlap)  in sorted_vec_overlap_reads.iter(){
            if !relax && *overlap < min_overlap{
                break;
            }
            let v_index = *vertex as i32;
            let edges = adj_list_edges.get(&v_index).unwrap();
            let mut max_dist = Vec::new();
            for _i in 0..ploidy{
                max_dist.push(i32::MIN);
            }
            for edge in edges.iter(){
                for (j,cluster) in clusters.iter().enumerate(){
                    if cluster.contains(&edge[1]){
                        let dist = edge[0];
                        if max_dist[j] < dist{
                            max_dist[j] = dist;
                        }
                    }
                }
            }
//            println!("{:?},{:?} max_dist, overlap",max_dist,overlap);

            //Find index of minimum distance cluster
            let mut min_index = 0;
            let mut min_score = i32::MAX;
            for i in 0..ploidy{
                if max_dist[i] < min_score{
                    min_index = i;
                    min_score = max_dist[i];
                }
            }

            let vertex_i32 = *vertex as i32;
            used_vertices.insert(vertex_i32);
            clusters[min_index].insert(vertex_i32);
        }
    }

    //Turn the vertex indices into actual fragments -- could probably come up with a more elegant
    //solution using some sort of map function...
    let mut partition = Vec::new();
    for cluster in clusters.iter(){
        let mut frag_set = FnvHashSet::default();
        for vertex in cluster.iter(){
            let vertex_usize = *vertex as usize;
            frag_set.insert(*vec_all_reads[vertex_usize]);
        }
        partition.push(frag_set)
    }
    partition

}

