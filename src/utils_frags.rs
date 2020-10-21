use crate::types_structs::Frag;
use crate::types_structs::HapBlock;
use fxhash::{FxHashMap, FxHashSet};

// Get the number # of different bases between the
// two fragments
pub fn distance(r1: &Frag, r2: &Frag) -> (i32,i32) {
    let mut diff = 0;
    let mut same = 0;

    for pos in r1.positions.intersection(&r2.positions) {
        if r1.seq_dict.get(pos) == r2.seq_dict.get(pos) {
            same += 1;
        } else {
            diff += 1;
        }
    }

    (same,diff)
}

pub fn distance_read_haplo(
    r1: &Frag,
    hap: &FxHashMap<usize, FxHashMap<usize, usize>>,
) -> (usize, usize) {
    let mut diff = 0;
    let mut same = 0;
    for pos in r1.positions.iter() {
        if !hap.contains_key(pos) {
            continue;
        }

        let frag_var = r1.seq_dict.get(pos).unwrap();
        let consensus_var = hap
            .get(pos)
            .unwrap()
            .iter()
            .max_by_key(|entry| entry.1)
            .unwrap()
            .0;
        if *frag_var == *consensus_var {
            same += 1;
        } else {
            diff += 1;
        }
    }

    (same, diff)
}

pub fn distance_read_haplo_range(
    r1: &Frag,
    hap: &FxHashMap<usize, FxHashMap<usize, usize>>,
    start : usize,
    end : usize
) -> (usize, usize) {
    let mut diff = 0;
    let mut same = 0;
    for pos in r1.positions.iter() {
        if !hap.contains_key(pos) || *pos > end || *pos < start {
            continue;
        }

        let frag_var = r1.seq_dict.get(pos).unwrap();
        let consensus_var = hap
            .get(pos)
            .unwrap()
            .iter()
            .max_by_key(|entry| entry.1)
            .unwrap()
            .0;
        if *frag_var == *consensus_var {
            same += 1;
        } else {
            diff += 1;
        }
    }

    (same, diff)
}

//Index each position by the set of fragments which overlap that position
pub fn get_all_overlaps(frags: &Vec<Frag>) -> FxHashMap<usize, FxHashSet<&Frag>> {
    let mut overlaps = FxHashMap::default();
    for frag in frags.iter() {
        for pos in frag.positions.iter() {
            let pos_set = overlaps.entry(*pos).or_insert(FxHashSet::default());
            pos_set.insert(frag);
        }
    }

    overlaps
}

//Find all distances between fragments.
//Assumes sorted fragments by first position. 
pub fn get_all_distances(frags: &Vec<Frag>) -> FxHashMap<&Frag, FxHashMap<&Frag, i32>> {
    let mut pairwise_distances = FxHashMap::default();

    for (i, frag1) in frags.iter().enumerate() {
        let mut from_frag_distance = FxHashMap::default();
        for j in i..frags.len() {
            let frag2 = &frags[j];

            if frag1.last_position < frag2.first_position {
                break;
            }

            from_frag_distance.insert(frag2, distance(frag1, frag2).1);
        }
        pairwise_distances.insert(frag1, from_frag_distance);
    }

    pairwise_distances
}

pub fn check_overlap(r1: &Frag, r2: &Frag) -> bool {
    if r1.last_position < r2.first_position {
        return false;
    }
    if r2.last_position < r1.first_position {
        return false;
    }
    let t: Vec<_> = r1.positions.intersection(&r2.positions).collect();
    if t.len() == 0 {
        return false;
    } else {
        return true;
    }
}

pub fn hap_block_from_partition(part: &Vec<FxHashSet<&Frag>>) -> HapBlock {
    let mut block_vec = Vec::new();
    for reads in part.iter() {
        let mut hap_map = FxHashMap::default();
        for frag in reads.iter() {
            for pos in frag.positions.iter() {
                let var_at_pos = frag.seq_dict.get(pos).unwrap();
                let sites = hap_map.entry(*pos).or_insert(FxHashMap::default());
                let site_counter = sites.entry(*var_at_pos).or_insert(0);
                *site_counter += 1;
            }
        }
        block_vec.push(hap_map);
    }
    HapBlock { blocks: block_vec }
}

pub fn get_avg_length(all_frags : &Vec<Frag>, quantile : f64) -> usize{
   let mut length_vec = Vec::new(); 
   for frag in all_frags.iter(){
       length_vec.push(frag.last_position - frag.first_position);
   }
   length_vec.sort();
   return length_vec[(length_vec.len() as f64 * quantile) as usize];
}

pub fn get_length_gn(all_frags : &Vec<Frag>) -> usize{
    let mut last_pos = 0;
    for frag in all_frags.iter(){
        if frag.last_position > last_pos{
            last_pos = frag.last_position;
        }
    }
    last_pos
}
