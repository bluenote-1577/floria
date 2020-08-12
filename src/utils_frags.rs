use crate::types_structs::Frag;
use crate::types_structs::HapBlock;
use fnv::{FnvHashMap,FnvHashSet};

// Get the number # of different bases between the
// two fragments
pub fn distance(r1 : &Frag, r2 : &Frag) -> i32{

    let mut diff = 0;
    let mut _same = 0;

    for pos in r1.positions.intersection(&r2.positions){
        if r1.seq_dict.get(pos) == r2.seq_dict.get(pos){
            _same += 1;
        }
        else{
            diff += 1;
        }
    }
    
    diff
}

//Index each position by the set of fragments which overlap that position
pub fn get_all_overlaps(frags : &Vec<Frag>) -> FnvHashMap<usize,FnvHashSet<&Frag>>{
    let mut overlaps = FnvHashMap::default();
    for frag in frags.iter(){
        for pos in frag.positions.iter(){
            let pos_set = overlaps.entry(*pos).or_insert(FnvHashSet::default());
            pos_set.insert(frag);
        }
    }

    overlaps
}

//Find all distances between fragments. 
//IMPORTANT!!! ASSUMES SORTED FRAGMENT .txt FILE BY STARTING POSITION!
pub fn get_all_distances(frags : &Vec<Frag>) -> FnvHashMap<&Frag,FnvHashMap<&Frag,i32>>{

    let mut pairwise_distances = FnvHashMap::default();

    for (i,frag1) in frags.iter().enumerate(){
        let mut from_frag_distance = FnvHashMap::default();
        for j in i..frags.len(){
            let frag2 = &frags[j];

            if frag1.last_position < frag2.first_position {
                break; 
            }

            from_frag_distance.insert(frag2,distance(frag1,frag2));
        }
        pairwise_distances.insert(frag1,from_frag_distance);
    }

    pairwise_distances

}

pub fn check_overlap(r1 : &Frag, r2 : &Frag) -> bool {
    if r1.last_position < r2.first_position{
        return false
    }
    if r2.last_position < r1.first_position{
        return false
    }
    let t : Vec<_> = r1.positions.intersection(&r2.positions).collect();
    if t.len()  == 0{
        return false
    }
    else {
        return true
    }
}

pub fn hap_block_from_partition(part : &Vec<FnvHashSet<&Frag>>) -> HapBlock{
    let mut block_vec = Vec::new();
    for reads in part.iter(){
        let mut hap_map = FnvHashMap::default();
        for frag in reads.iter(){
            for pos in frag.positions.iter(){
                let var_at_pos = frag.seq_dict.get(pos).unwrap();
                let sites = hap_map.entry(*pos).or_insert(FnvHashMap::default());
                let site_counter = sites.entry(*var_at_pos).or_insert(0);
                *site_counter += 1;
            }
        }
        block_vec.push(hap_map);
    }
    HapBlock{blocks : block_vec,}
}
