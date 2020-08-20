use std::time::Instant;
use haplotype_phaser::file_reader;
use haplotype_phaser::local_clustering;
use haplotype_phaser::utils_frags;
use fxhash::{FxHashSet,FxHashMap};
use std::collections::{HashSet,BTreeMap};

#[test]
fn frag_reader_test() {
    let frags = file_reader::get_frags_container("/home/jshaw/haplotype_phaser/tests/test_file.txt");
    assert_eq!(frags.len(),3);
    assert_eq!(frags[0].id,"t1");
    assert_eq!(frags[0].counter_id,0);
    assert_eq!(frags[0].positions.len(),2);
    assert_eq!(frags[2].id,"t3");
    assert_eq!(frags[2].positions.len(),4);
}

#[test]
fn utils_frags_test(){
    let frags = file_reader::get_frags_container("/home/jshaw/haplotype_phaser/tests/test_file.txt");
    assert_eq!(utils_frags::distance(&frags[0],&frags[1]),1);
    assert_eq!(utils_frags::distance(&frags[1],&frags[2]),0);

    let all_distances = utils_frags::get_all_distances(&frags);
    assert_eq!(*(all_distances.get(&frags[0]).unwrap().get(&frags[1]).unwrap()),1);
    assert_eq!(*(all_distances.get(&frags[1]).unwrap().get(&frags[2]).unwrap()),0);
    assert_eq!(*(all_distances.get(&frags[0]).unwrap().get(&frags[2]).unwrap()),1);
}

#[test]
fn frags_test(){
    let frags = file_reader::get_frags_container("/home/jshaw/haplotype_phaser/tests/test_file.txt");
    let mut hashset = FxHashSet::default();
    hashset.insert(&frags[0]);
    hashset.insert(&frags[0]);
    assert_eq!(hashset.len(),1);
    hashset.insert(&frags[1]);
    assert_eq!(hashset.len(),2);
    hashset.insert(&frags[2]);
    assert_eq!(hashset.len(),3);
}

#[test]
fn local_cluster_test(){
    let frags = file_reader::get_frags_container("/home/jshaw/haplotype_phaser/tests/test_file.txt");
    let indexed_reads = utils_frags::get_all_overlaps(&frags);
    let interval_reads_all  = local_clustering::find_reads_in_interval(1,100,&frags);
    let interval_reads_5= local_clustering::find_reads_in_interval(5,6,&frags);

    assert_eq!(interval_reads_all.len(),3);
    assert_eq!(interval_reads_5.len(),1);

}

#[test]
fn speed_test(){
    let frags = file_reader::get_frags_container("/home/jshaw/haplotype_phaser/tests/test_file.txt");
    let mut hashset = FxHashSet::default();
    let mut hashmap = BTreeMap::new();

    let start_t = Instant::now();
    hashset.insert(&frags[0]);
    hashset.insert(&frags[1]);
    println!("Time taken to insert 1 set frag {:?}", Instant::now()-start_t);

    let start_t = Instant::now();
    hashmap.insert(0,1);
    hashmap.insert(2,1);
    println!("Time taken to insert 1 map frag {:?}", Instant::now()-start_t);

    let start_t = Instant::now();
    hashset.get(&frags[0]);
    println!("Time taken to get frag {:?}", Instant::now()-start_t);

    let start_t = Instant::now();
    hashset.contains(&frags[0]);
    println!("Time taken to check frag {:?}", Instant::now()-start_t);

    let mut hashset_int = FxHashSet::default();;

    let start_t = Instant::now();
    hashset_int.insert(0);
    hashset_int.insert(1);
    hashset_int.insert(2);
    println!("Time taken to insert 2int {:?}", Instant::now()-start_t);

    let start_t = Instant::now();
    hashset_int.get(&0);
    println!("Time taken to get int {:?}", Instant::now()-start_t);

    let start_t = Instant::now();
    utils_frags::distance(&frags[2],&frags[1]);
    println!("Time taken to get distance {:?}", Instant::now()-start_t);

}
