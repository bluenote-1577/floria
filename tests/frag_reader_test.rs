use std::time::Instant;
use flopp::file_reader;
use flopp::local_clustering;
use flopp::utils_frags;
use fxhash::{FxHashSet,FxHashMap};
use std::collections::{HashSet,BTreeMap};

#[test]
fn frag_reader_test() {
    let flopp_dir = "/home/jshaw/practical_prob_2020_paper/flopp/";
    let frags_map = file_reader::get_frags_container(flopp_dir.to_owned() + "/tests/test_file.txt");
    let frags = frags_map.get("frag_contig").unwrap();
    assert_eq!(frags.len(),3);
    assert_eq!(frags[0].id,"t1");
    assert_eq!(frags[0].counter_id,0);
    assert_eq!(frags[0].positions.len(),2);
    assert_eq!(frags[2].id,"t3");
    assert_eq!(frags[2].positions.len(),4);
}

#[test]
fn utils_frags_test(){
    let flopp_dir = "/home/jshaw/practical_prob_2020_paper/flopp/";
    let frags_map = file_reader::get_frags_container(flopp_dir.to_owned() + "/tests/test_file.txt");
    let frags = frags_map.get("frag_contig").unwrap();

    assert_eq!(utils_frags::distance(&frags[0],&frags[1]).1,1);
    assert_eq!(utils_frags::distance(&frags[1],&frags[2]).1,0);

    let all_distances = utils_frags::get_all_distances(frags);
    assert_eq!(*(all_distances.get(&frags[0]).unwrap().get(&frags[1]).unwrap()),1);
    assert_eq!(*(all_distances.get(&frags[1]).unwrap().get(&frags[2]).unwrap()),0);
    assert_eq!(*(all_distances.get(&frags[0]).unwrap().get(&frags[2]).unwrap()),1);
}

#[test]
fn frags_test(){
    let flopp_dir = "/home/jshaw/practical_prob_2020_paper/flopp/";
    let frags_map = file_reader::get_frags_container(flopp_dir.to_owned() + "/tests/test_file.txt");
    let frags = frags_map.get("frag_contig").unwrap();

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
    let flopp_dir = "/home/jshaw/practical_prob_2020_paper/flopp/";
    let frags_map = file_reader::get_frags_container(flopp_dir.to_owned() + "/tests/test_file.txt");
    let frags = frags_map.get("frag_contig").unwrap();
    let indexed_reads = utils_frags::get_all_overlaps(frags);
    let interval_reads_all  = local_clustering::find_reads_in_interval(1,100,frags, usize::MAX);
    let interval_reads_5= local_clustering::find_reads_in_interval(5,6,frags, usize::MAX);

    assert_eq!(interval_reads_all.len(),3);
    assert_eq!(interval_reads_5.len(),1);

}

