use haplotype_phaser::file_reader;
use haplotype_phaser::local_clustering;
use haplotype_phaser::utils_frags;
use fnv::FnvHashSet;

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
    let mut hashset = FnvHashSet::default();
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
    let interval_reads_all  = local_clustering::find_reads_in_interval(1,100,&indexed_reads);
    let interval_reads_5= local_clustering::find_reads_in_interval(5,6,&indexed_reads);

    assert_eq!(interval_reads_all.len(),3);
    assert_eq!(interval_reads_5.len(),1);

}
