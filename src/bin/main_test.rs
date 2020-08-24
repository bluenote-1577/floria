use fxhash::{FxHashMap, FxHashSet};
use haplotype_phaser::file_reader;
use haplotype_phaser::local_clustering;
use haplotype_phaser::utils_frags;
use std::collections::{BTreeMap, HashSet,HashMap};
use std::time::Instant;

fn main() {
    let frags =
        file_reader::get_frags_container("/home/jshaw/haplotype_phaser/tests/test_file_long.txt");
    let mut hashset = FxHashSet::default();
    let mut hashmap = HashMap::new();

    let s = &frags[0];
    let start_t = Instant::now();
    hashset.insert(s);
    hashset.insert(s);
    hashset.insert(s);
    println!(
        "Time taken to insert 2 set frag {:?}",
        Instant::now() - start_t
    );

    let start_t = Instant::now();
    hashset.insert(s);
    println!(
        "Time taken to insert already frag {:?}",
        Instant::now() - start_t
    );



    let start_t = Instant::now();
    hashmap.insert(0, 1);
    hashmap.insert(2, 1);
    println!(
        "Time taken to insert 2 map int {:?}",
        Instant::now() - start_t
    );

    let start_t = Instant::now();
    hashset.get(&frags[0]).unwrap();
    println!("Time taken to get frag {:?}", Instant::now() - start_t);

    let start_t = Instant::now();
    hashset.contains(&frags[0]);
    println!("Time taken to check frag {:?}", Instant::now() - start_t);

    let mut hashset_int = FxHashSet::default();

    let start_t = Instant::now();
    hashset_int.insert(0);
    hashset_int.insert(0);
    hashset_int.insert(0);
    println!("Time taken to insert 2int {:?}", Instant::now() - start_t);

    let start_t = Instant::now();
    hashset_int.get(&0).unwrap();
    println!("Time taken to get int {:?}", Instant::now() - start_t);

    let start_t = Instant::now();
    utils_frags::distance(&frags[2], &frags[1]);
    println!("Time taken to get distance {:?}", Instant::now() - start_t);
}
