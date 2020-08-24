extern crate time;
use rayon::prelude::*;
use std::sync::{Mutex};
use haplotype_phaser::file_reader;
use haplotype_phaser::local_clustering;
use haplotype_phaser::utils_frags;
use haplotype_phaser::vcf_polishing;
use haplotype_phaser::types_structs::{HapBlock,Frag};
use fxhash::FxHashSet;
use std::env;
use std::time::Instant;

fn main() {
    let num_t = 20;
    rayon::ThreadPoolBuilder::new().num_threads(num_t).build_global().unwrap();
    let ploidy = 4;
    let iqr_factor = 1.5;
    let polish = true;
    let fill = true;
    let greedy = false;
    let bam = false;
    let estimate_epsilon = true;
    let use_manual_params = false;

    let args: Vec<String> = env::args().collect();
    println!("Reading frags");
    let all_frags;
    if bam{
        all_frags = file_reader::get_frags_from_bamvcf(&args[2],&args[1]);
    }
    else{
        all_frags = file_reader::get_frags_container(&args[1]);
    }
    let (genotype_dict, vcf_ploidy) = file_reader::get_genotypes_from_vcf(&args[2]);

    if vcf_ploidy != ploidy {
        panic!("VCF File ploidy doesn't match input ploidy");
    }

    let avg_read_length = utils_frags::get_avg_length(&all_frags);
    println!("Median read length is {}",avg_read_length);
    let binomial_factor = (avg_read_length as f64)/(ploidy as f64);
    let parts: Mutex<Vec<(Vec<FxHashSet<&Frag>>,usize)>> = Mutex::new(vec![]);
    let lengths: Mutex<Vec<usize>> = Mutex::new(vec![]);
    let scores: Mutex<Vec<(f64,usize)>> = Mutex::new(vec![]);
    let mut length_block = avg_read_length;
    let overlap = 0;
    let length_gn = utils_frags::get_length_gn(&all_frags);
    println!("Length of genome is {}",length_gn);
    let num_epsilon_attempts = 50/ploidy;
    let mut epsilon = 0.03;
    if use_manual_params{
        length_block = 100;
    }
    let num_iters = length_gn / length_block;
    if estimate_epsilon{
        epsilon = local_clustering::estimate_epsilon(num_iters,num_epsilon_attempts,ploidy,&all_frags,length_block);
    }
    if use_manual_params{
        epsilon = 0.015;
    }

    

    println!("Estimated epsilon is {}",epsilon);
//    let num_iters = 100;

    println!("Generating haplotype blocks");
    let start_t = Instant::now();

    (0..num_iters).collect::<Vec<usize>>().into_par_iter().for_each(|x| {
        println!("{} iteration number", x);
        let block_start = x*length_block + 1;

        let part = local_clustering::generate_hap_block(
            block_start,
            block_start + length_block + overlap,
            ploidy,
            &all_frags,
        );

        let (best_score, best_part, _best_block) =
            local_clustering::optimize_clustering(part, epsilon, &genotype_dict,polish , 10, binomial_factor);

        let mut locked_parts = parts.lock().unwrap();
        let mut locked_lengths = lengths.lock().unwrap();
        let mut locked_scores = scores.lock().unwrap();
        locked_parts.push((best_part,x));
        locked_lengths.push(length_block);
        locked_scores.push((best_score,x));

//        println!("UPEM Score : {}", best_score);
    });

    println!("Time taken local clustering {:?}", Instant::now() - start_t);

    let mut scores = scores.lock().unwrap().to_vec();
    scores.sort_by(|a,b| a.1.cmp(&b.1));
    let scores = scores.into_iter().map(|x| x.0).collect();
    let mut parts = parts.lock().unwrap().to_vec();
    parts.sort_by(|a,b| a.1.cmp(&b.1));
    let parts = parts.into_iter().map(|x| x.0).collect();

    let start_t = Instant::now();
    let part_filled;

    if fill{
        part_filled = vcf_polishing::replace_with_filled_blocks(&scores,parts,iqr_factor,length_block,&all_frags);
    }
    else{
        part_filled = parts;
    }
    println!("Time taken block filling {:?}", Instant::now() - start_t);

    let start_t = Instant::now();

    let final_part = vcf_polishing::link_blocks(&part_filled,greedy);
    let mut final_block = utils_frags::hap_block_from_partition(&final_part);
    if polish{
        final_block = vcf_polishing::polish_using_vcf(&genotype_dict,&final_block,&(1..length_gn+1).collect::<Vec<_>>());
    }

    println!("Time taken linking, polishing blocks {:?}", Instant::now() - start_t);

    file_reader::write_blocks_to_file("test.txt", &vec![final_block], &vec![length_gn]);
}
