extern crate time;
use fxhash::FxHashSet;
use haplotype_phaser::file_reader;
use haplotype_phaser::local_clustering;
use haplotype_phaser::types_structs::{Frag, HapBlock};
use haplotype_phaser::utils_frags;
use haplotype_phaser::vcf_polishing;
use rayon::prelude::*;
use std::env;
use std::sync::Mutex;
use std::time::Instant;

//TODO : Parse commandline arguments in a good way.
fn main() {
    let num_t = 4;
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_t)
        .build_global()
        .unwrap();
    let ploidy = 4;
    //The inter quantile range outlier factor. 1.5 is standard and corresponds to 3 stds for a
    //normal. Blocks outside of the range will get filled in.
    let iqr_factor = 1.5;
    //Whether or not we polish using genotyping information from VCF.
    let polish = true;
    //If we fill blocks with poor UPEM score.
    let fill = true;
    //If we use greedy merging of blocks or enumerate all permutations.
    let greedy = false;
    //If the user is getting frag files from BAM and VCF.
    let bam = false;
    //If we estimate the frag error rate by clustering a few random test blocks.
    let estimate_epsilon = true;
    //TODO this is for testing
    let use_manual_params = false;
    //Number of iterations for the iterative UPEM optimization
    let num_iters_optimizing = 10;
    let output_blocks_str = "test.txt";
    let output_frag_str = "fragtest.txt";

    let args: Vec<String> = env::args().collect();
    println!("Reading frags");
    let mut all_frags;
    if bam {
        all_frags = file_reader::get_frags_from_bamvcf(&args[2], &args[1]);
    } else {
        all_frags = file_reader::get_frags_container(&args[1]);
    }

    //We need frags sorted by first position to make indexing easier.
    all_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));
    let (genotype_dict, vcf_ploidy) = file_reader::get_genotypes_from_vcf(&args[2]);

    //If the VCF file is misformatted or has weird genotyping call we can catch that here. 
    if vcf_ploidy != ploidy {
        panic!("VCF File ploidy doesn't match input ploidy");
    }

    //We use the median # bases spanned by fragments as the length of blocks.
    let avg_read_length = utils_frags::get_avg_length(&all_frags);
    println!("Median read length is {}", avg_read_length);
    //The sample size correcton factor for the binomial test used on the bases/errors.
    let binomial_factor = (avg_read_length as f64) / (ploidy as f64);
    //Final partitions
    let parts: Mutex<Vec<(Vec<FxHashSet<&Frag>>, usize)>> = Mutex::new(vec![]);
    //The length of each local haplotype block. This is currently useless because all haplotype
    //blocks have the same fixed length, but we may change this in the future. 
    let lengths: Mutex<Vec<usize>> = Mutex::new(vec![]);
    //UPEM scores for each block.
    let scores: Mutex<Vec<(f64, usize)>> = Mutex::new(vec![]);
    let mut length_block = avg_read_length;
    //If we want blocks to overlap -- I don't think we actually want blocks to overlap but this may
    //be an optional parameter. TODO
    let overlap = 0;
    //Get last SNP on the genome covered over all fragments.
    let length_gn = utils_frags::get_length_gn(&all_frags);
    println!("Length of genome is {}", length_gn);
    //How many blocks we iterate through to estimate epsilon.
    let num_epsilon_attempts = 50 / ploidy;
    //TODO This is for testing, clean this up.
    let mut epsilon = 0.03;
    if use_manual_params {
        length_block = 100;
    }
    let num_iters = length_gn / length_block;
    if estimate_epsilon {
        epsilon = local_clustering::estimate_epsilon(
            num_iters,
            num_epsilon_attempts,
            ploidy,
            &all_frags,
            length_block,
        );
    }
    if use_manual_params {
        epsilon = 0.015;
    }

    println!("Estimated epsilon is {}", epsilon);

    println!("Generating haplotype blocks");
    let start_t = Instant::now();

    //Embarassing parallel building of local haplotype blocks using rayon crate.
    (0..num_iters)
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|x| {
            println!("{} iteration number", x);
            let block_start = x * length_block + 1;
            let part = local_clustering::generate_hap_block(
                block_start,
                block_start + length_block + overlap,
                ploidy,
                &all_frags,
            );

            let (best_score, best_part, _best_block) = local_clustering::optimize_clustering(
                part,
                epsilon,
                &genotype_dict,
                polish,
                num_iters_optimizing,
                binomial_factor,
            );

            let mut locked_parts = parts.lock().unwrap();
            let mut locked_lengths = lengths.lock().unwrap();
            let mut locked_scores = scores.lock().unwrap();
            locked_parts.push((best_part, x));
            locked_lengths.push(length_block);
            locked_scores.push((best_score, x));

            //        println!("UPEM Score : {}", best_score);
        });

    println!("Time taken local clustering {:?}", Instant::now() - start_t);

    //Sort the vectors which may be out of order due to multi-threading.  
    let mut scores = scores.lock().unwrap().to_vec();
    scores.sort_by(|a, b| a.1.cmp(&b.1));
    let scores = scores.into_iter().map(|x| x.0).collect();
    let mut parts = parts.lock().unwrap().to_vec();
    parts.sort_by(|a, b| a.1.cmp(&b.1));
    let parts = parts.into_iter().map(|x| x.0).collect();

    let start_t = Instant::now();
    let part_filled;

    //Fill blocks
    if fill {
        part_filled = vcf_polishing::replace_with_filled_blocks(
            &scores,
            parts,
            iqr_factor,
            length_block,
            &all_frags,
        );
    } else {
        part_filled = parts;
    }
    println!("Time taken block filling {:?}", Instant::now() - start_t);

    let start_t = Instant::now();

    //Link and polish all blocks. 
    let final_part = vcf_polishing::link_blocks(&part_filled, greedy);
    let mut final_block = utils_frags::hap_block_from_partition(&final_part);
    if polish {
        final_block = vcf_polishing::polish_using_vcf(
            &genotype_dict,
            &final_block,
            &(1..length_gn + 1).collect::<Vec<_>>(),
        );
    }

    println!(
        "Time taken linking, polishing blocks {:?}",
        Instant::now() - start_t
    );

    let start_t = Instant::now();
    //Write blocks to file. Write fragments to a file if the user chooses to do that instead. 
    file_reader::write_blocks_to_file(output_blocks_str, &vec![final_block], &vec![length_gn]);
    file_reader::write_frags_file(all_frags, output_frag_str.to_string());
    println!(
        "Time taken writing blocks to {} and fragments to {} : {:?}",
        output_blocks_str,
        output_frag_str,
        Instant::now() - start_t
    );
}
