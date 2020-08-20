extern crate time;
use haplotype_phaser::file_reader;
use haplotype_phaser::local_clustering;
use haplotype_phaser::utils_frags;
use haplotype_phaser::vcf_polishing;
use std::env;
use std::time::Instant;

fn main() {
    let ploidy = 4;
    let epsilon = 0.03;
    let iqr_factor = 1.5;
    let polish = true;
    let fill = true;
    let greedy = false;

    let args: Vec<String> = env::args().collect();
    println!("Reading frags");
    let all_frags = file_reader::get_frags_container(&args[1]);
    let (genotype_dict, vcf_ploidy) = file_reader::get_genotypes_from_vcf(&args[2]);

    if vcf_ploidy != ploidy {
        panic!("VCF File ploidy doesn't match input ploidy");
    }

    let mut blocks = Vec::new();
    let mut parts = Vec::new();
    let mut lengths = Vec::new();
    let mut scores = Vec::new();
    let length_block = 100;
    let overlap = 0;
    let num_iters = 44900/length_block;
//    let num_iters = 100;

    println!("Generating haplotype blocks");
    let mut block_start = 1;
    let start_t = Instant::now();

    for _i in 0..num_iters {
        println!("{} iteration number", _i);

        let part = local_clustering::generate_hap_block(
            block_start,
            block_start + length_block + overlap,
            ploidy,
            &all_frags,
        );

        let (best_score, best_part, best_block) =
            local_clustering::optimize_clustering(part, epsilon, &genotype_dict,polish , 10);

        let vector_to_polish: Vec<usize> = (block_start..block_start + length_block+overlap+1).collect();
        let polished_block =
            vcf_polishing::polish_using_vcf(&genotype_dict, &best_block, &vector_to_polish);

        block_start += length_block;

        if polish{
            blocks.push(polished_block);
        }
        else{
            blocks.push(best_block);
        }
        parts.push(best_part);
        lengths.push(length_block);
        scores.push(best_score);

        println!("UPEM Score : {}", best_score);
    }
    println!("Time taken local clustering {:?}", Instant::now() - start_t);

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
        final_block = vcf_polishing::polish_using_vcf(&genotype_dict,&final_block,&(1..num_iters*(length_block)).collect::<Vec<_>>());
    }

    println!("Time taken linking, polishing blocks {:?}", Instant::now() - start_t);

    file_reader::write_blocks_to_file("test.txt", &vec![final_block], &vec![num_iters*length_block]);
}
