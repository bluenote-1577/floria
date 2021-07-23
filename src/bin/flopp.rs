extern crate time;
use clap::{App, AppSettings, Arg};
use flopp::file_reader;
use flopp::local_clustering;
use flopp::types_structs::Frag;
use flopp::types_structs::HapBlock;
use flopp::utils_frags;
use flopp::vcf_polishing;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use std::sync::Mutex;
use std::time::Instant;

fn main() {
    let matches = App::new("flopp")
                          .version("0.1.0")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("flopp - fast local polyploid phasing from long read sequencing.\n\nExample usage :\nflopp -b bamfile.bam -v vcffile.vcf -p (ploidy) -o results.txt (with VCF and BAM)\nflopp -f fragfile.frags -p (ploidy) -o unpolished_results.txt (with fragment file)\nflopp -f fragfile.frags -v vcffile.vcf -p (ploidy) -o polished_results.txt (polished output with fragment file and VCF)")
                          .arg(Arg::with_name("frag")
                               .short("f")
                               .value_name("FRAGFILE")
                               .help("Input a fragment file.")
                               .takes_value(true))
                          .arg(Arg::with_name("vcf no polish")
                              .short("c")
                              .value_name("VCFFILE")
                               .help("Input a VCF: Does NOT use polishing. Use this if your VCF is not polyploid.")
                                .takes_value(true))
                          .arg(Arg::with_name("bam")
                              .short("b")
                              .value_name("BAMFILE")
                               .help("Input a bam file.")
                                .takes_value(true))
                          .arg(Arg::with_name("vcf")
                               .short("v")
                               .help("Input a VCF: Mandatory if using BAM file; Enables genotype polishing if using frag file.")
                               .value_name("VCFFILE")
                               .takes_value(true))
                          .arg(Arg::with_name("ploidy")
                              .short("p")
                              .help("Ploidy of organism.")
                              .value_name("PLOIDY")
                              .required(true)
                              .takes_value(true))
                          .arg(Arg::with_name("threads")
                              .short("t")
                              .help("Number of threads to use (default : 10).")
                              .value_name("THREADS")
                              .takes_value(true))
                            .arg(Arg::with_name("range")
                              .short("r")
                              .hidden(true)
                              .help("Haplotyping of only variants in a specified range, inclusive (default : 1-[last position of chromosome]).")
                              .value_name("RANGE")
                              .takes_value(true))
                          .arg(Arg::with_name("output")
                              .short("o")
                              .help("Name of output file (default : flopp_out.txt)")
                              .value_name("OUTPUT")
                              .takes_value(true))
                          .arg(Arg::with_name("partition output")
                              .short("h")
                              .help("Partition BAM based on read partition. (default : no BAM output. Specify file name when using -h.)")
                              .value_name("PARTITION OUTPUT BAM")
                              .takes_value(true))
                          .arg(Arg::with_name("epsilon")
                              .short("e")
                              .takes_value(true)
                              .hidden(true))
                          .arg(Arg::with_name("binomial factor")
                              .short("s")
                              .takes_value(true)
                              .value_name("BINOMIAL FACTOR")
                              .help("The normalizing factor for UPEM (sigma in the paper)"))
                          .arg(Arg::with_name("use_mec")
                              .short("m")
                              .help("Use MEC score instead of UPEM for cluster refinement. Use this when your haplotypes have unbalanced coverage."))
                          .get_matches();

    let num_t_str = matches.value_of("threads").unwrap_or("10");
    let num_t = match num_t_str.parse::<usize>() {
        Ok(num_t) => num_t,
        Err(_) => panic!("Number of threads must be positive integer"),
    };
    let ploidy = matches.value_of("ploidy").unwrap();
    let ploidy = match ploidy.parse::<usize>() {
        Ok(ploidy) => ploidy,
        Err(_) => panic!("Must input valid ploidy"),
    };

    let use_mec = matches.is_present("use_mec");

    let heuristic_multiplier_str = matches.value_of("binomial factor").unwrap_or("25.0");
    let heuristic_multiplier = match heuristic_multiplier_str.parse::<f64>() {
        Ok(heuristic_multiplier) => heuristic_multiplier,
        Err(_) => panic!("Must input valid binomial normalization."),
    };

    //If the user is splitting the bam file according to the output partition.
    let bam_part_out;
    let bam_part_out_dir = match matches.value_of("partition output") {
        None => {
            bam_part_out = false;
            "_"
        }
        Some(bam_part_out_dir) => {
            bam_part_out = true;
            bam_part_out_dir
        }
    };

    //If the user is getting frag files from BAM and VCF.
    let bam;
    let bam_file = match matches.value_of("bam") {
        None => {
            bam = false;
            "_"
        }
        Some(bam_file) => {
            bam = true;
            bam_file
        }
    };

    //If user is using a frag file.
    let frag;
    let frag_file = match matches.value_of("frag") {
        None => {
            frag = false;
            "_"
        }
        Some(frag_file) => {
            frag = true;
            frag_file
        }
    };

    //Whether or not we polish using genotyping information from VCF.
    let vcf;
    let mut vcf_file = match matches.value_of("vcf") {
        None => {
            vcf = false;
            "_"
        }
        Some(vcf_file) => {
            vcf = true;
            vcf_file
        }
    };

    //Use a VCF without polishing.
    let vcf_nopolish;
    let vcf_file_nopolish = match matches.value_of("vcf no polish") {
        None => {
            vcf_nopolish = false;
            "_"
        }
        Some(vcf_file_nopolish) => {
            vcf_nopolish = true;
            vcf_file_nopolish
        }
    };

    if vcf_nopolish{
        vcf_file = vcf_file_nopolish;
    }

    if vcf_nopolish && vcf{
        panic!("Only use one of the VCF options. -c if diploid VCF or choosing to polish, -v otherwise.\n");
    }


    //Only haplotype variants in a certain range : TODO
    let _range;
    let _range_string = match matches.value_of("range") {
        None => {
            _range = false;
            "_"
        }
        Some(_range_string) => {
            _range = true;
            _range_string
        }
    };

    let output_blocks_str = matches.value_of("output").unwrap_or("flopp_output.txt");

    if bam && frag {
        panic!("If using frag as input, BAM file should not be specified")
    }

    if bam && (!vcf && !vcf_nopolish){
        panic!("Must input VCF file if using BAM file");
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_t)
        .build_global()
        .unwrap();

    //CONSTANTS - Constants which users probably should not change.

    //The inter quantile range outlier factor. 3.0 is standard for detecting extreme outliers
    //Blocks outside of the range will get filled in.
    let iqr_factor = 3.0;
    let polish = vcf;
    //If we fill blocks with poor UPEM score.
    let fill = true;
    //If we estimate the frag error rate by clustering a few random test blocks.
    let estimate_epsilon = true;
    //Number of iterations for the iterative UPEM optimization
    let num_iters_optimizing = 10;

    println!("Reading inputs (BAM/VCF/frags).");
    let start_t = Instant::now();
    let mut all_frags_map;
    if bam {
        all_frags_map = file_reader::get_frags_from_bamvcf(vcf_file, bam_file);
    } else {
        all_frags_map = file_reader::get_frags_container(frag_file);
    }
    println!("Time taken reading inputs {:?}", Instant::now() - start_t);

    let mut genotype_dict_map: FxHashMap<String, FxHashMap<usize, FxHashMap<usize, usize>>> =
        FxHashMap::default();
    let mut snp_to_genome_pos_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    if vcf || vcf_nopolish {
        let (snp_to_genome_pos_t, genotype_dict_t, vcf_ploidy) =
            file_reader::get_genotypes_from_vcf_hts(vcf_file);
        snp_to_genome_pos_map = snp_to_genome_pos_t;
        genotype_dict_map = genotype_dict_t;

        //If the VCF file is misformatted or has weird genotyping call we can catch that here.
        if vcf_ploidy != ploidy {
            if polish{
                panic!("VCF File ploidy doesn't match input ploidy");
            }
        }
    }

    let mut first_iter = true;

    for (contig, all_frags) in all_frags_map.iter_mut() {
        if snp_to_genome_pos_map.contains_key(contig) || bam == false {
            let mut genotype_dict: &FxHashMap<usize, FxHashMap<usize, usize>> =
                &FxHashMap::default();
            let mut snp_to_genome_pos: &Vec<usize> = &Vec::new();

            if bam == true {
                snp_to_genome_pos = snp_to_genome_pos_map.get(contig).unwrap();
                if polish == true{
                    genotype_dict = genotype_dict_map.get(contig).unwrap();
                }
            } 
            //I think this is done because there is assumd to be only 1 contig,
            //so we just iterate over a container of size 1. Furthermore, 
            //genotype_dict_map is empty in the case of using frags. 
            else {
                for value in genotype_dict_map.values() {
                    genotype_dict = value;
                }
                for value in snp_to_genome_pos_map.values() {
                    snp_to_genome_pos = value;
                }
            }

            //We need frags sorted by first position to make indexing easier.
            all_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));

            //We use the median # bases spanned by fragments as the length of blocks.
            let avg_read_length = utils_frags::get_avg_length(&all_frags, 0.5);
            println!("Median read length is {}", avg_read_length);

            let block_len_quant = 0.33;

            //The sample size correction factor for the binomial test used on the bases/errors.
            let binomial_factor = (avg_read_length as f64) / heuristic_multiplier;
            println!("Binomial adjustment factor is {}", binomial_factor);

            //Final partitions
            let parts: Mutex<Vec<(Vec<FxHashSet<&Frag>>, usize)>> = Mutex::new(vec![]);

            //The length of each local haplotype block. This is currently useless because all haplotype
            //blocks have the same fixed length, but we may change this in the future.
            let lengths: Mutex<Vec<usize>> = Mutex::new(vec![]);

            //UPEM scores for each block.
            let scores: Mutex<Vec<(f64, usize)>> = Mutex::new(vec![]);
            let length_block = utils_frags::get_avg_length(&all_frags, block_len_quant);

            //If we want blocks to overlap -- I don't think we actually want blocks to overlap but this may
            //be an optional parameter for testing purposes.
            let overlap = 0;

            //Get last SNP on the genome covered over all fragments.
            let length_gn = utils_frags::get_length_gn(&all_frags);
            println!("Length of genome is {}", length_gn);
            println!("Length of each block is {}", length_block);
            //How many blocks we iterate through to estimate epsilon.
            let num_epsilon_attempts = 20;
            let mut epsilon = 0.03;
            let num_iters = length_gn / length_block;
            if estimate_epsilon {
                epsilon = local_clustering::estimate_epsilon(
                    num_iters,
                    num_epsilon_attempts,
                    ploidy,
                    &all_frags,
                    length_block,
                    epsilon,
                );
            }

            //TEST TODO
            //epsilon = epsilon * ploidy as f64 / 4.0;

            if epsilon == 0.0 {
                epsilon = 0.010;
            }

            match matches.value_of("epsilon") {
                None => {}
                Some(value) => {
                    epsilon = value.parse::<f64>().unwrap();
                }
            };

            println!("Estimated epsilon is {}", epsilon);

            println!("Generating haplotype blocks");
            let start_t = Instant::now();

            //Embarassing parallel building of local haplotype blocks using rayon crate.
            (0..num_iters)
                .collect::<Vec<usize>>()
                .into_par_iter()
                .for_each(|x| {
                    //            println!("{} iteration number", x);
                    let block_start = x * length_block + 1;
                    let part = local_clustering::generate_hap_block(
                        block_start,
                        block_start + length_block + overlap,
                        ploidy,
                        &all_frags,
                        epsilon,
                    );

                    let (best_score, best_part, _best_block) =
                        local_clustering::optimize_clustering(
                            part,
                            epsilon,
                            &genotype_dict,
                            polish,
                            num_iters_optimizing,
                            binomial_factor,
                            use_mec,
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
                    epsilon,
                );
            } else {
                part_filled = parts;
            }
            println!("Time taken block filling {:?}", Instant::now() - start_t);

            let start_t = Instant::now();

            //Link and polish all blocks.
            //    let mut final_part = vcf_polishing::link_blocks(&part_filled);
            let mut final_part = vcf_polishing::link_blocks_greedy(&part_filled, &all_frags);
            //    let mut final_part = vcf_polishing::link_blocks_heur(&part_filled,4);

            //    for i in 0..ploidy{
            //        let inter =  &final_part[i].intersection(&final_part2[i]).collect::<Vec<_>>();
            //        dbg!(inter.len());
            //        let diff =  &final_part[i].symmetric_difference(&final_part2[i]).collect::<Vec<_>>();
            //        dbg!(diff.len());
            //        dbg!(final_part[i].len());
            //        dbg!(final_part2[i].len());
            //    }
            let final_block_unpolish = utils_frags::hap_block_from_partition(&final_part);
            let mut final_block_polish = HapBlock { blocks: Vec::new() };
            if polish {
                final_block_polish = vcf_polishing::polish_using_vcf(
                    &genotype_dict,
                    &final_block_unpolish,
                    &(1..length_gn + 1).collect::<Vec<_>>(),
                );
            }

            //TEST SECTION
            //
            vcf_polishing::remove_duplicate_reads(
                &mut final_part,
                &all_frags,
                &final_block_unpolish,
            );

            let (f_binom_vec, f_freq_vec) =
                local_clustering::get_partition_stats(&final_part, &final_block_unpolish);
            let final_score = local_clustering::get_mec_score(&f_binom_vec, &f_freq_vec, 0.0, 0.0);
            println!(
                "Final MEC score for the partition is {:?}.",
                -1.0 * final_score
            );

            //    vcf_polishing::map_reads_against_hap_errors(&final_part,&final_block_polish,epsilon);
            //    let final_block_unpolish = utils_frags::hap_block_from_partition(&final_part);
            //    let mut final_block_polish = HapBlock { blocks: Vec::new() };
            //    if polish {
            //        final_block_polish = vcf_polishing::polish_using_vcf(
            //            &genotype_dict,
            //            &final_block_unpolish,
            //            &(1..length_gn + 1).collect::<Vec<_>>(),
            //        );
            //    }
            //

            println!(
                "Time taken linking, polishing blocks {:?}",
                Instant::now() - start_t
            );

            let start_t = Instant::now();
            //Write blocks to file. Write fragments to a file if the user chooses to do that instead.
            if polish {
                file_reader::write_blocks_to_file(
                    output_blocks_str,
                    &vec![final_block_polish],
                    &vec![length_gn],
                    &snp_to_genome_pos,
                    &final_part,
                    first_iter,
                    contig,
                );
            } else {
                file_reader::write_blocks_to_file(
                    output_blocks_str,
                    &vec![final_block_unpolish],
                    &vec![length_gn],
                    &snp_to_genome_pos,
                    &final_part,
                    first_iter,
                    contig,
                );
            }

            if bam_part_out {
                file_reader::write_output_partition_to_file(&final_part, bam_part_out_dir, contig);
            }

            first_iter = false;

            println!(
                "Time taken writing blocks for contig {} to {} : {:?}",
                contig,
                output_blocks_str,
                Instant::now() - start_t
            );
        }
    }
}
