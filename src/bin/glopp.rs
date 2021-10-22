extern crate time;
use clap::{App, AppSettings, Arg};
use flopp::file_reader;
use flopp::global_clustering;
use flopp::local_clustering;
use flopp::types_structs::Frag;
use flopp::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use std::sync::Mutex;
use std::time::Instant;
fn main() {
    let matches = App::new("glopp")
                          .version("0.1.0")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("glopp - polyploid phasing from read sequencing.\n\nExample usage :\nglopp -b bamfile.bam -c vcffile.vcf -p (ploidy) -o results (with VCF and BAM)\nglopp -f fragfile.frags -o unpolished_results.txt (with fragment file and with automatic ploidy detection)")
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
                               .takes_value(true)
                               .hidden(true))
                          .arg(Arg::with_name("supp_aln_anchor")
                               .short("a")
                               .help("Input names of neighbouring contigs when only phasing one contig. The adjacent contigs in an assembly graph will be used to anchor the initial partition when using glopp. Usage: -a contig_1,contig_2")
                               .value_name("CONTIG1,CONTIG2")
                               .takes_value(true))
                          .arg(Arg::with_name("ploidy")
                              .short("p")
                              .help("Ploidy of organism. If not given, glopp will estimate the ploidy.")
                              .value_name("PLOIDY")
                              .takes_value(true))
                          .arg(Arg::with_name("threads")
                              .short("t")
                              .help("Number of threads to use (default : 10).")
                              .value_name("THREADS")
                              .takes_value(true))
                          .arg(Arg::with_name("output")
                              .short("xxxx")
                              .help("Name of output file (default : glopp_out.txt)")
                              .value_name("OUTPUT")
                              .takes_value(true)
                              .hidden(true))
                          .arg(Arg::with_name("partition output")
                              .short("o")
                              .help("Partition BAM based on read partition. (default : no BAM output. Specify file name when using -h.)")
                              .value_name("PARTITION OUTPUT BAM")
                              .takes_value(true))
                          .arg(Arg::with_name("epsilon")
                              .short("e")
                              .takes_value(true)
                              .help("Estimated allele call error rate. (default: 0.05. If using short reads, make sure to adjust this)"))
                          .arg(Arg::with_name("max_number_solns")
                              .short("n")
                              .takes_value(true)
                              .value_name("NUMBER OF SOLNS")
                              .help("Maximum number of solutions for beam search."))
                          .arg(Arg::with_name("use_mec")
                              .short("m")
                              .help("Use MEC score instead instead of a probabilistic objective function."))
                            .arg(Arg::with_name("use_ref_bias")
                              .short("R")
                              .help("Use reference bias adjustment (in progress)."))
                          .arg(Arg::with_name("verbose")
                              .short("r")
                              .help("Verbose output."))
                          .arg(Arg::with_name("dont_filter_supplementary")
                              .short("S")
                              .help("Use all supplementary alignments from the BAM file without filtering (filtering by default)."))
                          .get_matches();

    let mut estimate_ploidy = false;
    let max_number_solns_str = matches.value_of("max_number_solns").unwrap_or("10");
    let max_number_solns = max_number_solns_str.parse::<usize>().unwrap();
    let num_t_str = matches.value_of("threads").unwrap_or("10");
    let num_t = match num_t_str.parse::<usize>() {
        Ok(num_t) => num_t,
        Err(_) => panic!("Number of threads must be positive integer"),
    };

    let large_numb = 300;
    let ploidy = matches.value_of("ploidy").unwrap_or("300");
    let mut ploidy = ploidy.parse::<usize>().unwrap();
    if ploidy == large_numb {
        estimate_ploidy = true;
    }

    let use_mec = matches.is_present("use_mec");
    let use_ref_bias = matches.is_present("use_ref_bias");
    let filter_supplementary = !matches.is_present("dont_filter_supplementary");
    // Set up our logger if the user passed the debug flag
    if matches.is_present("verbose") {
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Trace)
            .init()
            .unwrap();
    } else {
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Debug)
            .init()
            .unwrap();
    }

    //If the user is splitting the bam file according to the output partition.
    let part_out_dir = matches
        .value_of("partition output")
        .unwrap_or("glopp_out_dir");

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
            panic!("Don't allow -v option for now.");
            vcf = true;
            vcf_file
        }
    };

    //Using supplemental alignments/assembly graph
    let use_supp_anchor;
    let contig_supp_string = match matches.value_of("supp_aln_anchor") {
        None => {
            use_supp_anchor = false;
            "_"
        }
        Some(contig_names) => {
            use_supp_anchor = true;
            contig_names
        }
    };

    let contig_anchors: Vec<String> = contig_supp_string
        .split(',')
        .into_iter()
        .map(|s| s.to_owned())
        .collect();
    if contig_anchors.len() < 2 && use_supp_anchor {
        panic!("Number of supplementary contig anchors must exceed 1. Exiting");
    }

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

    if vcf_nopolish {
        vcf_file = vcf_file_nopolish;
    }

    if vcf_nopolish && vcf {
        panic!("Only use one of the VCF options. -c if diploid VCF or choosing to polish, -v otherwise.\n");
    }

    if !bam && !frag{
        panic!("Must input a BAM file.")
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

    let _output_blocks_str = matches.value_of("output").unwrap_or("glopp_output.txt");

    if bam && frag {
        panic!("If using frag as input, BAM file should not be specified")
    }

    if bam && (!vcf && !vcf_nopolish) {
        panic!("Must input VCF file if using BAM file");
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_t)
        .build_global()
        .unwrap();

    //CONSTANTS - Constants which users probably should not change.

    let polish = vcf;
    //If we estimate the frag error rate by clustering a few random test blocks.
    //    let estimate_epsilon = true;

    println!("Reading inputs (BAM/VCF/frags).");
    let start_t = Instant::now();
    let mut all_frags_map;
    if bam {
        all_frags_map =
            file_reader::get_frags_from_bamvcf(vcf_file, bam_file, filter_supplementary);
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
            if polish {
                panic!("VCF File ploidy doesn't match input ploidy");
            }
        }
    }

    let first_iter = true;

    for (contig, all_frags) in all_frags_map.iter_mut() {
        let mut prev_expected_score = f64::MAX;
        println!("Number of fragments {}", all_frags.len());
        for frag in all_frags.iter() {
            if let Some(sup_cont) = &frag.supp_aln {
                log::trace!(
                    "{}, {}, {}, {}",
                    sup_cont,
                    frag.first_position,
                    frag.last_position,
                    frag.id
                );
            }
        }
        if snp_to_genome_pos_map.contains_key(contig) || bam == false {
            let mut _genotype_dict: &FxHashMap<usize, FxHashMap<usize, usize>> =
                &FxHashMap::default();
            let mut snp_to_genome_pos: &Vec<usize> = &Vec::new();

            if bam == true {
                snp_to_genome_pos = snp_to_genome_pos_map.get(contig).unwrap();
                if polish == true {
                    _genotype_dict = genotype_dict_map.get(contig).unwrap();
                }
            }
            //I think this is done because there is assumd to be only 1 contig,
            //so we just iterate over a container of size 1. Furthermore,
            //genotype_dict_map is empty in the case of using frags.
            else {
                for value in genotype_dict_map.values() {
                    _genotype_dict = value;
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
            //let mut binomial_factor = (avg_read_length as f64) / heuristic_multiplier;
            let mut _binomial_factor = 1.0;
            println!("Binomial adjustment factor is {}", _binomial_factor);

            //Final partitions
            let _parts: Mutex<Vec<(Vec<FxHashSet<&Frag>>, usize)>> = Mutex::new(vec![]);

            //The length of each local haplotype block. This is currently useless because all haplotype
            //blocks have the same fixed length, but we may change this in the future.
            let _lengths: Mutex<Vec<usize>> = Mutex::new(vec![]);

            //UPEM scores for each block.
            let _scores: Mutex<Vec<(f64, usize)>> = Mutex::new(vec![]);
            let length_block = utils_frags::get_avg_length(&all_frags, block_len_quant);

            //If we want blocks to overlap -- I don't think we actually want blocks to overlap but this may
            //be an optional parameter for testing purposes.
            let _overlap = 0;
            //let cutoff_value = (1.0 / (ploidy + 1) as f64).ln();
            let cutoff_value = f64::MIN;

            //Get last SNP on the genome covered over all fragments.
            let length_gn = utils_frags::get_length_gn(&all_frags);
            println!("Length of genome is {}", length_gn);
            println!("Length of each block is {}", length_block);
            //How many blocks we iterate through to estimate epsilon.
            let mut epsilon = 0.05;
            //            let num_iters = length_gn / length_block;
            //            if estimate_epsilon {
            //                epsilon = local_clustering::estimate_epsilon(
            //                    num_iters,
            //                    num_epsilon_attempts,
            //                    ploidy,
            //                    &all_frags,
            //                    length_block,
            //                    epsilon,
            //                );
            //            }

            //TEST TODO
            //epsilon = epsilon * ploidy as f64 / 4.0;

            if epsilon < 0.01 {
                epsilon = 0.010;
            }

            match matches.value_of("epsilon") {
                None => {}
                Some(value) => {
                    epsilon = value.parse::<f64>().unwrap();
                }
            };

            println!("Epsilon is {}", epsilon);

            let num_estimate_tries = 20;
            if estimate_ploidy {
                ploidy = local_clustering::estimate_ploidy(
                    length_gn,
                    num_estimate_tries,
                    &all_frags,
                    epsilon,
                );
            }
            println!("Ploidy is {}", ploidy);
            //Phasing occurs here
            let start_t = Instant::now();
            let initial_part;
            //If first_pos = last_pos, then the initial_part is empty and we rely on the beam
            //search to determine the correct initial partition.
            let first_pos = 1;
            let last_pos = 1;
            if !use_supp_anchor {
                initial_part = global_clustering::get_initial_clique(
                    all_frags, ploidy, epsilon, first_pos, last_pos,
                );
            } else {
                //                initial_part = global_clustering::get_initial_clique(all_frags, ploidy, epsilon);
                initial_part = global_clustering::get_initial_from_anchor(
                    all_frags,
                    ploidy,
                    epsilon,
                    &contig_anchors,
                );
            }
            let (break_positions, final_part) = global_clustering::beam_search_phasing(
                initial_part,
                all_frags,
                epsilon,
                _binomial_factor,
                cutoff_value,
                max_number_solns,
                use_mec,
                use_supp_anchor,
                use_ref_bias,
            );
            println!("Time taken for phasing {:?}", Instant::now() - start_t);

            let final_block_unpolish = utils_frags::hap_block_from_partition(&final_part);
            let (f_binom_vec, f_freq_vec) =
                local_clustering::get_partition_stats(&final_part, &final_block_unpolish);
            let final_score =
                -1.0 * local_clustering::get_mec_score(&f_binom_vec, &f_freq_vec, 0.0, 0.0);
            println!("Final MEC score for the partition is {:?}.", final_score);

            file_reader::write_output_partition_to_file(
                &final_part,
                part_out_dir,
                contig,
                &break_positions,
            );

            file_reader::write_blocks_to_file(
                part_out_dir,
                &vec![final_block_unpolish],
                &vec![length_gn],
                &snp_to_genome_pos,
                &final_part,
                first_iter,
                contig,
                &break_positions,
            );
        }
    }
}
