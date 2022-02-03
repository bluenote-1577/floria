extern crate time;
use clap::{App, AppSettings, Arg};
use flopp::file_reader;
use flopp::global_clustering;
use flopp::graph_processing;
use flopp::local_clustering;
use flopp::types_structs::Frag;
use flopp::utils_frags;
use fxhash::{FxHashMap};
use std::sync::Mutex;
use std::time::Instant;
fn main() {
    let matches = App::new("glopp")
                          .version("0.1.0")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("glopp - polyploid phasing from read sequencing.\n\nExample usage :\nglopp -b bamfile.bam -c vcffile.vcf -o results \n")
                          .arg(Arg::with_name("frag")
                               .short("f")
                               .value_name("FILE")
                               .help("Input a fragment file.")
                               .hidden(true)
                               .takes_value(true))
                          .arg(Arg::with_name("vcf no polish")
                              .short("c")
                              .value_name("FILE")
                               .help("Input a VCF.")
                                .takes_value(true))
                          .arg(Arg::with_name("bam")
                              .short("b")
                              .value_name("FILE")
                               .help("Input a bam file.")
                                .takes_value(true))
                          .arg(Arg::with_name("vcf")
                               .short("v")
                               .help("Input a VCF: Mandatory if using BAM file; Enables genotype polishing if using frag file.")
                               .value_name("FILE")
                               .takes_value(true)
                               .hidden(true))
                          .arg(Arg::with_name("ploidy") //Useful for testing. 
                              .short("p")
                              .help("Ploidy of organism. If not given, glopp will estimate the ploidy.")
                              .value_name("INT")
                              .takes_value(true)
                              .hidden(true))
                          .arg(Arg::with_name("threads")
                              .short("t")
                              .help("Number of threads to use. Not implemented yet. (default: 10).")
                              .value_name("INT")
                              .takes_value(true))
                          .arg(Arg::with_name("partition output")
                              .short("o")
                              .help("Output folder. (default: glopp_out_dir)")
                              .value_name("STRING")
                              .takes_value(true))
                          .arg(Arg::with_name("epsilon")
                              .short("e")
                              .takes_value(true)
                              .value_name("FLOAT")
                              .help("Estimated allele call error rate. (default: 0.04. If using short reads, make sure to adjust this)"))
                          .arg(Arg::with_name("max_number_solns")
                              .short("n")
                              .takes_value(true)
                              .value_name("INT")
                              .help("Maximum number of solutions for beam search. (default: 10)"))
                          .arg(Arg::with_name("num_iters_ploidy_est")
                              .short("q")
                              .takes_value(true)
                              .value_name("NUMBER BLOCKS")
                              .hidden(true)
                              .help("The number of blocks for flow graph construction when using fragments. (default 10)"))
                          .arg(Arg::with_name("bam_block_length")
                              .short("l")
                              .takes_value(true)
                              .value_name("INT")
                              .help("Length of blocks (in nucleotides) for flow graph construction when using bam file. (default: 15000)"))
                          .arg(Arg::with_name("dont_use_mec")
                              .short("u")
                              .help("")
                              .hidden(true))
                            .arg(Arg::with_name("use_ref_bias")
                              .short("R")
                              .hidden(true)
                              .help("Use reference bias adjustment (in progress)."))
                          .arg(Arg::with_name("verbose")
                              .short("r")
                              .help("Verbose output."))
                          .arg(Arg::with_name("dont_filter_supplementary")
                              .short("S")
                              .help("Use all supplementary alignments from the BAM file without filtering (filtering by default if using supp alignments)."))
                          .arg(Arg::with_name("use_supplementary")
                              .short("X")
                              .help("Use supplementary alignments (default: don't use; have not tested fully yet)."))
                          .get_matches();

    //Parse command line args.
    let max_number_solns_str = matches.value_of("max_number_solns").unwrap_or("10");
    let max_number_solns = max_number_solns_str.parse::<usize>().unwrap();
    let num_t_str = matches.value_of("threads").unwrap_or("10");
    let num_t = match num_t_str.parse::<usize>() {
        Ok(num_t) => num_t,
        Err(_) => panic!("Number of threads must be positive integer"),
    };

    let mut estimate_ploidy = false;
    let large_numb = 300;
    let ploidy = matches.value_of("ploidy").unwrap_or("300");
    let ploidy = ploidy.parse::<usize>().unwrap();
    if ploidy == large_numb {
        estimate_ploidy = true;
    }

    let block_length = matches.value_of("bam_block_length").unwrap_or("15000");
    let block_length = block_length.parse::<usize>().unwrap();
    //    let use_mec = matches.is_present("use_mec");
    let use_mec = true;
    let use_ref_bias = matches.is_present("use_ref_bias");
    let filter_supplementary = !matches.is_present("dont_filter_supplementary");
    let use_supplementary = matches.is_present("use_supplementary");

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
        .unwrap_or("glopp_out_dir").to_string();

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
        Some(_vcf_file) => {
            panic!("Don't allow -v option for now.");
            vcf = true;
            _vcf_file
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

    if vcf_nopolish {
        vcf_file = vcf_file_nopolish;
    }

    if vcf_nopolish && vcf {
        panic!("Only use one of the VCF options. -c if diploid VCF or choosing to polish, -v otherwise.\n");
    }

    if !bam && !frag {
        panic!("Must input a BAM file.")
    }

    //Only haplotype variants in a certain range : TODO not implemented yet
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

    //Polishing is a relic from flopp... I won't remove for now because it may be useful for
    //testing.
    let polish = vcf;

    println!("Reading inputs (BAM/VCF/frags).");
    let start_t = Instant::now();
    let mut all_frags_map;
    if bam {
        all_frags_map = file_reader::get_frags_from_bamvcf(
            vcf_file,
            bam_file,
            filter_supplementary,
            use_supplementary,
        );
    } else {
        all_frags_map = file_reader::get_frags_container(frag_file);
    }
    println!("Time taken reading inputs {:?}", Instant::now() - start_t);

    let mut snp_to_genome_pos_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    if vcf || vcf_nopolish {
        let (snp_to_genome_pos_t, _genotype_dict_t, vcf_ploidy) =
            file_reader::get_genotypes_from_vcf_hts(vcf_file);
        snp_to_genome_pos_map = snp_to_genome_pos_t;

        //If the VCF file is misformatted or has weird genotyping call we can catch that here.
        if vcf_ploidy != ploidy {
            if polish {
                panic!("VCF File ploidy doesn't match input ploidy");
            }
        }
    }

    let first_iter = true;

    for (contig, all_frags) in all_frags_map.iter_mut() {
        if all_frags.len() == 0 {
            println!("Contig {} has no fragments", contig);
            continue;
        }

        println!("Number of fragments {}", all_frags.len());
        if snp_to_genome_pos_map.contains_key(contig) || bam == false {
            let contig_out_dir = format!("{}/{}", part_out_dir, contig);
            let mut snp_to_genome_pos: &Vec<usize> = &Vec::new();

            if bam == true {
                snp_to_genome_pos = snp_to_genome_pos_map.get(contig).unwrap();
            }

            //We need frags sorted by first position to make indexing easier. We want the
            //counter_id to reflect the position in the vector.
            all_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));
            for (i, frag) in all_frags.iter_mut().enumerate() {
                frag.counter_id = i;
            }

            //We use the median # bases spanned by fragments as the length of blocks.
            let avg_read_length = utils_frags::get_avg_length(&all_frags, 0.5);
            println!("Median read length is {}", avg_read_length);

            let block_len_quant = 0.33;

            //Final partitions
            //let _parts: Mutex<Vec<(Vec<FxHashSet<&Frag>>, usize)>> = Mutex::new(vec![]);

            //The length of each local haplotype block. This is currently useless because all haplotype
            //blocks have the same fixed length, but we may change this in the future.
            //let _lengths: Mutex<Vec<usize>> = Mutex::new(vec![]);

            //UPEM scores for each block.
            //let _scores: Mutex<Vec<(f64, usize)>> = Mutex::new(vec![]);
            let _median_length_block = utils_frags::get_avg_length(&all_frags, block_len_quant);

            //let cutoff_value = (1.0 / (ploidy + 1) as f64).ln();
            let cutoff_value = f64::MIN;

            //Get last SNP on the genome covered over all fragments.
            let length_gn = utils_frags::get_length_gn(&all_frags);
            println!("Length of genome is {}", length_gn);
            println!("Length of each block is {}", block_length);
            let mut epsilon = 0.04;

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

            if estimate_ploidy {
                let num_locs_string = matches.value_of("num_iters_ploidy_est").unwrap_or("10");
                let num_locs = num_locs_string.parse::<usize>().unwrap();
                let mut hap_graph = graph_processing::generate_hap_graph(
                    length_gn,
                    num_locs,
                    &all_frags,
                    epsilon,
                    &snp_to_genome_pos,
                    max_number_solns,
                    block_length,
                    contig_out_dir.to_string(),
                );
                let flow_up_vec = graph_processing::solve_lp_graph(&hap_graph, contig_out_dir.to_string());
                graph_processing::get_disjoint_paths_rewrite(
                    &mut hap_graph,
                    flow_up_vec,
                    epsilon,
                    contig_out_dir.to_string(),
                    &snp_to_genome_pos,
                );

            } 
            //We don't actually use this path anymore, but it can be useful for testing purposes. 
            else {
                println!("Ploidy is {}", ploidy);
                //Phasing occurs here
                let start_t = Instant::now();
                let initial_part;
                //If first_pos = last_pos, then the initial_part is empty and we rely on the beam
                //search to determine the correct initial partition.
                let first_pos = 1;
                let last_pos = 1;
                initial_part = global_clustering::get_initial_clique(
                    all_frags, ploidy, epsilon, first_pos, last_pos,
                );
                let binom_factor = 1.;
                let all_frags_refs: Vec<&Frag> = all_frags.iter().collect();
                let (break_positions, final_part) = global_clustering::beam_search_phasing(
                    initial_part,
                    &all_frags_refs,
                    epsilon,
                    binom_factor,
                    cutoff_value,
                    max_number_solns,
                    use_mec,
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
                    vec![],
                    contig_out_dir.to_string(),
                    contig,
                    &snp_to_genome_pos
                );

                file_reader::write_blocks_to_file(
                    contig_out_dir.to_string(),
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
}
