extern crate time;
use clap::{AppSettings, Arg, Command};
use fxhash::FxHashMap;
use sheaf::file_reader;
use sheaf::graph_processing;
use sheaf::utils_frags;
use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::time::Instant;
use std::fs;

#[allow(deprecated)]
fn main() {
    let input_options = "INPUT";
    let output_options = "OUTPUT";
    let alg_options = "ALGORITHM";
    let mandatory_options = "REQUIRED";
    let matches = Command::new("glopp")
                          .version("0.1.0")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("glopp - polyploid phasing from read sequencing.\n\nExample usage :\nglopp -b bamfile.bam -c vcffile.vcf -o results \n")
                          .arg(Arg::new("frag")
                               .short('f')
                               .value_name("FILE")
                               .help("Input a fragment file.")
                               .hide(true)
                               .takes_value(true))
                          .arg(Arg::new("bam")
                              .short('b')
                              .value_name("BAM FILE")
                              .required(true)
                              .help("Indexed and sorted primary bam file.")
                              .takes_value(true)
                              .help_heading(mandatory_options)
                              .display_order(1))
                          .arg(Arg::new("vcf")
                              .short('c')
                              .value_name("VCF FILE")
                              .required(true)
                              .help("A VCF file; must have contig header information present. See README for generating a VCF file with such information.")
                              .takes_value(true)
                              .help_heading(mandatory_options))
                          .arg(Arg::new("reference_fasta")
                              .short('R')
                              .takes_value(true)
                              .value_name("FASTA FILE")
                              .help("RECOMMENDED: Improve calls by realigning onto a reference with alternate alleles.")
                              .help_heading(input_options)
                              .display_order(1))
                          .arg(Arg::new("ploidy") //Useful for testing. 
                              .short('p')
                              .help("Ploidy of organism. If not given, glopp will estimate the ploidy.")
                              .value_name("INT")
                              .takes_value(true)
                              .hide(true))
                          .arg(Arg::new("threads")
                              .short('t')
                              .help("Number of threads to use. (default: 10).")
                              .value_name("INT")
                              .takes_value(true)
                              )
                          .arg(Arg::new("partition output")
                              .short('o')
                              .help("Output folder. (default: glopp_out_dir)")
                              .value_name("STRING")
                              .takes_value(true)
                              .help_heading(output_options))
                          .arg(Arg::new("epsilon")
                              .short('e')
                              .takes_value(true)
                              .value_name("FLOAT")
                              .help("Estimated allele call error rate. (default: 0.04 without -H, 0.02 with -H. If using short reads, make sure to adjust this)")
                              .help_heading(alg_options)
                              .display_order(1))
                          .arg(Arg::new("max_number_solns")
                              .short('n')
                              .takes_value(true)
                              .value_name("INT")
                              .help("Maximum number of solutions for beam search. Increasing may improve accuracy slightly. (default: 10)")
                              .help_heading(alg_options))
                          .arg(Arg::new("num_iters_ploidy_est")
                              .short('q')
                              .takes_value(true)
                              .value_name("NUMBER BLOCKS")
                              .hide(true)
                              .help("The number of blocks for flow graph construction when using fragments. (default 10)"))
                          .arg(Arg::new("bam_block_length")
                              .short('l')
                              .takes_value(true)
                              .value_name("INT")
                              .help("Length of blocks (in nucleotides) for flow graph construction when using bam file. (default: 15000)")
                              .help_heading(alg_options)
                              .display_order(1))
                          .arg(Arg::new("dont_use_mec")
                              .short('u')
                              .help("")
                              .hide(true))
                          .arg(Arg::new("verbose")
                              .short('r')
                              .help("Verbose output."))
                          .arg(Arg::new("snp_density")
                              .short('d')
                              .takes_value(true)
                              .value_name("FLOAT")
                              .help("Minimum SNP density to phase. Blocks with SNP density less than this value will not be phased. (Default: 0.001 i.e. 1 SNP per 1000 bases)")
                              .help_heading(alg_options))
                          .arg(Arg::new("use_supplementary")
                              .short('X')
                              .help("Use supplementary alignments (default: don't use).")
                              .help_heading(input_options))
                          .arg(Arg::new("hybrid")
                              .short('H')
                              .takes_value(true)
                              .value_name("BAM FILE")
                              .help("RECOMMENDED: Use short aligned short reads to polish long-read SNPs.")
                              .help_heading(input_options)
                              .display_order(1))
                          .arg(Arg::new("reassign_short")
                              .long("reassign-short")
                              .help("Reassign short reads when using the -H option to the best haplotigs. (Default: no reassignment)")
                              .help_heading(output_options))
                          .arg(Arg::new("do_binning")
                              .long("bin-by-cov")
                              .help("Increase contiguity by binning haplogroups using coverage (testing in progress).")
                              .help_heading(alg_options)
                              .display_order(2))
                          .arg(Arg::new("extend_read_clipping")
                              .long("extend-trimming")
                              .help("Trim less carefully against the reference for supplementary reads. May allow downstream assembly of sequences absent from reference with tradeoff for assembly quality. (testing in progress).")
                              .help_heading(output_options))
                          .arg(Arg::new("use_gaps")
                              .hide(true)
                              .long("use-gaps")
                              .help("Use gap information between SNPs while phasing (testing in progress)."))
                          .arg(Arg::new("list_to_phase")
                              .short('G')
                              .multiple(true)
                              .value_name("LIST OF CONTIGS")
                              .takes_value(true)
                              .help("Phase only contigs in this argument. Usage: -G contig1 contig2 contig3 ...")
                              .help_heading(input_options))
                          .get_matches();

    //Parse command line args.
    let max_number_solns_str = matches.value_of("max_number_solns").unwrap_or("10");
    let max_number_solns = max_number_solns_str.parse::<usize>().unwrap();
    let num_t_str = matches.value_of("threads").unwrap_or("10");
    let num_t = match num_t_str.parse::<usize>() {
        Ok(num_t) => num_t,
        Err(_) => panic!("Number of threads must be positive integer"),
    };

    let hybrid = matches.is_present("hybrid");
    let reassign_short = matches.is_present("reassign_short");
    let do_binning = matches.is_present("do_binning");
    let extend_read_clipping = matches.is_present("extend_read_clipping");
    let short_bam_file = matches.value_of("hybrid").unwrap_or("");
    let list_to_phase: Vec<&str>;
    if let Some(values) = matches.values_of("list_to_phase") {
        list_to_phase = values.collect();
    } else {
        list_to_phase = vec![];
    }

    let block_length = matches.value_of("bam_block_length").unwrap_or("15000");
    let block_length = block_length.parse::<usize>().unwrap();
    //    let use_mec = matches.is_present("use_mec");
    let reference_fasta = matches.value_of("reference_fasta").unwrap_or("");
    let filter_supplementary = true;
    let use_supplementary = matches.is_present("use_supplementary");
    let use_gaps = matches.is_present("use_gaps");

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
        .unwrap_or("glopp_out_dir")
        .to_string();
    let snp_density = matches
        .value_of("snp_density")
        .unwrap_or("0.001")
        .parse::<f64>()
        .unwrap();

    if Path::new(&part_out_dir).exists() {
        panic!("Output directory exists; output directory must not be an existing directory");
    }

    std::fs::create_dir_all(&part_out_dir).unwrap();
    let mut cmd_file =
        File::create(format!("{}/cmd.log", part_out_dir)).expect("Can't create file");
    for arg in env::args() {
        write!(cmd_file, "{} ", arg).unwrap();
    }

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

    let vcf_file = matches.value_of("vcf").unwrap();

    if !bam && !frag {
        panic!("Must input a BAM file.")
    }

    if bam && frag {
        panic!("If using frag as input, BAM file should not be specified")
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_t)
        .build_global()
        .unwrap();

    println!("Preprocessing inputs");
    let start_t = Instant::now();
    let contigs_to_phase;
    if bam {
        contigs_to_phase = file_reader::get_contigs_to_phase(&bam_file)
    } else {
        contigs_to_phase = vec!["frag_contig".to_string()];
    }
    println!("Read BAM header successfully.");

    let mut chrom_seqs = FxHashMap::default();
    let (snp_to_genome_pos_t, _genotype_dict_t, _vcf_ploidy) =
        file_reader::get_genotypes_from_vcf_hts(vcf_file);
    let snp_to_genome_pos_map = snp_to_genome_pos_t;
    let vcf_profile = file_reader::get_vcf_profile(&vcf_file, &contigs_to_phase);
    println!("Read VCF successfully.");
    if reference_fasta != "" {
        chrom_seqs = file_reader::get_fasta_seqs(&reference_fasta);
        println!("Read reference fasta successfully.");
    }
    println!("Finished preprocessing in {:?}", Instant::now() - start_t);

    for contig in contigs_to_phase.iter() {
        if !list_to_phase.contains(&&contig[..]) && !list_to_phase.is_empty() {
            continue;
        } else if !vcf_profile.vcf_pos_allele_map.contains_key(contig.as_str())
            || vcf_profile.vcf_pos_allele_map[contig.as_str()].len() < 100
        {
            log::trace!(
                "Contig {} not present or has < 500 variants. Continuing",
                contig
            );
            continue;
        }

        let start_t = Instant::now();
        println!("-----{}-----", contig);
        println!("Reading inputs for contig {} (BAM/VCF/frags).", contig);
        let mut all_frags;
        if frag {
            let mut all_frags_map = file_reader::get_frags_container(frag_file);
            all_frags = all_frags_map.remove(contig).unwrap();
        } else {
            all_frags = file_reader::get_frags_from_bamvcf_rewrite(
                &vcf_profile,
                bam_file,
                short_bam_file,
                filter_supplementary,
                use_supplementary,
                &chrom_seqs,
                &contig,
                use_gaps,
            );
        }
        if all_frags.len() == 0 {
            println!("Contig {} has no fragments", contig);
            continue;
        }

        println!("Time taken reading inputs {:?}", Instant::now() - start_t);

        let start_t = Instant::now();
        if snp_to_genome_pos_map.contains_key(contig) || bam == false {
            let contig_out_dir = format!("{}/{}", part_out_dir, contig);
            fs::create_dir_all(&contig_out_dir).unwrap();

            let mut snp_to_genome_pos: &Vec<usize> = &Vec::new();

            if bam == true {
                snp_to_genome_pos = snp_to_genome_pos_map.get(contig).unwrap();
            }

            //We need frags sorted by first position to make indexing easier. We want the
            //counter_id to reflect the position in the vector.
            all_frags.sort();
            //            all_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));
            for (i, frag) in all_frags.iter_mut().enumerate() {
                frag.counter_id = i;
            }


            //Get last SNP on the genome covered over all fragments.
            let length_gn = utils_frags::get_length_gn(&all_frags);
            println!("Length of genome is {} SNPs", length_gn);
            let mut epsilon = 0.04;
            if hybrid {
                epsilon = 0.02;
            }

            //Do hybrid error correction
            let mut final_frags;
            let mut short_frags = vec![];
            if hybrid {
                let ff_sf = utils_frags::hybrid_correction(all_frags);
                final_frags = ff_sf.0;
                short_frags = ff_sf.1;
                //                final_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));
                final_frags.sort();
                for (i, frag) in final_frags.iter_mut().enumerate() {
                    frag.counter_id = i;
                }
            } else {
                final_frags = all_frags;
            }

            let avg_read_length = utils_frags::get_avg_length(&final_frags, 0.5);
            println!("Median read length is {} SNPs", avg_read_length);

            println!("Number of fragments {}", final_frags.len());

            match matches.value_of("epsilon") {
                None => {}
                Some(value) => {
                    epsilon = value.parse::<f64>().unwrap();
                }
            };

            println!("Epsilon is {}", epsilon);

            let num_locs_string = matches.value_of("num_iters_ploidy_est").unwrap_or("10");
            let num_locs = num_locs_string.parse::<usize>().unwrap();
            let mut hap_graph = graph_processing::generate_hap_graph(
                length_gn,
                num_locs,
                &final_frags,
                epsilon,
                &snp_to_genome_pos,
                max_number_solns,
                block_length,
                contig_out_dir.to_string(),
                snp_density,
            );
            let flow_up_vec =
                graph_processing::solve_lp_graph(&hap_graph, contig_out_dir.to_string());
            let (all_path_parts, path_parts_snp_endpoints) =
                graph_processing::get_disjoint_paths_rewrite(
                    &mut hap_graph,
                    flow_up_vec,
                    epsilon,
                    contig_out_dir.to_string(),
                    &short_frags,
                    reassign_short,
                    &vcf_profile,
                    contig,
                    block_length,
                    do_binning,
                    &snp_to_genome_pos
                );

            file_reader::write_outputs(
                &all_path_parts,
                &path_parts_snp_endpoints,
                contig_out_dir.to_string(),
                &format!("all"),
                &contig, 
                &snp_to_genome_pos,
                extend_read_clipping,
                epsilon
            );
        }

        println!("Total time taken is {:?}", Instant::now() - start_t);
    }
}
