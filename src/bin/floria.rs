extern crate time;
use clap::{AppSettings, Arg, Command};
use floria::file_reader;
use floria::file_writer;
use floria::graph_processing;
use floria::parse_cmd_line;
use floria::part_block_manip;
use floria::solve_flow;
use floria::utils_frags;
use std::fs;
use std::path::Path;
use std::time::Instant;

//This makes statically compiled musl library
//much much faster. Set to default for x86 systems...
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[allow(deprecated)]
fn main() {
    let input_options = "INPUT";
    let output_options = "OUTPUT";
    let alg_options = "ALGORITHM";
    let mandatory_options = "REQUIRED";
    let matches = Command::new("floria")
                          .version("0.0.1")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("floria - polyploid phasing from read sequencing.\n\nExample usage :\nfloria -b bamfile.bam -c vcffile.vcf -o results \n")
                          .arg(Arg::new("bam")
                              .short('b')
                              .value_name("BAM FILE")
                              .required(true)
                              .help("Indexed and sorted bam file to phase.")
                              .takes_value(true)
                              .help_heading(mandatory_options)
                              .display_order(1))
                          .arg(Arg::new("vcf")
                              .short('v')
                              .value_name("VCF FILE")
                              .required(true)
                              .help("VCF file with contig header information; see README if contig header information is not present.")
                              .takes_value(true)
                              .help_heading(mandatory_options))
                          .arg(Arg::new("reference_fasta")
                              .short('r')
                              .takes_value(true)
                              .value_name("FASTA FILE")
                              .help("Reference fasta for the BAM file.")
                              .help_heading(mandatory_options)
                              .required(true))
                          .arg(Arg::new("threads")
                              .long("threads")
                              .short('t')
                              .help("Number of threads to use. (default: 10)")
                              .value_name("INT")
                              .takes_value(true)
                              )
                          .arg(Arg::new("snp_count_filter")
                              .long("snp-count-filter") 
                              .takes_value(true)
                              .help("Skip contigs with less than --snp-count-filter SNPs. (default: 100)")
                              .help_heading(input_options))
                          .arg(Arg::new("output dir")
                              .short('o')
                              .long("output-dir")
                              .help("Output folder. (default: floria_out_dir)")
                              .value_name("STRING")
                              .takes_value(true)
                              .help_heading(output_options))
                          .arg(Arg::new("epsilon")
                              .short('e')
                              .long("epsilon")
                              .takes_value(true)
                              .value_name("FLOAT")
                              .help("Estimated allele call error rate; we recommend: noisy long-reads, 0.02-0.06; hifi, 0.01; short-reads, 0.01. (default: if no -e provided, -e is estimated from data.)")
                              .help_heading(alg_options)
                              .display_order(1))
                          .arg(Arg::new("max_number_solns")
                              .short('n')
                              .long("beam-solns")
                              .takes_value(true)
                              .value_name("INT")
                              .help("Maximum number of solutions for beam search. (default: 10)")
                              .help_heading(alg_options))
                          .arg(Arg::new("max_ploidy")
                              .short('p')
                              .long("max-ploidy")
                              .takes_value(true)
                              .value_name("INT")
                              .help("Maximum ploidy to try to phase up to. (default: 5)")
                              .help_heading(alg_options))
                          .arg(Arg::new("bam_block_length")
                              .short('l')
                              .long("block-length")
                              .takes_value(true)
                              .value_name("INT")
                              .help("Length of blocks (in bp) for flow graph construction when using bam file. (default: 66th percentile of read-lengths, minimum 500 bp)")
                              .help_heading(alg_options)
                              .display_order(1))
                        .arg(Arg::new("debug")
                              .long("debug")
                              .help("Debugging output."))
                        .arg(Arg::new("trace")
                              .long("trace")
                              .help("Trace output."))

                          .arg(Arg::new("snp_density")
                              .short('d')
                              .long("snp-density")
                              .takes_value(true)
                              .value_name("FLOAT")
                              .help("Blocks with fraction of SNPs per base less than -d won't be phased. (default: 0.0005)")
                              .help_heading(alg_options))
                          //Always use qual scores right now.. hidden option
                          .arg(Arg::new("use_qual_scores")
                              .short('q')
                              .takes_value(false)
                              .hidden(true)
                              .help("Use quality scores for reads in MEC optimization. (default: use quality scores)")
                              .help_heading(alg_options))
                            .arg(Arg::new("no stop heuristic")
                              .long("no-stop-heuristic")
                              .takes_value(false)
                              .help("Do not use MEC stopping heuristic for local ploidy computation; only stop phasing when SNP error rate in block is < epsilon. (default: use stopping heuristic)")
                              .help_heading(alg_options))
                           .arg(Arg::new("ploidy sensitivity")
                              .short('s')
                              .long("ploidy-sensitivity")
                              .takes_value(true)
                              .value_name("1,2,3")
                              .help("Sensitivity for the local ploidy stopping heuristic; higher values try to phase more haplotypes, but may give more spurious results. (default: 2)")
                              .help_heading(alg_options))
                          .arg(Arg::new("ignore monomorphic")
                              .long("ignore-monomorphic") 
                              .help("Ignore SNPs that have minor allele frequency less than -e.")
                              .help_heading(input_options))
                          .arg(Arg::new("dont use supplementary")
                              .short('X')
                              .long("no-supp")
                              .help("Do not use supplementary alignments. (default: use supp. alignments with MAPQ = 60)")
                              .help_heading(input_options))
                          .arg(Arg::new("hybrid")
                              .short('H')
                              .long("hybrid")
                              .takes_value(true)
                              .value_name("BAM FILE")
                              .help("UNSTABLE BETA MODE: use short aligned short reads from the same sample to polish long-read SNPs.")
                              .help_heading(input_options))
                          .arg(Arg::new("overwrite")
                              .long("overwrite")
                              .help("Force overwrite for output directory.")
                              .help_heading(output_options))
                            .arg(Arg::new("gzip-reads")
                              .long("gzip-reads")
                              .help("Output gzipped reads instead of raw fastq if using --output-reads.")
                              .help_heading(output_options))
                            .arg(Arg::new("output reads")
                              .long("output-reads")
                              .help("Output reads in fastq format for the resulting haplosets.")
                              .help_heading(output_options))
                          .arg(Arg::new("reassign_short")
                              .long("reassign-short")
                              .help("Reassign short-reads used for polishing to the best haploset. (only when using the -H option)")
                              .help_heading(output_options))
                          .arg(Arg::new("do_binning")
                              .long("bin-by-cov")
                              .help("Increase contiguity by binning haplogroups using coverage (testing in progress).")
                              .help_heading(alg_options)
                              .hidden(true)
                              .display_order(2))
                          .arg(Arg::new("trim")
                              .long("extra-trimming")
                              .help("Trim the reads extra carefully against the reference genome. May cause more fragmented but accurate assemblies.")
                              .help_heading(output_options))
                          .arg(Arg::new("list_to_phase")
                              .short('G')
                              .long("contigs")
                              .multiple(true)
                              .value_name("LIST OF CONTIGS")
                              .takes_value(true)
                              .help("Phase only contigs in this argument. Usage: -G contig1 contig2 contig3 ...")
                              .help_heading(input_options))
                            .arg(Arg::new("mapq_cutoff")
                              .short('m')
                              .long("mapq-cutoff")
                              .takes_value(true)
                              .value_name("INT")
                              .help("Primary MAPQ cutoff. (default: 15)")
                              .help_heading(alg_options))
                            .arg(Arg::new("supp_aln_dist_cutoff")
                              .long("supp-aln-dist-cutoff")
                              .takes_value(true)
                              .value_name("INT")
                              .help("Maximum allowed distance between supp. alignments. (default: 40000)")
                              .help_heading(alg_options))

                          .get_matches();

    let options = parse_cmd_line::parse_cmd_line(matches);

    let start_t_initial = Instant::now();
    log::info!("Preprocessing VCF/Reference");
    let start_t = Instant::now();
    let contigs_to_phase;
    contigs_to_phase = file_reader::get_contigs_to_phase(&options.bam_file);
    let (mut main_bam, mut short_bam) = file_reader::get_bam_readers(&options);
    log::debug!("Read BAM file successfully.");

    let mut chrom_seqs = None;
    //I have two instances of data structures that go from snp positions to genome positions, for
    //some reason. This is bad. snp_to_genome_pos_map and also another map in vcf_profile. TODO fix
    //this. Also, snp_to_genome_pos_map requires a -1 to get the index correct...
    //
    //
    //let snp_to_genome_pos_map = file_reader::get_genotypes_from_vcf_hts(options.vcf_file.clone());
    let vcf_profile = file_reader::get_vcf_profile(&options.vcf_file, &contigs_to_phase);
    let snp_to_genome_pos_map = file_reader::get_genotypes_from_vcf_hts(options.vcf_file.clone());
    log::debug!("Read VCF successfully.");
    if options.reference_fasta != "" {
        chrom_seqs = Some(file_reader::get_fasta_seqs(&options.reference_fasta));
        log::debug!("Read reference fasta successfully.");
    }
    log::info!("Finished preprocessing in {:?}", Instant::now() - start_t);

    let mut warn_first_length = true;
    for contig in contigs_to_phase.iter() {
        if !options.list_to_phase.contains(&contig.to_string()) && !options.list_to_phase.is_empty()
        {
            continue;
        } else if !vcf_profile.vcf_pos_allele_map.contains_key(contig.as_str())
            || vcf_profile.vcf_pos_allele_map[contig.as_str()].len() < options.snp_count_filter
        {
            if warn_first_length {
                log::warn!(
                    "A contig ({}) is not present or has < {} variants. This warning will not be shown from now on. Make sure to change --snp-count-filter if you want to phase small contigs.",
                    contig,
                    options.snp_count_filter,
                );
            }
            warn_first_length = false;
            continue;
        }

        let start_t = Instant::now();
        //log::info!("-----{}-----", contig);
        log::info!(
            "Reading and realigning inputs for contig {} (BAM/VCF).",
            contig
        );
        let mut all_frags;
        let frags_without_snps;
        (all_frags, frags_without_snps) = file_reader::get_frags_from_bamvcf_rewrite(
            &mut main_bam,
            &mut short_bam,
            &vcf_profile,
            &options,
            &mut chrom_seqs,
            &contig,
        );
        log::info!("Number of reads passing filtering: {}", all_frags.len());
        if all_frags.len() == 0 {
            log::debug!("Contig {} has no fragments", contig);
            continue;
        }

        if snp_to_genome_pos_map.contains_key(contig) {
            let contig_out_dir = format!("{}/{}", options.out_dir, contig);

            if Path::new(&contig_out_dir).exists() && options.overwrite {
                let res = fs::remove_dir_all(&contig_out_dir);
                if res.is_err() {
                    log::warn!(
                        "Could not remove {} successfully. Proceeding ...",
                        &contig_out_dir
                    );
                }
            }

            fs::create_dir_all(&contig_out_dir).unwrap();

            let snp_to_genome_pos: &Vec<usize>;
            snp_to_genome_pos = snp_to_genome_pos_map.get(contig).unwrap();

            //We need frags sorted by first position to make indexing easier. We want the
            //counter_id to reflect the position in the vector.
            all_frags.sort();
            //            all_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));
            for (i, frag) in all_frags.iter_mut().enumerate() {
                frag.counter_id = i;
            }

            //Get last SNP on the genome covered over all fragments.
            let length_gn = utils_frags::get_length_gn(&all_frags);
            log::info!("Contig {} has {} SNPs", contig, length_gn);

            //Do hybrid error correction
            let mut final_frags;
            let mut short_frags = vec![];
            if options.hybrid {
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

            if options.ignore_monomorphic {
                final_frags = utils_frags::remove_monomorphic_allele(final_frags, options.epsilon);
            }

            log::info!(
                "Reading inputs, realigning time taken {:?}",
                Instant::now() - start_t
            );

            let avg_read_length = utils_frags::get_avg_length(&final_frags, 0.5);
            log::debug!("Median number of SNPs in a read is {}", avg_read_length);
            log::debug!("Number of fragments {}", final_frags.len());
            log::debug!("Epsilon is {}", options.epsilon);

            log::info!("Local phasing with {} threads...", options.num_threads);
            let phasing_t = Instant::now();
            let mut hap_graph = graph_processing::generate_hap_graph(
                &final_frags,
                &snp_to_genome_pos,
                contig_out_dir.to_string(),
                &options,
            );
            log::info!("Phasing time taken {:?}", Instant::now() - phasing_t);

            log::info!("Solving flow problem...");
            let highs_t = Instant::now();
            let flow_up_vec = solve_flow::solve_lp_graph(&hap_graph);
            log::info!("Flow solved in time {:?}", Instant::now() - highs_t);

            //            let minilp_t = Instant::now();
            //            let flow_up_vec =
            //                solve_flow::solve_lp_graph_minilp(&hap_graph, contig_out_dir.to_string());
            //            log::debug!("minilp time taken {:?}", Instant::now() - minilp_t);

            let (all_path_parts, path_parts_snp_endpoints) =
                graph_processing::get_disjoint_paths_rewrite(
                    &mut hap_graph,
                    flow_up_vec,
                    contig_out_dir.to_string(),
                    &vcf_profile,
                    contig,
                    &options,
                );

            let (sorted_path_parts, sorted_snp_endpoints) =
                part_block_manip::process_reads_for_final_parts(
                    all_path_parts,
                    &short_frags,
                    path_parts_snp_endpoints,
                    &options,
                    &snp_to_genome_pos,
                );

            //Scan over path_parts_snp_endpoints and only put in frags_without_snps that end there.
            let snpless_frags_between_gaps = part_block_manip::get_frags_in_snpless_gaps(
                &sorted_snp_endpoints,
                &snp_to_genome_pos,
                &frags_without_snps,
                options.block_length,
                &final_frags,
            );

            file_writer::write_outputs(
                &sorted_path_parts,
                &sorted_snp_endpoints,
                contig_out_dir.to_string(),
                &format!("{}", &contig),
                &contig,
                &snp_to_genome_pos,
                &options,
                &snpless_frags_between_gaps,
                &mut chrom_seqs.as_mut().unwrap()
            );
        }
    }
    log::info!("Total time taken is {:?}", Instant::now() - start_t_initial);
}
