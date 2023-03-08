extern crate time;
use clap::{AppSettings, Arg, Command};
use fxhash::FxHashMap;
use glopp::file_reader;
use glopp::solve_flow;
use glopp::parse_cmd_line;
use glopp::graph_processing;
use glopp::part_block_manip;
use glopp::utils_frags;
use std::fs;
use std::time::Instant;

//This makes statically compiled musl library
//much much faster. Set to default for x86 systems...
#[cfg(target_env="musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[allow(deprecated)]
fn main() {
    let input_options = "INPUT";
    let output_options = "OUTPUT";
    let alg_options = "ALGORITHM";
    let mandatory_options = "REQUIRED";
    let matches = Command::new("glopp")
                          .version("0.0.1")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("glopp - polyploid phasing from read sequencing.\n\nExample usage :\nglopp -b bamfile.bam -c vcffile.vcf -o results \n")
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
                          .arg(Arg::new("threads")
                              .short('t')
                              .help("Number of threads to use. (default: 10).")
                              .value_name("INT")
                              .takes_value(true)
                              )
                          .arg(Arg::new("output dir")
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
                          .arg(Arg::new("max_ploidy")
                              .short('p')
                              .takes_value(true)
                              .value_name("INT")
                              .help("Maximum ploidy to try to phase up to. (default: 5)")
                              .help_heading(alg_options))
                          .arg(Arg::new("bam_block_length")
                              .short('l')
                              .takes_value(true)
                              .value_name("INT")
                              .help("Length of blocks (in nucleotides) for flow graph construction when using bam file. (default: 15000)")
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
                              .takes_value(true)
                              .value_name("FLOAT")
                              .help("Minimum SNP density to phase. Blocks with SNP density less than this value will not be phased. (Default: 0.001 i.e. 1 SNP per 1000 bases)")
                              .help_heading(alg_options))
                          //Always use qual scores right now.. hidden option
                          .arg(Arg::new("use_qual_scores")
                              .short('q')
                              .takes_value(false)
                              .hidden(true)
                              .help("Use quality scores for reads in MEC optimization. (Default: do not use quality scores)")
                              .help_heading(alg_options))
                            .arg(Arg::new("no stop heuristic")
                              .long("no-stop-heuristic")
                              .takes_value(false)
                              .help("Do not use stopping heuristic for local ploidy computation, only epsilon stopping criterion. (Default: use stopping heuristic)")
                              .help_heading(alg_options))
                          .arg(Arg::new("snp_count_filter")
                              .long("snp-count-filter") 
                              .takes_value(true)
                              .help("Skip contigs with less than --snp-count-filter SNPs (Default: 100)")
                              .help_heading(input_options))
                          .arg(Arg::new("use monomorphic")
                              .long("use-monomorphic") 
                              .help("Use SNPs that have minor allele frequency less than epsilon (Default: ignore monomorphic SNPs)")
                              .help_heading(input_options))
                          .arg(Arg::new("use_supplementary")
                              .short('X')
                              .help("Use supplementary alignments (Default: don't use).")
                              .help_heading(input_options))
                          .arg(Arg::new("hybrid")
                              .short('H')
                              .takes_value(true)
                              .value_name("BAM FILE")
                              .help("RECOMMENDED: Use short aligned short reads to polish long-read SNPs.")
                              .help_heading(input_options)
                              .display_order(1))
                            .arg(Arg::new("gzip-reads")
                              .long("gzip-reads")
                              .help("output gzipped reads. ")
                              .help_heading(output_options))
                            .arg(Arg::new("no output reads")
                              .long("no-output-reads")
                              .help("do not output reads in fastq format. (default: reads are output)")
                              .help_heading(output_options))
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
                          .arg(Arg::new("list_to_phase")
                              .short('G')
                              .multiple(true)
                              .value_name("LIST OF CONTIGS")
                              .takes_value(true)
                              .help("Phase only contigs in this argument. Usage: -G contig1 contig2 contig3 ...")
                              .help_heading(input_options))
                            .arg(Arg::new("mapq_cutoff")
                              .short('m')
                              .takes_value(true)
                              .value_name("INT")
                              .help("MAPQ cutoff. (default: 15)")
                              .help_heading(alg_options))
                          .get_matches();

    let options = parse_cmd_line::parse_cmd_line(matches);

    let start_t_initial = Instant::now();
    log::info!("Preprocessing VCF/Reference");
    let start_t = Instant::now();
    let contigs_to_phase;
    contigs_to_phase = file_reader::get_contigs_to_phase(&options.bam_file);
    log::debug!("Read BAM header successfully.");

    let mut chrom_seqs = FxHashMap::default();
    let (snp_to_genome_pos_t, _genotype_dict_t, _vcf_ploidy) =
        file_reader::get_genotypes_from_vcf_hts(options.vcf_file.clone());
    let snp_to_genome_pos_map = snp_to_genome_pos_t;
    let vcf_profile = file_reader::get_vcf_profile(&options.vcf_file, &contigs_to_phase);
    log::debug!("Read VCF successfully.");
    if options.reference_fasta != "" {
        chrom_seqs = file_reader::get_fasta_seqs(&options.reference_fasta);
        log::info!("Read reference fasta successfully.");
    }
    log::info!("Finished preprocessing in {:?}", Instant::now() - start_t);

    for contig in contigs_to_phase.iter() {
        if !options.list_to_phase.contains(&contig.to_string()) && !options.list_to_phase.is_empty() {
            continue;
        } else if !vcf_profile.vcf_pos_allele_map.contains_key(contig.as_str())
            || vcf_profile.vcf_pos_allele_map[contig.as_str()].len() < options.snp_count_filter
        {
            log::warn!(
                "Contig '{}' not present or has < {} variants. Continuing (change --snp-count-filter to phase small contigs)",
                contig,
                options.snp_count_filter,
            );
            continue;
        }

        let start_t = Instant::now();
        //log::info!("-----{}-----", contig);
        log::info!("Reading inputs for contig {} (BAM/VCF).", contig);
        let mut all_frags;
        all_frags = file_reader::get_frags_from_bamvcf_rewrite(
            &vcf_profile,
            &options,
            &chrom_seqs,
            &contig,
        );
        if all_frags.len() == 0 {
            log::warn!("Contig {} has no fragments", contig);
            continue;
        }

        if snp_to_genome_pos_map.contains_key(contig){
            let contig_out_dir = format!("{}/{}", options.out_dir, contig);
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

            if !options.use_monomorphic{
                final_frags = utils_frags::remove_monomorphic_allele(final_frags, options.epsilon);
            }

            log::info!("Time taken reading inputs {:?}", Instant::now() - start_t);

            let avg_read_length = utils_frags::get_avg_length(&final_frags, 0.5);
            log::debug!("Median number of SNPs in a read is {}", avg_read_length);
            log::debug!("Number of fragments {}", final_frags.len());
            log::debug!("Epsilon is {}", options.epsilon);

            log::info!("Local phasing with {} threads...", options.num_threads);
            let phasing_t= Instant::now();
            let mut hap_graph = graph_processing::generate_hap_graph(
                &final_frags,
                &snp_to_genome_pos,
                contig_out_dir.to_string(),
                &options
            );
            log::info!("Phasing time taken {:?}", Instant::now() - phasing_t);

            log::info!("Solving flow problem...");
            let highs_t = Instant::now();
            let flow_up_vec =
                solve_flow::solve_lp_graph(&hap_graph);
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

            let (sorted_path_parts, sorted_snp_endpoints) = part_block_manip::process_reads_for_final_parts(
                all_path_parts,
                &short_frags,
                path_parts_snp_endpoints,
                &options,
                &snp_to_genome_pos,
            );

            file_reader::write_outputs(
                &sorted_path_parts,
                &sorted_snp_endpoints,
                contig_out_dir.to_string(),
                &format!("all"),
                &contig,
                &snp_to_genome_pos,
                &options,
            );
        }

        log::info!("Total time taken is {:?}", Instant::now() - start_t_initial);
    }
}
