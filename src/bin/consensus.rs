extern crate time;
use clap::{App, AppSettings, Arg};
use flopp::file_reader;
use flopp::local_clustering;
use flopp::types_structs::Frag;
use flopp::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
fn main() {
    let matches = App::new("glopp")
        .version("0.1.0")
        .setting(AppSettings::ArgRequiredElseHelp)
        .about("Haplotype consensus module. ")
        .arg(
            Arg::with_name("vcf")
                .short("c")
                .value_name("VCFFILE")
                .help("Input a VCF.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("bam")
                .short("b")
                .multiple(true)
                .value_name("BAMFILE1 BAMFILE2 ...")
                .help("Input bam files.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .help("Name of output file (default : consensus_out.txt)")
                .value_name("OUTPUT")
                .takes_value(true),
        )
        .get_matches();

    //If the user is getting frag files from BAM and VCF.
    let bam;
    let bam_files: Vec<&str>;
    let _bm = match matches.value_of("bam") {
        None => {
            bam = false;
            "_"
        }
        Some(bam_files) => {
            bam = true;
            bam_files
        }
    };

    if !bam {
        panic!("No bam file input found");
    }

    bam_files = matches.values_of("bam").unwrap().collect();
    dbg!(&bam_files);

    //Whether or not we polish using genotyping information from VCF.
    let vcf;
    let vcf_file = match matches.value_of("vcf") {
        None => {
            vcf = false;
            "_"
        }
        Some(vcf_file) => {
            vcf = true;
            vcf_file
        }
    };

    if !vcf {
        panic!("No VCF file input found");
    }

    let part_out_dir = matches.value_of("output").unwrap_or("glopp_out_dir").to_string();
    let mut snp_to_genome_pos_map: FxHashMap<String, Vec<usize>> = FxHashMap::default();


    if vcf{
        let (snp_to_genome_pos_t, _genotype_dict_t, vcf_ploidy) =
            file_reader::get_genotypes_from_vcf_hts(vcf_file);
        snp_to_genome_pos_map = snp_to_genome_pos_t;

        //If the VCF file is misformatted or has weird genotyping call we can catch that here.
    }



    let mut final_part_owned = vec![];
    let mut final_part_reference = vec![];
    let mut length_gn = 0;
    let mut contig_n = String::from("");
    for bam_file in bam_files {
        let all_frags_map = file_reader::get_frags_from_bamvcf(vcf_file, bam_file, true, true);
        for (contig, bam_fragments) in all_frags_map.iter() {
            let length_gn_bam = utils_frags::get_length_gn(&bam_fragments);
            if length_gn_bam > length_gn {
                length_gn = length_gn_bam;
            }
            contig_n = contig.clone();
            let as_map: FxHashSet<Frag> = bam_fragments.iter().cloned().collect();
            final_part_owned.push(as_map);
        }
    }
    let snp_to_genome_pos = snp_to_genome_pos_map.get(&contig_n).unwrap();


    println!("Max length of genome is {} SNPS.", length_gn);


    for part in final_part_owned.iter() {
        let set_ref: FxHashSet<&Frag> = part.iter().collect();
        final_part_reference.push(set_ref);
    }

    let first_iter = true;
    let final_block_unpolish = utils_frags::hap_block_from_partition(&final_part_reference);
    let (f_binom_vec, f_freq_vec) =
        local_clustering::get_partition_stats(&final_part_reference, &final_block_unpolish);

    let (f_binom_vec_rf, f_freq_vec_rf) =
        local_clustering::get_partition_stats_ref_wild(&final_part_reference, &final_block_unpolish);


    let final_score = local_clustering::get_mec_score(&f_binom_vec, &f_freq_vec, 0.0, 0.0);
    let mut total_num_alleles = 0;
    for (good,bad) in f_binom_vec{
        total_num_alleles += good;
        total_num_alleles += bad;
    }

    println!("Expected MEC error with error rate 0.05 if only 1 haplotype present is {}",total_num_alleles as f64 * 0.05);
    println!("Good/Bad for Ref/Wild: {:?}", f_binom_vec_rf);
    println!("Error rate for Ref: {}", f_binom_vec_rf[2].0.1 as f64/ (f_binom_vec_rf[2].0.1 as f64 + f_binom_vec_rf[2].0.0 as f64));
    println!("Error rate for Wild: {}", f_binom_vec_rf[2].1.1 as f64/ (f_binom_vec_rf[2].1.1 as f64 + f_binom_vec_rf[2].1.0 as f64));

    println!(
        "Final MEC score for the partition is {:?}.",
        -1.0 * final_score
    );

    file_reader::write_output_partition_to_file(
        &final_part_reference,
        vec![],
        part_out_dir.clone(),
        &String::from("cons"),
    );


    file_reader::write_blocks_to_file(
        part_out_dir.clone(),
        &vec![final_block_unpolish],
        &vec![length_gn],
        &snp_to_genome_pos,
        &final_part_reference,
        first_iter,
        &String::from("cons"),
        &FxHashMap::default()
    );
}
