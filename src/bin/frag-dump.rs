extern crate time;
use clap::{App, AppSettings, Arg};
use flopp::file_reader;
use std::time::Instant;

fn main() {
    let matches = App::new("frag-dump")
                          .version("0.1.0")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("Turn VCF + BAM -> Fragment files. Output can be used to debug or to input into other haplotype phasing algorithms.\n\nExample usage : frag-dump -b bamfile.bam -v vcffile.vcf -o output.txt") 
                          .arg(Arg::with_name("bam")
                              .short("b")
                              .value_name("BAMFILE")
                               .help("Input a bam file.")
                                .takes_value(true)
                                .required(true))
                          .arg(Arg::with_name("vcf")
                               .short("v")
                               .help("Input a VCF : Mandatory if using BAM file; Enables genotype polishing if using frag file.")
                               .value_name("VCFFILE")
                               .takes_value(true)
                               .required(true))
                          .arg(Arg::with_name("output")
                              .short("o")
                              .help("Name of output file (default : flopp_out.txt)")
                              .value_name("OUTPUT")
                              .takes_value(true)
                              .required(true))
                          .get_matches();

    let bam_file = matches.value_of("bam").unwrap();

    //Whether or not we polish using genotyping information from VCF.
    let vcf_file = matches.value_of("vcf").unwrap();
    let start_t = Instant::now();
    let output_frag_str = matches.value_of("output").unwrap_or("flopp_frags.txt");

    //CONSTANTS - Constants which users probably should not change.

    println!("Reading frags");
    let mut all_frags = file_reader::get_frags_from_bamvcf(vcf_file, bam_file);

    //We need frags sorted by first position to make indexing easier.
    all_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));
    file_reader::write_frags_file(all_frags, output_frag_str.to_string());
    println!(
        "Time taken reading and writing fragments to {} : {:?}",
        output_frag_str,
        Instant::now() - start_t
    );
}
