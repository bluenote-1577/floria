extern crate time;
use clap::{App, AppSettings, Arg};
use glopp::file_reader;
use glopp::file_writer;
use glopp::types_structs::*;
use std::time::Instant;

fn main() {
    let matches = App::new("vartig-dump")
                          .version("0.1.0")
                          .setting(AppSettings::ArgRequiredElseHelp)
                          .about("Turn VCF + BAM -> Fragment files. Output can be used to debug or to input into other haplotype phasing algorithms.\n\nExample usage : frag-dump -b bamfile.bam -v vcffile.vcf -o output.txt") 
                          .arg(Arg::with_name("bam")
                              .short('b')
                              .value_name("BAMFILE")
                               .help("Input a bam file.")
                                .takes_value(true)
                                .required(true))
                          .arg(Arg::with_name("vcf")
                               .short('v')
                               .help("Input a VCF : Mandatory if using BAM file; Enables genotype polishing if using frag file.")
                               .value_name("VCFFILE")
                               .takes_value(true)
                               .required(true))
                          .arg(Arg::with_name("output")
                              .short('o')
                              .help("Name of output file (default : BAMFILE_vartigs.txt)")
                              .value_name("OUTPUT")
                              .takes_value(true)
                              .required(true))
                          .get_matches();

    let mut options = Options::default();
    options.supp_aln_dist_cutoff = 10000000000;
    let bam_file = matches.value_of("bam").unwrap();
    options.bam_file = bam_file.to_string();
    options.mapq_cutoff = 30;
    let contigs_to_phase = file_reader::get_contigs_to_phase(&bam_file);

    let vcf_file = matches.value_of("vcf").unwrap();
    let snp_to_genome_pos_t =
        file_reader::get_genotypes_from_vcf_hts(vcf_file.clone());
    let snp_to_genome_pos_map = snp_to_genome_pos_t;
    let vcf_profile = file_reader::get_vcf_profile(&vcf_file, &contigs_to_phase);

    let vtig_string = format!("{}_vartigs.txt", bam_file);
    let output_frag_str = matches.value_of("output").unwrap_or(&vtig_string);

    //CONSTANTS - Constants which users probably should not change.

    let (mut main_bam, mut short_bam) = file_reader::get_bam_readers(&options);
    let mut chrom_seqs = None;
    for contig in contigs_to_phase.iter(){
        let mut all_frags = file_reader::get_frags_from_bamvcf_rewrite(&mut main_bam, &mut short_bam, &vcf_profile, &options, &mut chrom_seqs, &contig );
        all_frags.sort_by(|a, b| a.first_position.cmp(&b.first_position));
        file_writer::write_alignment_as_vartig(&all_frags, output_frag_str, contig, &snp_to_genome_pos_map[contig], 1, snp_to_genome_pos_map[contig].len() as SnpPosition, output_frag_str);
   }
}
