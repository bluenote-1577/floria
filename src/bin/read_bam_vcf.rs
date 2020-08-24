use haplotype_phaser::file_reader;
use std::env;

fn main(){
    let args: Vec<String> = env::args().collect();
    let vcf_file = &args[1];
    let bam_file = &args[2];

    let frags = file_reader::get_frags_from_bamvcf(vcf_file,bam_file);
    dbg!(frags.len());
//    for frag in frags{
//        dbg!(frag.first_position,frag.last_position,frag.id,frag.seq_dict.len());
//    }
}
