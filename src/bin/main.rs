use haplotype_phaser::file_reader;
use haplotype_phaser::local_clustering;
use haplotype_phaser::utils_frags;
use haplotype_phaser::vcf_polishing;
use std::env;

fn main() {
    let ploidy = 4;

    let args: Vec<String> = env::args().collect();
    println!("Reading frags");
    let all_frags = file_reader::get_frags_container(&args[1]);
    let genotype_dict = file_reader::get_genotypes_from_vcf(&args[2]);

    println!("Indexing reads");
    let pos_read_sets = utils_frags::get_all_overlaps(&all_frags);

    println!("Processing distances");
    let all_distances = utils_frags::get_all_distances(&all_frags);

    println!("Generating haplotype blocks");
    let part = local_clustering::generate_hap_block(1,100,&pos_read_sets,&all_distances,ploidy);
    let hap_block = utils_frags::hap_block_from_partition(&part);
    vcf_polishing::polish_using_vcf(&genotype_dict,&hap_block,vec![1,5,10,15]);


    //for map in hap_block.blocks.iter(){
    //    let mut v: Vec<_> = map.into_iter().collect();
    //    v.sort_by(|x,y| x.0.cmp(&y.0));
    //    for thing in v.iter(){
    //        println!("{:?}",thing);
    //    }
    //}
    

//    for read_set in pos_read_sets.iter(){
//        println!("{:?}",read_set)
//    }
}


