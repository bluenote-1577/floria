use crate::types_structs::Frag;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use fxhash::FxHashMap;
use crate::types_structs::{Genotype, SnpPosition, GnPosition};
//Only do realign around SNPs, will add around deletions later
pub fn realign(
    ref_gn: &[u8],
    frag: &mut Frag,
    var_to_gn_pos: &FxHashMap<SnpPosition, GnPosition>,
    gn_pos_to_allele: &FxHashMap<GnPosition, Vec<Genotype>>,
) {
    let flank = 16;
    let block_size = 8;
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };
    let mut a = Block::<false, false>::new(2 * flank, 2 * flank , 2 * flank);
    for (snp_pos, orig_geno) in frag.seq_dict.iter_mut() {
        let snp_gn_pos = var_to_gn_pos[snp_pos] as usize;
        let snp_q_pos = frag.snp_pos_to_seq_pos[&(snp_pos)].1 as usize;
        if !(flank > snp_gn_pos
            || flank + snp_gn_pos >= ref_gn.len()
            || flank > snp_q_pos
            || flank + snp_q_pos >= frag.seq_string[0].len())
        {
            let mut ref_str = ref_gn[snp_gn_pos - flank..snp_gn_pos + flank].to_vec();
            let alleles = &gn_pos_to_allele[&snp_gn_pos];
            let mut best_score = i32::MIN;
            let mut best_geno = 0;
            let q = PaddedBytes::from_bytes::<NucMatrix>(
                &frag.seq_string[0].slice(snp_q_pos - flank, snp_q_pos + flank).ascii(),
                block_size,
            );
            for i in 0..alleles.len() {
                ref_str[flank] = alleles[i] as u8;
                if ref_str.iter().any(|x| x.to_ascii_uppercase() < b'A' || x.to_ascii_uppercase() > b'Z'){
                    dbg!(&ref_str);
                    panic!();
                }
                let r = PaddedBytes::from_bytes::<NucMatrix>(&ref_str, block_size);

                //let a = Block::<false, false>::align(&q, &r, &NW1, gaps, block_size..=block_size, 0);

                // Align with traceback, but no x drop threshold.
//                let a =
//                    Block::<false,false>::align(&q, &r, &NW1, gaps, block_size..=block_size, 0);
                a.align(&q, &r, &NW1, gaps, block_size..=block_size, 0);
                let res = a.res();
                let score = res.score;
                if score > best_score {
                    best_score = score;
                    best_geno = i as Genotype;
                }

            }
            if *orig_geno != best_geno {
//                println!("Called geno {}, best geno realign {} at {}", orig_geno, best_geno, snp_pos);
            }
            *orig_geno = best_geno;
        }
    }
}
