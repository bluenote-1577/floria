use crate::types_structs::Frag;
use block_aligner::cigar::Operation;
use block_aligner::cigar::*;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use fxhash::FxHashMap;

//Only do realign around SNPs, will add around deletions later
pub fn realign(
    ref_gn: &[u8],
    frag: &mut Frag,
    var_to_gn_pos: &FxHashMap<i64, i64>,
    gn_pos_to_allele: &FxHashMap<i64, Vec<u8>>,
) {
    let flank = 16;
    let block_size = 8;
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };
    for (snp_pos, orig_geno) in frag.seq_dict.iter_mut() {
        let snp_gn_pos = var_to_gn_pos[&(*snp_pos as i64)] as usize;
        let snp_q_pos = frag.snp_pos_to_seq_pos[&(snp_pos)].1 as usize;
        if !(flank > snp_gn_pos
            || flank + snp_gn_pos >= ref_gn.len()
            || flank > snp_q_pos
            || flank + snp_q_pos >= frag.seq_string[0].len())
        {
            let mut ref_str = ref_gn[snp_gn_pos - flank..snp_gn_pos + flank].to_vec();
            let alleles = &gn_pos_to_allele[&(snp_gn_pos as i64)];
            let mut best_score = i32::MIN;
            let mut best_geno = 0;
            let q = PaddedBytes::from_bytes::<NucMatrix>(
                &frag.seq_string[0].to_ascii_vec()[snp_q_pos - flank..snp_q_pos + flank],
                block_size,
            );
            for i in 0..alleles.len() {
                ref_str[flank] = alleles[i];
                let r = PaddedBytes::from_bytes::<NucMatrix>(&ref_str, block_size);

                // Align with traceback, but no x drop threshold.
                let a =
                    Block::<_, true, false>::align(&q, &r, &NW1, gaps, block_size..=block_size, 0);
                let res = a.res();
                let score = res.score;
                if score > best_score {
                    best_score = score;
                    best_geno = i;
                }

            }
            if *orig_geno != best_geno {
//                println!("Called geno {}, best geno realign {} at {}", orig_geno, best_geno, snp_pos);
            }
            *orig_geno = best_geno;
        }
    }
}
