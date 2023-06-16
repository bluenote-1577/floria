use ordered_float::*;
use derivative::Derivative;
use debruijn::dna_string::DnaString;
use crate::utils_frags;
use fxhash::{FxHashMap, FxHashSet};
use rust_htslib::bam::Record;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::rc::Rc;

pub type GnPosition = usize;
pub type SnpPosition = u32;
pub type Genotype = u8;
pub type GenotypeCount = OrderedFloat<f64>;
pub type Haplotype = FxHashMap<SnpPosition, FxHashMap<Genotype, GenotypeCount>>;
pub static GAP_CHAR: Genotype = 9;

pub type FlowUpVec = Vec<((usize, usize), (usize, usize), f64)>;

#[derive(Debug, Clone, Derivative)]
#[derivative(Default)]
pub struct Options{
    pub bam_file: String,
    pub vcf_file: String,
    pub use_qual_scores: bool,
    pub gzip: bool,
    pub output_reads: bool,
    pub mapq_cutoff: u8,
    pub epsilon: f64,
    pub dont_use_supp_aln: bool,
    pub reassign_short: bool,
    pub do_binning: bool,
    pub max_number_solns: usize,
    pub snp_density: f64,
    pub max_ploidy: usize,
    pub out_dir: String,
    pub hybrid: bool,
    pub list_to_phase: Vec<String>,
    pub block_length: usize,
    pub reference_fasta: String,
    pub trim_reads: bool,
    pub short_bam_file: String,
    pub snp_count_filter: usize,
    pub stopping_heuristic: bool,
    pub ignore_monomorphic: bool,
    pub num_threads: usize,
    pub overwrite: bool,
    pub ploidy_sensitivity: u8,
    #[derivative(Default(value="40000"))]
    pub supp_aln_dist_cutoff: i64
}

#[derive(Debug, Clone, Default)]
pub struct VcfProfile<'a> {
    pub vcf_pos_allele_map: FxHashMap<&'a str, FxHashMap<GnPosition, Vec<Genotype>>>,
    pub vcf_pos_to_snp_counter_map: FxHashMap<&'a str, FxHashMap<GnPosition, SnpPosition>>,
    pub vcf_snp_pos_to_gn_pos_map: FxHashMap<&'a str, FxHashMap<SnpPosition, GnPosition>>,
}

#[derive(Debug, Clone)]
pub struct TraceBackNode {
    pub score: f64,
    pub prev_ind: Option<usize>,
    pub is_sink: bool,
    pub is_source: bool,
}
//Positions are inclusive
#[derive(Eq, Debug, Clone, Default)]
pub struct Frag {
    pub id: String,
    pub counter_id: usize,
    pub seq_dict: FxHashMap<SnpPosition, Genotype>,
    pub qual_dict: FxHashMap<SnpPosition, u8>,
    pub first_position: SnpPosition,
    pub last_position: SnpPosition,
    pub positions: FxHashSet<SnpPosition>,
    pub seq_string: Vec<DnaString>,
    pub qual_string: Vec<Vec<u8>>, //Needs to be PHRED scaled i.e. +33 from value
    pub is_paired: bool,
    pub snp_pos_to_seq_pos: FxHashMap<SnpPosition, (u8, GnPosition)>,
    pub first_pos_base: GnPosition,
    pub last_pos_base: GnPosition

    
}

impl Ord for Frag{
    fn cmp(&self, other: &Frag) -> Ordering {
        return (self.first_position,other.last_position,self.counter_id).cmp(&(other.first_position,self.last_position,other.counter_id));
        //I tried the below ordering. Gives similar results. 
        //return (self.first_position,other.seq_dict.len(),self.counter_id).cmp(&(other.first_position,self.seq_dict.len(),other.counter_id));
    }
}

impl PartialOrd for Frag{
    fn partial_cmp(&self, other: &Frag) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Hash for Frag {
    fn hash<H: Hasher>(&self, _state: &mut H) {
        //        self.counter_id;
        self.counter_id.hash(_state);
    }
}

impl PartialEq for Frag {
    fn eq(&self, other: &Self) -> bool {
        self.counter_id == other.counter_id
    }
}

#[derive(PartialEq)]
pub struct SearchNode<'a> {
    pub read: &'a Frag,
    pub part: usize,
    pub score: f64,
    pub freqs: Vec<usize>,
    pub error_vec: Vec<(f64, f64)>,
    pub block_id: usize,
    pub parent_node: Option<Rc<SearchNode<'a>>>,
    pub current_pos: SnpPosition,
    pub broken_blocks: FxHashSet<usize>,
}

impl Ord for SearchNode<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score.partial_cmp(&other.score).unwrap()
    }
}

impl PartialOrd for SearchNode<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for SearchNode<'_> {}

impl Drop for SearchNode<'_> {
    fn drop(&mut self) {
        let mut prev = self.parent_node.take();

        while let Some(rc) = prev {
            if let Ok(mut graph) = Rc::try_unwrap(rc) {
                prev = graph.parent_node.take();
            } else {
                break;
            }
        }
    }
}

pub struct HapNode<'a> {
    pub frag_set: FxHashSet<&'a Frag>,
    pub out_edges: Vec<(usize, f64)>,
    pub in_edges: Vec<(usize, f64)>,
    pub column: usize,
    pub row: usize,
    cov: f64,
    pub id: usize,
    pub out_flows: Vec<(usize, f64)>,
    pub hap_map: Haplotype,
    pub snp_endpoints: (SnpPosition, SnpPosition),
}

impl<'a> HapNode<'a> {
    pub fn new(frag_set: FxHashSet<&'a Frag>, snp_endpoints: (SnpPosition, SnpPosition)) -> HapNode<'a> {
        let mut hap_map = FxHashMap::default();
        for frag in frag_set.iter() {
            for pos in frag.seq_dict.keys(){
                if *pos <= snp_endpoints.1 && *pos >= snp_endpoints.0 {
                    let var_at_pos = frag.seq_dict.get(pos).unwrap();
                    let sites = hap_map.entry(*pos).or_insert(FxHashMap::default());
                    let site_counter = sites.entry(*var_at_pos).or_insert(OrderedFloat(0.));
                    *site_counter += utils_frags::phred_scale(frag, pos);
                }
            }
        }
        let mut allele_cov_list = vec![];
        for values1 in hap_map.values() {
            for allele_count in values1.values() {
                allele_cov_list.push(*allele_count);
            }
        }

        allele_cov_list.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let cov;
        if allele_cov_list.is_empty() {
            cov = 0.;
        } else {
            cov = *allele_cov_list[allele_cov_list.len() * 2 / 3];
        }
        //        let cov = allele_cov_list.last().unwrap_or(&0.);
        let toret = HapNode {
            frag_set: frag_set,
            out_edges: vec![],
            in_edges: vec![],
            column: usize::MAX,
            row: usize::MAX,
            id: usize::MAX,
            cov: cov,
            out_flows: vec![],
            hap_map: hap_map,
            snp_endpoints: snp_endpoints,
        };
        return toret;
    }

    pub fn cov(&self) -> f64 {
        return self.cov;
    }
}

pub fn build_child_node<'a>(
    read: &'a Frag,
    part: usize,
    block_id: usize,
    parent_node_option: Option<Rc<SearchNode<'a>>>,
    error_vec: Vec<(f64, f64)>,
    score: f64,
    current_pos: SnpPosition,
) -> SearchNode<'a> {
    let mut new_freqs = vec![];
    let mut new_parent_node_option: Option<Rc<SearchNode>> = None;
    if let Some(parent_node) = parent_node_option {
        for freq in parent_node.freqs.iter() {
            new_freqs.push(*freq);
        }
        new_parent_node_option = Some(parent_node);
        new_freqs[part] += 1;
    } else {
        new_freqs = vec![1; error_vec.len()];
    }
    let updated_freqs = new_freqs;

    let toret = SearchNode {
        read: read,
        part: part,
        score: score,
        freqs: updated_freqs,
        error_vec: error_vec,
        block_id: block_id,
        parent_node: new_parent_node_option,
        current_pos: current_pos,
        broken_blocks: FxHashSet::default(),
    };

    toret
}

#[derive(Debug, PartialEq, Eq)]
pub struct HapBlock {
    pub blocks: Vec<Haplotype>,
}

impl Ord for HapBlock {
    fn cmp(&self, other: &Self) -> Ordering {
        self.blocks.len().partial_cmp(&other.blocks.len()).unwrap()
    }
}

impl PartialOrd for HapBlock {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn build_frag(id: String, counter_id: usize, is_paired: bool) -> Frag {
    let toret = Frag {
        id: id,
        counter_id: counter_id,
        seq_dict: FxHashMap::default(),
        qual_dict: FxHashMap::default(),
        first_position: SnpPosition::MAX,
        last_position: SnpPosition::MIN,
        positions: FxHashSet::default(),
        seq_string: vec![DnaString::new(); 2],
        qual_string: vec![vec![]; 2],
        is_paired: is_paired,
        snp_pos_to_seq_pos: FxHashMap::default(),
        first_pos_base: GnPosition::MAX,
        last_pos_base: GnPosition::MAX,
    };

    toret
}
#[inline]
pub fn update_frag(
    frag: &mut Frag,
    geno: Genotype,
    snp_pos: SnpPosition,
    qual: u8,
    pair_number: u8,
    is_supp: bool,
    record: &Record,
    qpos: GnPosition,
) {
    frag.seq_dict.insert(snp_pos, geno);
    frag.qual_dict.insert(snp_pos, qual);
    let mut seq_pos = qpos;
    if is_supp {
        let clipping_offset = record.cigar().leading_hardclips();
        seq_pos = qpos + clipping_offset as usize;
    }
    frag.snp_pos_to_seq_pos
        .insert(snp_pos, (pair_number, seq_pos as usize));
    if !is_supp && frag.seq_string[pair_number as usize].len() == 0 {
        frag.seq_string[pair_number as usize] =
            DnaString::from_acgt_bytes(&record.seq().as_bytes());
        if qual <= 255 - 33 {
            frag.qual_string[pair_number as usize] = record.qual().iter().map(|x| x + 33).collect();
        } else {
            frag.qual_string[pair_number as usize] = record.qual().iter().map(|_x| 255).collect();
        }
    }
    if snp_pos < frag.first_position {
        frag.first_position = snp_pos;
    }
    if snp_pos > frag.last_position {
        frag.last_position = snp_pos
    }
}

#[inline]
pub fn build_truncated_hap_block(
    block: &HapBlock,
    frag: &Frag,
    part: usize,
    current_startpos: SnpPosition,
) -> (FxHashSet<usize>, HapBlock) {
    let ploidy = block.blocks.len();
    let mut blocks_broken = FxHashSet::default();

    //manual deepcopy
    let mut block_vec = block.blocks.clone();
    //count the number of SNPs after the current SNP. This is to check for breaks.
    //this hack needs to be done to circumvent cirular contigs, where weird stuff can happen.
    let mut num_after = vec![0; ploidy];
    let mut num_before = vec![0; ploidy];
    let boundary_end = current_startpos + 50;

    //Add the +50 restriciton for circularity to avoid the case where a fragment covers snps
    //1,2,3,... and n,n-1,n-2,... where n is the number of snps.
    for i in 0..ploidy {
        for pos in block.blocks[i].keys() {
            let boundary_start = *pos + 50;
            if *pos >= current_startpos && *pos < boundary_end {
                num_after[i] += 1;
            }

            if *pos < current_startpos && boundary_start > current_startpos {
                num_before[i] += 1;
            }
            if *pos < current_startpos {
                block_vec[i].remove(pos);
            }
        }
    }

    for i in 0..ploidy {
        if num_after[i] == 0 && num_before[i] != 0 {
            blocks_broken.insert(i);
        }
    }

    for pos in frag.seq_dict.keys() {
        let var_at_pos = frag.seq_dict.get(pos).unwrap();
        let sites = block_vec[part].entry(*pos).or_insert(FxHashMap::default());
        let site_counter = sites.entry(*var_at_pos).or_insert(OrderedFloat(0.));
        *site_counter += utils_frags::phred_scale(frag,pos);
    }

    return (blocks_broken, HapBlock { blocks: block_vec });
}
