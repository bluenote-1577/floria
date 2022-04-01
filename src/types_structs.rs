use fxhash::{FxHashMap, FxHashSet};
use rust_htslib::bam::Record;
use std::hash::{Hash, Hasher};
use std::rc::Rc;
use std::cmp::Ordering;
use debruijn::dna_string::DnaString;

type GnPosition = usize;
type SnpPosition = usize;
type Genotype = usize;

#[derive(Debug, Clone, Default)]
pub struct VcfProfile<'a>{
    pub vcf_pos_allele_map : FxHashMap<&'a str, FxHashMap<i64, Vec<u8>>>,
    pub vcf_pos_to_snp_counter_map : FxHashMap<&'a str, FxHashMap<i64, i64>>,
    pub vcf_snp_pos_to_gn_pos_map : FxHashMap<&'a str, FxHashMap< i64, i64>>,
}

#[derive(Debug, Clone)]
pub struct TraceBackNode{
    pub score: f64,
    pub prev_ind: Option<usize>,
    pub is_sink: bool,
    pub is_source: bool
}
//Positions are inclusive
#[derive(Eq, Debug, Clone, Default)]
pub struct Frag {
    pub id: String,
    pub counter_id: usize,
    pub seq_dict: FxHashMap<SnpPosition, Genotype>,
    pub qual_dict: FxHashMap<SnpPosition, u8>,
    pub positions: FxHashSet<SnpPosition>,
    pub first_position: SnpPosition,
    pub last_position: SnpPosition,
    pub seq_string: Vec<DnaString>,
    pub qual_string: Vec<Vec<u8>>, //Needs to be PHRED scaled i.e. +33 from value
    pub is_paired :bool, 
    pub snp_pos_to_seq_pos: FxHashMap<SnpPosition,(u8, GnPosition)>,
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
    pub current_pos: usize,
    pub broken_blocks: FxHashSet<usize>,
}

impl Ord for SearchNode<'_>{
    fn cmp(&self, other: &Self) -> Ordering {
        self.score.partial_cmp(&other.score).unwrap()
    }
}

impl PartialOrd for SearchNode<'_>{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for SearchNode<'_>{}

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
    pub out_edges: Vec<(usize,f64)>,
    pub in_edges: Vec<(usize,f64)>,
    pub column: usize,
    pub row: usize,
    cov: f64,
    pub id: usize,
    pub out_flows: Vec<(usize,f64)>,
    pub hap_map : FxHashMap<usize, FxHashMap<usize, usize>>,
    pub snp_endpoints : (usize,usize)
}

impl<'a> HapNode<'a> {
    pub fn new(frag_set: FxHashSet<&'a Frag>, snp_endpoints: (usize,usize)) -> HapNode<'a> {
        
        let mut hap_map = FxHashMap::default();
        for frag in frag_set.iter() {
            for pos in frag.positions.iter() {
                if *pos as usize <= snp_endpoints.1 && *pos as usize >= snp_endpoints.0{
                    let var_at_pos = frag.seq_dict.get(pos).unwrap();
                    let sites = hap_map.entry(*pos).or_insert(FxHashMap::default());
                    let site_counter = sites.entry(*var_at_pos).or_insert(0);
                    *site_counter += 1;
                }
            }
        }
        let mut allele_cov_list = vec![];
        for values1 in hap_map.values(){
            for allele_count in values1.values(){
                allele_cov_list.push(*allele_count as f64);
            }
        }

        allele_cov_list.sort_by(|a,b| a.partial_cmp(&b).unwrap());
        let cov;
        if allele_cov_list.is_empty(){
            cov = 0.;
        }
        else{
            cov = allele_cov_list[allele_cov_list.len()*2/3];
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
            snp_endpoints: snp_endpoints
        };
        return toret;

    }

    pub fn cov(&self) -> f64{
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
    current_pos: usize,
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

#[derive(Debug, PartialEq,Eq)]
pub struct HapBlock {
    pub blocks: Vec<FxHashMap<usize, FxHashMap<usize, usize>>>,
}

impl Ord for HapBlock{
    fn cmp(&self, other: &Self) -> Ordering {
        self.blocks.len().partial_cmp(&other.blocks.len()).unwrap()
    }
}

impl PartialOrd for HapBlock{
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
        positions: FxHashSet::default(),
        first_position: usize::MAX,
        last_position: usize::MIN,
        seq_string: vec![DnaString::new();2],
        qual_string: vec![vec![];2],
        is_paired: is_paired,
        snp_pos_to_seq_pos: FxHashMap::default(),
    };

    toret
}
#[inline]
pub fn update_frag(frag: &mut Frag, geno: usize, snp_pos: usize, qual: u8, pair_number: u8, is_supp: bool,  record: &Record, qpos : usize) {
    frag.seq_dict.insert(snp_pos, geno);
    frag.qual_dict.insert(snp_pos, qual);
    frag.positions.insert(snp_pos);
    let mut seq_pos = qpos;
    if is_supp{
        let clipping_offset = record.cigar().leading_hardclips();
        seq_pos = qpos + clipping_offset as usize;
    }
    frag.snp_pos_to_seq_pos.insert(snp_pos, (pair_number, seq_pos as usize));
    if !is_supp && frag.seq_string[pair_number as usize].len() == 0{
        frag.seq_string[pair_number as usize] = DnaString::from_acgt_bytes(&record.seq().as_bytes());
        if qual <= 255 - 33{
        frag.qual_string[pair_number as usize] = record.qual().iter().map(|x| x + 33).collect();
        }
        else{
            frag.qual_string[pair_number as usize] = record.qual().iter().map(|x| 255).collect();
        }
    }
    if snp_pos < frag.first_position {
        frag.first_position = snp_pos;
    }
    if snp_pos > frag.last_position {
        frag.last_position = snp_pos
    }
}

pub fn build_truncated_hap_block(
    block: &HapBlock,
    frag: &Frag,
    part: usize,
    current_startpos: usize,
) -> (FxHashSet<usize>, HapBlock) {
    let ploidy = block.blocks.len();
    let mut blocks_broken = FxHashSet::default();

    //manual deepcopy
    let mut block_vec = block.blocks.clone();
    //count the number of SNPs after the current SNP. This is to check for breaks.
    //this hack needs to be done to circumvent cirular contigs, where weird stuff can happen.
    let mut num_after = vec![0; ploidy];
    let mut num_before = vec![0; ploidy];

    //Add the +50 restriciton for circularity to avoid the case where a fragment covers snps
    //1,2,3,... and n,n-1,n-2,... where n is the number of snps.
    for i in 0..ploidy {
        for pos in block.blocks[i].keys() {
            if *pos >= current_startpos && *pos < current_startpos + 50 {
                num_after[i] += 1;
            }

            if *pos < current_startpos && *pos + 50 > current_startpos {
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

    for pos in frag.positions.iter() {
        let var_at_pos = frag.seq_dict.get(pos).unwrap();
        let sites = block_vec[part].entry(*pos).or_insert(FxHashMap::default());
        let site_counter = sites.entry(*var_at_pos).or_insert(0);
        *site_counter += 1;
    }

    return (blocks_broken, HapBlock { blocks: block_vec });
}
