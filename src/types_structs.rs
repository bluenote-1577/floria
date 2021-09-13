use fxhash::{FxHashMap, FxHashSet};
use std::hash::{Hash, Hasher};
use std::ptr;
use std::rc::Rc;

//Positions are inclusive
#[derive(Eq, Debug, Clone)]
pub struct Frag {
    pub id: String,
    pub counter_id: usize,
    pub seq_dict: FxHashMap<usize, usize>,
    pub qual_dict: FxHashMap<usize, u8>,
    pub positions: FxHashSet<usize>,
    pub first_position: usize,
    pub last_position: usize,
    pub supp_aln: Option<String>
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

pub struct SearchNode<'a> {
    pub read: &'a Frag,
    pub part: usize,
    pub score: f64,
    pub freqs: Vec<usize>,
    pub error_vec: Vec<(usize,usize)>,
    pub block_id: usize,
    pub parent_node: Option<Rc<SearchNode<'a>>>
}

pub fn build_child_node<'a>(
    read: &'a Frag,
    part: usize,
    block_id: usize,
    parent_node_option: Option<Rc<SearchNode<'a>>>,
    error_vec: Vec<(usize,usize)>,
    score: f64,
) -> SearchNode<'a> {
    let mut new_freqs = vec![];
    let mut new_parent_node_option: Option<Rc<SearchNode>> = None;
    if let Some(parent_node) = parent_node_option{
        for freq in parent_node.freqs.iter() {
            new_freqs.push(*freq);
        }
        new_parent_node_option = Some(parent_node);
        new_freqs[part] += 1;
    }
    else{
        new_freqs = vec![1;error_vec.len()];
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
    };

    toret
}

pub struct HapBlock {
    pub blocks: Vec<FxHashMap<usize, FxHashMap<usize, usize>>>,
}

pub fn build_frag(id: String, counter_id: usize, supp_aln: Option<String>) -> Frag {
    let toret = Frag {
        id: id,
        counter_id: counter_id,
        seq_dict: FxHashMap::default(),
        qual_dict: FxHashMap::default(),
        positions: FxHashSet::default(),
        first_position: usize::MAX,
        last_position: usize::MIN,
        supp_aln: supp_aln
    };

    toret
}

pub fn update_frag(frag: &mut Frag, geno: usize, qual: u8, snp_pos: usize) {
    frag.seq_dict.insert(snp_pos, geno);
    frag.qual_dict.insert(snp_pos, qual);
    frag.positions.insert(snp_pos);
    if snp_pos < frag.first_position {
        frag.first_position = snp_pos;
    }
    if snp_pos > frag.last_position {
        frag.last_position = snp_pos
    }
}

pub fn build_truncated_hap_block(block: &HapBlock, frag: &Frag, part: usize, current_startpos: usize) -> HapBlock{
    let ploidy = block.blocks.len();
    //manual deepcopy

    let mut block_vec = block.blocks.clone();
    for i in 0..ploidy{
        for pos in block.blocks[i].keys(){
            if *pos < current_startpos{
                block_vec[i].remove(pos);
            }
        }
    }

    for pos in frag.positions.iter(){
        let var_at_pos = frag.seq_dict.get(pos).unwrap();
        let sites = block_vec[part].entry(*pos).or_insert(FxHashMap::default());
        let site_counter = sites.entry(*var_at_pos).or_insert(0);
        *site_counter += 1;
    }

    return HapBlock{blocks : block_vec};

}

