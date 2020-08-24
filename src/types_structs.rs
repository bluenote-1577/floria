use fxhash::{FxHashMap,FxHashSet};
use std::hash::{Hash, Hasher};

//Positions are inclusive 
#[derive(Eq,Debug,Clone)]
pub struct Frag{
    pub id : String,
    pub counter_id : usize,
    pub seq_dict : FxHashMap<usize,usize>,
    pub qual_dict : FxHashMap<usize,u8>,
    pub positions : FxHashSet<usize>,
    pub first_position : usize,
    pub last_position : usize,
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

pub struct HapBlock{
    pub blocks: Vec<FxHashMap<usize,FxHashMap<usize,usize>>>,
}

pub fn build_frag(id : String, counter_id : usize) -> Frag{

    let toret = Frag
    {
        id : id,
        counter_id : counter_id,
        seq_dict : FxHashMap::default(),
        qual_dict : FxHashMap::default(),
        positions : FxHashSet::default(),
        first_position : usize::MAX,
        last_position : usize::MIN,
    };

    toret
}

pub fn update_frag(frag : &mut Frag, geno : usize, qual : u8, snp_pos : usize ){
    frag.seq_dict.insert(snp_pos,geno);
    frag.qual_dict.insert(snp_pos,qual);
    frag.positions.insert(snp_pos);
    if snp_pos < frag.first_position{
        frag.first_position = snp_pos;
    }
    if snp_pos > frag.last_position{
        frag.last_position = snp_pos
    }

}
