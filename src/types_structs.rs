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
