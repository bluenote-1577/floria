use fnv::FnvHashMap;
use fnv::FnvHashSet;
use std::hash::{Hash, Hasher};

//Positions are inclusive 
#[derive(PartialEq,Eq,Debug)]
pub struct Frag{
    pub id : String,
    pub counter_id : usize,
    pub seq_dict : FnvHashMap<usize,usize>,
    pub qual_dict : FnvHashMap<usize,u8>,
    pub positions : FnvHashSet<usize>,
    pub first_position : usize,
    pub last_position : usize,
}

impl Hash for Frag {
    fn hash<H: Hasher>(&self, _state: &mut H) {
        self.counter_id;
    }
}

pub struct HapBlock{
    pub blocks: Vec<FnvHashMap<usize,FnvHashMap<usize,usize>>>,
}
