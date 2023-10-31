use crate::types_structs::{GenotypeCount};
use ordered_float::OrderedFloat;
pub const NUM_ITER_OPTIMIZE:usize = 20;
pub const MIN_SHARED_READS_UNAMBIG: f64 = 2.;
pub const DIV_FACTOR: f64 = 0.25;
pub const PROB_CUTOFF: f64 = 0.01;
//pub const DIV_FACTOR: f64 = 0.05;
//pub const PROB_CUTOFF: f64 = 0.0001;

pub const HAPQ_CUTOFF: u8 = 0;
pub const MERGE_CUTOFF: f64 = 0.95;
//pub const SMALL_HAPLOGROUP_CUTOFF: usize = 20;
pub const SAME_SNP_DENSITY_CUTOFF: f64 = 1. / 10000.;
pub const DIST_COV_CUTOFF: GenotypeCount = OrderedFloat(0.5);
pub const USE_QUAL_SCORES: bool = true;
pub const MERGE_SIMILAR_HAPLOGROUPS: bool = false;
pub const SEPARATE_BROKEN_HAPLOGROUPS: bool = true;
pub const WEIRD_SPLIT: bool = false;
pub const FLOW_CUTOFF_MULT: f64 = 100.;
pub const HAPQ_CONSTANT: f64 = 40.;
pub const MINIMUM_BLOCK_SIZE: usize = 500;
pub const EXTENSION_BASES: usize = 25;

pub const CONTIG_PLOIDY_HEADER: &str = "contig\taverage_straincount\twhole_contig_multiplicity\tapproximate_coverage_ignoring_indels\ttotal_vartig_bases_covered\taverage_straincount_min15hapq\taverage_straincount_min30hapq\taverage_straincount_min45hapq\tavg_err\n";
