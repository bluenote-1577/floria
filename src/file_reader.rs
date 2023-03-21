use crate::alignment;
use std::ffi::{OsString};
use std::process::Command;
use crate::constants;
use crate::types_structs::{
    build_frag, Frag, Genotype, GnPosition, HapBlock, SnpPosition, VcfProfile, Options
};
use bio_types::genome::AbstractInterval;
use debruijn::dna_string::DnaString;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::IndexedReader;
use bio::io::fasta::{IndexedReader as FastaIndexedReader};
use rust_htslib::{bam, bam::Read as DUMMY_NAME1};
use rust_htslib::{bcf, bcf::Read as DUMMY_NAME2};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::str;
use std::sync::Mutex;

// The output is wrapped in a Result to allow matching on errors
// returns an Iterator to the Reader of the lines of the file.
//
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

// Given a frags.txt file specified as in H-PoP, we return a collection
// (vector) of fragments after processing it.
//
pub fn get_frags_container<P>(filename: P) -> FxHashMap<String, Vec<Frag>>
where
    P: AsRef<Path>,
{
    let mut all_frags = Vec::new();
    let mut counter = 0;

    //Make sure file is able to be read
    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(l) = line {
                let v: Vec<&str> = l.split('\t').collect();

                //First column is the # of blocks
                if let Ok(num_blocks) = v[0].parse::<i32>() {
                    //                    println!("{}",num_blocks);
                    let mut seqs = FxHashMap::default();
                    let mut quals = FxHashMap::default();
                    let mut list_of_positions = Vec::new();
                    let mut first_position = 1;
                    let mut last_position = 1;

                    // For each block, read it into a dictionary with corresp. base
                    for i in 0..num_blocks {
                        let index = i as usize;
                        let start_pos = v[2 * index + 2].parse::<SnpPosition>().unwrap();
                        if i == 0 {
                            first_position = start_pos;
                        }
                        for (j, c) in v[2 * index + 3].chars().enumerate() {
                            let j = j as SnpPosition;
                            seqs.insert(start_pos + j, c.to_digit(10).unwrap() as Genotype);
                            list_of_positions.push(start_pos + j);
                            last_position = start_pos + j
                        }
                    }

                    let qual_string = v.last().unwrap().as_bytes();
                    for (i, key) in list_of_positions.iter().enumerate() {
                        //We usually have a 33 offset for phred qualities. Rust should throw an
                        //error here if this result is negative.
                        quals.insert(*key, qual_string[i] - 33);
                    }

                    let new_frag = Frag {
                        id: v[1].to_string(),
                        counter_id: counter,
                        positions: seqs.keys().map(|x| *x).collect::<FxHashSet<SnpPosition>>(),
                        seq_dict: seqs,
                        qual_dict: quals,
                        first_position: first_position,
                        last_position: last_position,
                        seq_string: vec![DnaString::new(); 2],
                        qual_string: vec![vec![]; 2],
                        is_paired: false,
                        snp_pos_to_seq_pos: FxHashMap::default(),
                    };

                    all_frags.push(new_frag);
                    counter += 1
                } else {
                    panic!("Not a number found in first column");
                }
            }
        }
    }

    let mut frags_map = FxHashMap::default();
    frags_map.insert(String::from("frag_contig"), all_frags);
    frags_map
}



//Read a vcf file to get the genotypes. We read genotypes into a dictionary of keypairs where the
//keys are positions, and the values are dictionaries which encode the genotypes. E.g. the genotype
//1 1 0 0 at position 5 would be (5,{1 : 2, 0 : 2}).
pub fn get_genotypes_from_vcf_hts<P>(
    vcf_file: P,
) -> FxHashMap<String, Vec<usize>>
where
    P: AsRef<Path>,
{
    let mut vcf = match bcf::Reader::from_path(vcf_file) {
        Ok(vcf) => vcf,
        Err(_) => panic!("rust_htslib had an error while reading the VCF file. Exiting."),
    };
    let mut map_positions_vec = FxHashMap::default();
    //let mut positions_vec = Vec::new();
    //let mut genotype_dict = FxHashMap::default();
    let header = vcf.header().clone();

    if header.sample_count() > 1 {
        panic!("More than 1 sample detected in header of vcf file; please use only 1 sample");
    }

    //    if header.contig_count() > 1 {
    //        panic!("More than 1 contig detected in header of vcf file; please use only 1 contig/reference per vcf file.");
    //    }

    let mut last_ref_chrom: &[u8] = &[];

    for rec in vcf.records() {
        let unr = rec.unwrap();
        let alleles = unr.alleles();
        let mut is_snp = true;
        let record_rid = unr.rid().unwrap();
        let ref_chrom_vcf = header.rid2name(record_rid).unwrap();
        if last_ref_chrom != ref_chrom_vcf {
            last_ref_chrom = ref_chrom_vcf;
        }

        for allele in alleles.iter() {
            if allele.len() > 1 {
                is_snp = false;
                break;
            }
        }

        if !is_snp {
            //            println!(
            //                "VCF : Variant at position {} is not a snp. Ignoring.",
            //                unr.pos()
            //            );
            continue;
        }


        let positions_vec = map_positions_vec
            .entry(String::from_utf8(ref_chrom_vcf.to_vec()).unwrap())
            .or_insert(Vec::new());
        //+1 because htslib is 0 index by default
        positions_vec.push(unr.pos() as usize + 1);
    }

    map_positions_vec
}



fn alignment_passed_check(
    flags: u16,
    mapq: u8,
    use_supplementary: bool,
    filter_supplementary: bool,
    mapq_cutoff: u8
) -> (bool, bool) {
    let errors_mask = 1796;
    let secondary_mask = 256;
    let supplementary_mask = 2048;
    let mapq_supp_cutoff = 60;
    let mapq_normal_cutoff = mapq_cutoff;
    let first_in_pair_mask = 64;
    let second_in_pair_mask = 128;

    let is_supp;
    let is_paired = flags & first_in_pair_mask > 0 || flags & second_in_pair_mask > 0;
    if flags & supplementary_mask > 0 {
        is_supp = true;
        //Don't use supplementary alignments for short reads.
        //Increases complexity; maybe fix this in the future.
        if is_paired && is_supp {
            return (false, true);
        } else if !use_supplementary {
            return (false, true);
        } else if filter_supplementary {
            if mapq < mapq_supp_cutoff {
                return (false, true);
            }
        }
    } else {
        is_supp = false;
    }

    if mapq < mapq_normal_cutoff {
        return (false, is_supp);
    }
    //Erroneous alignment, skip
    if flags & errors_mask > 0 {
        //dbg!(&flags,&id_string);
        return (false, is_supp);
    }

    //Secondary alignment, skip
    if flags & secondary_mask > 0 {
        //dbg!(&flags,&id_string);
        return (false, is_supp);
    }

    return (true, is_supp);
}



pub fn get_vcf_profile<'a>(vcf_file: &str, ref_chroms: &'a Vec<String>) -> VcfProfile<'a> {
    let mut vcf_prof = VcfProfile::default();
    let mut vcf = match bcf::Reader::from_path(vcf_file) {
        Ok(vcf) => vcf,
        Err(_) => panic!("rust_htslib had an error reading the VCF file. Exiting."),
    };
    let mut snp_counter = 1;
    let mut vcf_pos_allele_map = FxHashMap::default();
    let mut vcf_pos_to_snp_counter_map = FxHashMap::default();
    let mut vcf_snp_pos_to_gn_pos_map = FxHashMap::default();
    let mut chrom_to_index_map = FxHashMap::default();
    for (i, chrom) in ref_chroms.iter().enumerate() {
        chrom_to_index_map.insert(chrom.as_bytes(), i);
    }

    let vcf_header = vcf.header().clone();

    let mut last_ref_chrom = &String::default();
    for rec in vcf.records() {
        let unr = rec.unwrap();
        let alleles = unr.alleles();
        let mut al_vec = Vec::new();
        let mut is_snp = true;

        let record_rid = unr.rid().unwrap();
        let ref_chrom_vcf =
            //String::from_utf8(vcf_header.rid2name(record_rid).unwrap().to_vec()).unwrap();
            vcf_header.rid2name(record_rid).unwrap();
        let result = chrom_to_index_map.get(&ref_chrom_vcf);
        if result.is_none() {
            continue;
        }
        let contig_name = &ref_chroms[*result.unwrap()];
        //dbg!(String::from_utf8_lossy(ref_chrom_vcf));
        if last_ref_chrom != contig_name {
            snp_counter = 1;
            last_ref_chrom = contig_name;
        }
        let pos_allele_map = vcf_pos_allele_map
            .entry(contig_name.as_str())
            .or_insert(FxHashMap::default());
        let pos_to_snp_counter_map = vcf_pos_to_snp_counter_map
            .entry(contig_name.as_str())
            .or_insert(FxHashMap::default());
        let snp_pos_to_gn_pos_map = vcf_snp_pos_to_gn_pos_map
            .entry(contig_name.as_str())
            .or_insert(FxHashMap::default());

        for allele in alleles.iter() {
            if allele.len() > 1 {
                is_snp = false;
                break;
            }
            al_vec.push(allele[0] as Genotype);
        }

        if !is_snp {
            continue;
        }

        snp_pos_to_gn_pos_map.insert(snp_counter, unr.pos() as GnPosition);
        pos_to_snp_counter_map.insert(unr.pos() as GnPosition, snp_counter);
        snp_counter += 1;
        pos_allele_map.insert(unr.pos() as GnPosition, al_vec);
    }

    vcf_prof.vcf_pos_allele_map = vcf_pos_allele_map;
    vcf_prof.vcf_pos_to_snp_counter_map = vcf_pos_to_snp_counter_map;
    vcf_prof.vcf_snp_pos_to_gn_pos_map = vcf_snp_pos_to_gn_pos_map;
    return vcf_prof;
}

pub fn get_bam_readers(
    options: &Options,
) -> (bam::IndexedReader, Option<bam::IndexedReader>){
    let long_bam_file = &options.bam_file;
    let short_bam_file = &options.short_bam_file;
    let short_bam_read;
    if short_bam_file != ""{
        let short_bam = match bam::IndexedReader::from_path(short_bam_file) {
            Ok(short_bam) => short_bam,
            Err(_) => {
                panic!("rust_htslib had an error while reading the short-read BAM file. Exiting")
            }
        };
        short_bam_read = Some(short_bam);
    }
    else{
        short_bam_read = None;
    }
    let long_bam = match bam::IndexedReader::from_path(long_bam_file) {
        Ok(long_bam) => long_bam,
        Err(_) => panic!("rust_htslib had an error while reading BAM file. Exiting"),
    };

    return (long_bam, short_bam_read);
}

pub fn get_frags_from_bamvcf_rewrite(
    main_bam: &mut bam::IndexedReader,
    short_bam: &mut Option<bam::IndexedReader>,
    vcf_profile: &VcfProfile,
    options: &Options,
    chrom_seqs: &mut Option<FastaIndexedReader<std::fs::File>>,
    contig: &str,
) -> Vec<Frag>
{

    let filter_supplementary = true;
    let use_supplementary = options.use_supp_aln;
    let vcf_pos_allele_map = &vcf_profile.vcf_pos_allele_map;
    let vcf_pos_to_snp_counter_map = &vcf_profile.vcf_pos_to_snp_counter_map;
    let vcf_snp_pos_to_gn_pos_map = &vcf_profile.vcf_snp_pos_to_gn_pos_map;

    let long_bam = main_bam;
    long_bam.fetch(contig).unwrap();

    let mut record_vec_long = vec![];
    for record in long_bam.records() {
        if record.is_ok() {
            record_vec_long.push(record.unwrap());
        }
    }

    let mut record_vec_short = vec![];
    if short_bam.is_some(){
        let short_bam = short_bam.as_mut().unwrap();
        short_bam.fetch(contig).unwrap();
        for record in short_bam.records() {
            if record.is_ok() {
                record_vec_short.push(record.unwrap());
            }
        }
    }
    let mut seq = Vec::new(); 
    if chrom_seqs.is_some(){
        chrom_seqs.as_mut().unwrap().fetch_all(contig).expect("Error reading fasta file.");
        chrom_seqs.as_mut().unwrap().read(&mut seq).expect("Error reading fasta file.");
    }

    let ref_id_to_frag_map: Mutex<FxHashMap<_, _>> = Mutex::new(FxHashMap::default());
    let rec_vecs = vec![record_vec_short, record_vec_long];
    for record_vec in rec_vecs {
        record_vec
            .into_par_iter()
            //            .into_iter()
            .enumerate()
            .for_each(|(count, record)| {
                if record.tid() < 0 {
                } else {
                    let passed_check = alignment_passed_check(
                        record.flags(),
                        record.mapq(),
                        use_supplementary,
                        filter_supplementary,
                        options.mapq_cutoff
                    );

                    //                    log::trace!(
                    //                        "{},{:?}",
                    //                        str::from_utf8(&record.qname()).unwrap(),
                    //                        passed_check
                    //                    );

                    if passed_check.0 {
                        let rec_name: Vec<u8> = record.qname().iter().cloned().collect();
                        let snp_positions_contig = &vcf_pos_to_snp_counter_map[contig];
                        let pos_allele_map = &vcf_pos_allele_map[contig];
                        let snp_to_gn_map = &vcf_snp_pos_to_gn_pos_map[contig];
                        //                        if str::from_utf8(&record.qname()) == Ok("pa1_4940"){
                        //                            dbg!(record.flags(),record.mapq(), passed_check, str::from_utf8(&rec_name));
                        //                        }
                        let mut frag =
                            frag_from_record(&record, snp_positions_contig, pos_allele_map, count);

                        if frag.seq_dict.keys().len() > 0 {
                            if !chrom_seqs.is_none() {
                                alignment::realign(
                                    &seq,
                                    &mut frag,
                                    &snp_to_gn_map,
                                    &pos_allele_map,
                                );
                            }
                            let mut locked = ref_id_to_frag_map.lock().unwrap();
                            let bucket = locked.entry(rec_name).or_insert(vec![]);
                            bucket.push((record.flags(), frag));
                        }
                    }
                }
            });
    }

    let ref_vec_frags = combine_frags(
        ref_id_to_frag_map.into_inner().unwrap(),
        &vcf_profile,
        contig,
    );

    //    for frag in ref_vec_frags.iter(){
    //        if frag.id.contains("485041"){
    //            dbg!(&frag);
    //        }
    //    }
    ref_vec_frags
}

pub fn get_fasta_seqs(fasta_file: &str) -> FastaIndexedReader<std::fs::File> {
    let mut os_string: OsString = Path::new(fasta_file).into();
    os_string.push(".");
    os_string.push("fai");
    let fai_path: &str = os_string.to_str().unwrap();
    if !Path::new(fai_path).exists(){
        log::warn!(".fai index not detected. Trying to index using samtools faidx if it is in PATH.");
        let status = Command::new("samtools")
            .arg("faidx")
            .arg(fasta_file)
            .status()
            .expect("Failed to run 'samtools index' on fasta file.");
        if !status.success() {
            log::error!("samtools faidx failed. Exiting");
            std::process::exit(1);
        }
    }
    let reader = FastaIndexedReader::from_file(&fasta_file.to_string());
    if !reader.is_err() {
        return reader.unwrap();
    }
    else{
        log::error!("Could not read fasta file. Exiting.");
        std::process::exit(1);
    }
}

fn combine_frags(
    id_to_frag_map: FxHashMap<Vec<u8>, Vec<(u16, Frag)>>,
    vcf_profile: &VcfProfile,
    contig: &str,
) -> Vec<Frag> {
    let first_in_pair_mask = 64;
    let second_in_pair_mask = 128;
    let supplementary_mask = 2048;

    let mut ref_frags = vec![];
    for (_id, mut frags) in id_to_frag_map {
        //        dbg!(&str::from_utf8(&_id), frags.len(), frags.iter().map(|x| x.0).collect::<Vec<u16>>());
        //paired
        if frags.len() == 2 && frags[0].1.is_paired && frags[1].1.is_paired {
            //            log::trace!(
            //                "{}, {},{},{},{}",
            //                frags[0].1.id,
            //                frags[0].0,
            //                frags[0].1.first_position,
            //                frags[1].0,
            //                frags[1].1.first_position
            //            );
            frags.sort();
            let first = std::mem::take(&mut frags[0]);
            let second = std::mem::take(&mut frags[1]);
            //            log::trace!(
            //                "{},{},{},{}",
            //                first.1.id,
            //                second.1.id,
            //                first.1.counter_id,
            //                second.1.counter_id
            //            );

            let mut first_frag;
            let mut sec_frag;

            if first.0 & first_in_pair_mask == first_in_pair_mask {
                first_frag = first.1;
                sec_frag = second.1
            } else if first.0 & second_in_pair_mask == second_in_pair_mask {
                first_frag = second.1;
                sec_frag = first.1
            } else {
                log::warn!("Read {} is not paired and has more than one primary alignment; something went wrong.", first.1.id);
                continue;
            }

            first_frag.seq_dict.extend(sec_frag.seq_dict);
            first_frag.qual_dict.extend(sec_frag.qual_dict);
            first_frag.positions.extend(sec_frag.positions);

            first_frag.first_position =
                SnpPosition::min(first_frag.first_position, sec_frag.first_position);
            first_frag.last_position =
                SnpPosition::max(first_frag.last_position, sec_frag.last_position);

            let mut temp = DnaString::new();
            std::mem::swap(&mut sec_frag.seq_string[0], &mut temp);
            first_frag.seq_string[1] = temp;
            first_frag.qual_string[1] = std::mem::take(&mut sec_frag.qual_string[0]);

            for (_snp_pos, (read_pair, _pos)) in sec_frag.snp_pos_to_seq_pos.iter_mut() {
                *read_pair = 1;
            }

            first_frag
                .snp_pos_to_seq_pos
                .extend(sec_frag.snp_pos_to_seq_pos);
            ref_frags.push(first_frag);
        } else if frags.len() == 1 && frags[0].0 & supplementary_mask == 0 {
            ref_frags.push(std::mem::take(&mut frags[0].1));
        } else {
            //Arbitrary cutoff for reference suppl. alignment distance
            let supp_aln_dist_cutoff = constants::SUPPL_ALN_DIST_CUTOFF;
            //2 or more fragments and no paired indicates a long supplementary alignment.
            for frag in frags.iter() {
                //                dbg!(
                //                    &frag.1.id,
                //                    &frag.1.first_position,
                //                    &frag.1.last_position,
                //                    &frag.1.snp_pos_to_seq_pos[&frag.1.first_position],
                //                    &frag.1.snp_pos_to_seq_pos[&frag.1.last_position]
                //                );
                if frag.1.is_paired{
                    log::warn!("Fragment {} is paired but appears in more than two mappings -- possibly a supplementary alignment. Careful.", &frag.1.id);
                }
            }

            let mut supp_intervals = vec![];

            for frag in frags.iter() {
                supp_intervals.push((frag.1.first_position, frag.1.last_position));
            }
            supp_intervals.sort();

            let snp_to_gn = &vcf_profile.vcf_snp_pos_to_gn_pos_map[contig];
            let mut take_primary_only = false;
            for i in 0..supp_intervals.len() - 1 {
                if snp_to_gn[&supp_intervals[i + 1].0] as i64
                    - snp_to_gn[&supp_intervals[i].1] as i64
                    > supp_aln_dist_cutoff
                {
                    take_primary_only = true;
                    break;
                }
            }

            let mut primary_alignment_index = None;
            for (i, frag) in frags.iter().enumerate() {
                if frag.0 & supplementary_mask != supplementary_mask {
                    if !primary_alignment_index.is_none() {
                        log::warn!("More than one primary alignment for read {}. Only one primary alignment allowed
                            per read unless paired. Using arbitrary primary alignment for read.", frag.1.id);
                    }
                    primary_alignment_index = Some(i);
                }
            }
            //Only suppl. alignments. Primary probably
            //got filtered out, but suppl didn't. Don't do anything.
            if primary_alignment_index.is_none() {
                //                for frag in frags.iter() {
                //                    dbg!(frag.0, &frag.1.id);
                //                }
                continue;
            }
            if take_primary_only {
                ref_frags.push(std::mem::take(
                    &mut frags[primary_alignment_index.unwrap()].1,
                ));
            } else {
                let mut primary_frag =
                    std::mem::take(&mut frags[primary_alignment_index.unwrap()].1);
                for i in 0..frags.len() {
                    if i == primary_alignment_index.unwrap() {
                        continue;
                    }
                    let frag = std::mem::take(&mut frags[i].1);
                    primary_frag.seq_dict.extend(frag.seq_dict);
                    primary_frag.qual_dict.extend(frag.qual_dict);
                    primary_frag.positions.extend(frag.positions);

                    primary_frag.first_position =
                        SnpPosition::min(primary_frag.first_position, frag.first_position);
                    primary_frag.last_position =
                        SnpPosition::max(primary_frag.last_position, frag.last_position);

                    primary_frag
                        .snp_pos_to_seq_pos
                        .extend(frag.snp_pos_to_seq_pos);
                }

                ref_frags.push(primary_frag);
            }
        }
    }
    return ref_frags;
}

fn frag_from_record(
    record: &bam::Record,
    snp_positions: &FxHashMap<GnPosition, SnpPosition>,
    pos_allele_map: &FxHashMap<GnPosition, Vec<Genotype>>,
    counter_id: usize,
) -> Frag {
    let first_in_pair_mask = 64;
    let second_in_pair_mask = 128;
    let supplementary_mask = 2048;
    let mut leading_hardclips = 0;
    let paired =
        (record.flags() & first_in_pair_mask > 0) || (record.flags() & second_in_pair_mask > 0);
    let aligned_pairs = record.aligned_pairs_full();
    let mut _last_read_aligned_pos = 0;
    let mut frag = build_frag(
        String::from_utf8(record.qname().to_vec()).unwrap(),
        counter_id,
        paired,
    );
    if record.flags() & supplementary_mask > 0 {
        leading_hardclips = record.cigar().leading_hardclips();
    }

    for pair in aligned_pairs {
        if pair[1].is_none() {
            continue;
        }
        let genome_pos = pair[1].unwrap() as GnPosition;
        if !snp_positions.contains_key(&genome_pos) {
            if !pair[0].is_none() {
                _last_read_aligned_pos = pair[0].unwrap() as GnPosition;
            }
            continue;
        } else {
            //Deletion
            if pair[0].is_none() {
            } else {
                let seq_pos = pair[0].unwrap() as GnPosition;
                let readbase = record.seq()[seq_pos] as Genotype;
                for (i, allele) in pos_allele_map
                    .get(&(genome_pos))
                    .unwrap()
                    .iter()
                    .enumerate()
                {
                    if readbase == *allele {
                        let snp_pos = snp_positions[&genome_pos] as SnpPosition;
                        frag.seq_dict.insert(snp_pos, i as Genotype);
                        frag.qual_dict.insert(snp_pos, record.qual()[seq_pos]);
                        if snp_pos < frag.first_position {
                            frag.first_position = snp_pos;
                        }
                        if snp_pos > frag.last_position {
                            frag.last_position = snp_pos
                        }
                        //Long read assumption.
                        frag.snp_pos_to_seq_pos
                            .insert(snp_pos, (0, seq_pos + leading_hardclips as usize));
                        break;
                    }
                }
            }
        }
    }

    frag.seq_string[0] = DnaString::from_acgt_bytes(&record.seq().as_bytes());
    frag.positions = frag
        .seq_dict
        .keys()
        .map(|x| *x)
        .collect::<FxHashSet<SnpPosition>>();
    frag.qual_string[0] = record.qual().iter().map(|x| x + 33).collect();
    return frag;
}

pub fn get_contigs_to_phase(bam_file: &str) -> Vec<String> {
    let bam = IndexedReader::from_path(bam_file).unwrap();
    return bam
        .header()
        .target_names()
        .iter()
        .map(|x| String::from_utf8(x.to_vec()).unwrap())
        .collect();
}

