use crate::constants;
use crate::part_block_manip;
use crate::types_structs::*;
use crate::utils_frags;
use bio::alphabets::dna::revcomp;
use bio::io::fastq;
use bio::io::fastq::Writer;
use debruijn::*;
use flate2::write::GzEncoder;
use flate2::Compression;
use fxhash::{FxHashMap, FxHashSet};
use std::collections::BTreeMap;
use std::fs;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::LineWriter;
use std::io::Write;
use std::path::Path;

pub fn write_outputs(
    part: &Vec<FxHashSet<&Frag>>,
    snp_range_parts_vec: &Vec<(SnpPosition, SnpPosition)>,
    out_bam_part_dir: String,
    prefix: &String,
    contig: &String,
    snp_pos_to_genome_pos: &Vec<usize>,
    options: &Options,
    snpless_frags: &Vec<&Frag>,
) {
    let trim_reads = options.trim_reads;
    let gzip = options.gzip;
    fs::create_dir_all(&out_bam_part_dir).unwrap();

    if !snp_range_parts_vec.is_empty() {
        fs::create_dir_all(&format!("{}/local_parts", out_bam_part_dir)).unwrap();
        fs::create_dir_all(&format!("{}/short_reads", out_bam_part_dir)).unwrap();
        fs::create_dir_all(&format!("{}/long_reads", out_bam_part_dir)).unwrap();
        fs::create_dir_all(&format!("{}/vartig_info", out_bam_part_dir)).unwrap();
    }

    let (hapqs, rel_err, avg_err) =
        part_block_manip::get_hapq(&part, snp_pos_to_genome_pos, snp_range_parts_vec, options);
    write_haplotypes(
        part,
        contig,
        snp_range_parts_vec,
        &out_bam_part_dir,
        snp_pos_to_genome_pos,
        &hapqs,
        &rel_err,
        &options.out_dir,
        avg_err
    );
    write_all_parts_file(
        part,
        contig,
        snp_range_parts_vec,
        &out_bam_part_dir,
        prefix,
        snp_pos_to_genome_pos,
        &hapqs,
        &rel_err,
    );
    write_nosnp_reads_parts(
        &out_bam_part_dir,
        &snpless_frags,
    );
    if options.output_reads {
        write_reads(
            part,
            snp_range_parts_vec,
            &out_bam_part_dir,
            !trim_reads,
            &hapqs,
            gzip,
        );
        write_nosnp_reads(
            &out_bam_part_dir,
            &snpless_frags,
            gzip,
        );

    }
}

fn write_nosnp_reads(out_bam_part_dir: &str, snpless_frags:&Vec<&Frag>, gzip: bool){
    let gz = if gzip { ".gz" } else { "" };
    let part_fastq_reads = format!("{}/long_reads/snpless.fastq{}", out_bam_part_dir, gz);
    let part_fastq_reads_paired1 = format!(
        "{}/short_reads/snpless_paired1.fastq{}",
        out_bam_part_dir, gz
    );
    let part_fastq_reads_paired2 = format!(
        "{}/short_reads/snpless_paired2.fastq{}",
        out_bam_part_dir, gz
    );

    let fastq_file = File::create(part_fastq_reads).expect("Can't create file");
    let fastq_file1 = File::create(part_fastq_reads_paired1).expect("Can't create file");
    let fastq_file2 = File::create(part_fastq_reads_paired2).expect("Can't create file");
    let mut fastq_writer;
    let mut fastq_writer_paired1;
    let mut fastq_writer_paired2;

    if gzip {
        let gz_encoder: Box<dyn Write> =
            Box::new(GzEncoder::new(fastq_file, Compression::default()));
        let gz_encoder1: Box<dyn Write> =
            Box::new(GzEncoder::new(fastq_file1, Compression::default()));
        let gz_encoder2: Box<dyn Write> =
            Box::new(GzEncoder::new(fastq_file2, Compression::default()));

        fastq_writer = fastq::Writer::new(gz_encoder);
        fastq_writer_paired1 = fastq::Writer::new(gz_encoder1);
        fastq_writer_paired2 = fastq::Writer::new(gz_encoder2);
    } else {
        let enc: Box<dyn Write> = Box::new(fastq_file);
        let enc1: Box<dyn Write> = Box::new(fastq_file1);
        let enc2: Box<dyn Write> = Box::new(fastq_file2);
        fastq_writer = fastq::Writer::new(enc);
        fastq_writer_paired1 = fastq::Writer::new(enc1);
        fastq_writer_paired2 = fastq::Writer::new(enc2);
    }

    for frag in snpless_frags{
        if frag.is_paired{
            write_paired_reads_no_trim(&mut fastq_writer_paired1, &mut fastq_writer_paired2, &frag);
        }
        else{
            if frag.seq_string[0].len() == 0{
                fastq_writer.write(&format!("{}",frag.id), None, &vec![78], &vec![33]).unwrap();
            }
            else{
                fastq_writer.write(&format!("{}",frag.id), None, &frag.seq_string[0].to_ascii_vec(), &frag.qual_string[0]).unwrap();
            }
        }
    }
}

fn write_nosnp_reads_parts(out_bam_part_dir: &str, snpless_frags:&Vec<&Frag>){
    let part_path = &format!("{}/reads_without_snps.tsv", out_bam_part_dir);
    let file = File::create(part_path).expect("Can't create file");

    let mut file = LineWriter::new(file);

    write!(file, "READ_NAME\tREAD_LENGTH_IN_BASES\n").unwrap();
    for frag in snpless_frags{
        let mut len = 0;
        for dna_string in frag.seq_string.iter(){
            len += dna_string.len();
        }
        write!(file, "{}\t{}\n", &frag.id, len).unwrap();
    }

}

fn write_paired_reads_no_trim<W: Write>(
    fastq_writer_paired1: &mut Writer<W>,
    fastq_writer_paired2: &mut Writer<W>,
    frag: &Frag,
) {
    if frag.seq_string[0].len() == 0 {
        fastq_writer_paired1
            .write(
                &format!("{}/1", frag.id),
                None,
                //Write N instead
                &vec![78],
                &vec![33],
            )
            .unwrap();
    } else {
        fastq_writer_paired1
            .write(
                &format!("{}/1", frag.id),
                None,
                &frag.seq_string[0].to_ascii_vec(),
                &frag.qual_string[0],
            )
            .unwrap();
    }
    if frag.seq_string[1].len() == 0 {
        fastq_writer_paired2
            .write(
                &format!("{}/2", frag.id),
                None,
                //Write N instead
                &vec![78],
                &vec![33],
            )
            .unwrap();
    } else {
        fastq_writer_paired2
            .write(
                &format!("{}/2", frag.id),
                None,
                &frag.seq_string[1].rc().to_ascii_vec(),
                &frag.qual_string[1],
            )
            .unwrap();
    }
}

fn _write_paired_reads<W: Write>(
    fastq_writer_paired1: &mut Writer<W>,
    fastq_writer_paired2: &mut Writer<W>,
    left_read_pair: u8,
    right_read_pair: u8,
    left_seq_pos: usize,
    right_seq_pos: usize,
    frag: &Frag,
) {
    if left_read_pair == right_read_pair {
        let writer;
        let other_writer;
        let read_pair;
        let other_read_pair;
        if left_read_pair == 0 {
            read_pair = 0;
            other_read_pair = 1;
            writer = fastq_writer_paired1;
            other_writer = fastq_writer_paired2;
        } else {
            read_pair = 1;
            other_read_pair = 0;
            writer = fastq_writer_paired2;
            other_writer = fastq_writer_paired1;
        }
        writer
            .write(
                &format!("{}/{}", frag.id, read_pair),
                None,
                &frag.seq_string[read_pair as usize].to_ascii_vec()
                    [left_seq_pos..right_seq_pos + 1],
                &frag.qual_string[read_pair as usize].as_slice()[left_seq_pos..right_seq_pos + 1],
            )
            .unwrap();
        other_writer
            .write(
                &format!("{}/{}", frag.id, other_read_pair),
                None,
                //Write N instead
                &vec![78],
                &vec![33],
            )
            .unwrap();
    } else {
        if frag.seq_string[left_read_pair as usize].len() == 0 {
            fastq_writer_paired1
                .write(
                    &format!("{}/1", frag.id),
                    None,
                    //Write N instead
                    &vec![78],
                    &vec![33],
                )
                .unwrap();
        } else {
            fastq_writer_paired1
                .write(
                    &format!("{}/1", frag.id),
                    None,
                    &frag.seq_string[left_read_pair as usize].to_ascii_vec()[left_seq_pos..],
                    &frag.qual_string[left_read_pair as usize].as_slice()[left_seq_pos..],
                )
                .unwrap();
        }
        if frag.seq_string[right_read_pair as usize].len() == 0 {
            fastq_writer_paired2
                .write(
                    &format!("{}/2", frag.id),
                    None,
                    //Write N instead
                    &vec![78],
                    &vec![33],
                )
                .unwrap();
        } else {
            let qual_cut_string =
                &frag.qual_string[right_read_pair as usize].as_slice()[..right_seq_pos];
            let rev_quals: Vec<u8> = qual_cut_string.into_iter().rev().map(|x| *x).collect();
            fastq_writer_paired2
                .write(
                    &format!("{}/2", frag.id),
                    None,
                    &revcomp(
                        &frag.seq_string[right_read_pair as usize].to_ascii_vec()[..right_seq_pos],
                    ),
                    //TODO Do we need to flip this as well?
                    rev_quals.as_slice(),
                )
                .unwrap();
        }
    }
}

fn write_fragset_haplotypes(
    frags: &FxHashSet<&Frag>,
    name: &str,
    dir: &str,
    snp_pos_to_genome_pos: &Vec<GnPosition>,
    _append: bool,
    left_snp_pos: SnpPosition,
    right_snp_pos: SnpPosition,
) -> Vec<u8> {
    let filename = format!("{}/vartig_info/{}_hap.txt", dir, name);
    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(filename)
        .unwrap();

    let hap_map = utils_frags::set_to_seq_dict(&frags, false);
    let emptydict = FxHashMap::default();
    let title_string = format!(">HAP{}.{}\tSNPRANGE:{}-{}\n", name, dir, left_snp_pos, right_snp_pos);
    write!(file, "{}", title_string).unwrap();
    let positions: Vec<&SnpPosition> = hap_map.keys().collect();
    if positions.len() == 0 {
        return vec![];
    }
    let mut vec_of_alleles = vec![];
    for pos in left_snp_pos..right_snp_pos + 1 {
        if snp_pos_to_genome_pos.len() == 0 {
            write!(file, "{}:NA\t", pos).unwrap();
        } else {
            write!(
                file,
                "{}:{}\t",
                pos,
                snp_pos_to_genome_pos[(pos - 1) as usize]
            )
            .unwrap();
        }
        let allele_map = hap_map.get(&pos).unwrap_or(&emptydict);
        //If a block has no coverage at a position, we write ?.
        if *allele_map == emptydict {
            file.write_all(b"?\t").unwrap();
            //This prints ?
            vec_of_alleles.push(15);
        } else {
            let best_allele = allele_map.iter().max_by_key(|entry| entry.1).unwrap().0;
            write!(file, "{}\t", best_allele).unwrap();
            vec_of_alleles.push(*best_allele as u8);
        }

        if *allele_map == emptydict {
            write!(file, "NA\t").unwrap();
        } else {
            let mut first = true;
            for (site, count) in allele_map {
                if !first {
                    write!(file, "|").unwrap();
                }
                if first {
                    first = false;
                }
                write!(file, "{}:{}", site, count.round() as usize).unwrap();
            }
            write!(file, "\t").unwrap();
        }
        write!(file, "\n").unwrap();
    }
    return vec_of_alleles;
}

fn write_reads(
    part: &Vec<FxHashSet<&Frag>>,
    snp_range_parts_vec: &Vec<(SnpPosition, SnpPosition)>,
    out_bam_part_dir: &String,
    extend_read_clipping: bool,
    hapqs: &Vec<u8>,
    gzip: bool,
) {
    for (i, set) in part.iter().enumerate() {
        if set.is_empty() {
            continue;
        }
        if snp_range_parts_vec.is_empty() {
            continue;
        }
        if hapqs[i] < constants::HAPQ_CUTOFF {
            continue;
        }
        let left_snp_pos = snp_range_parts_vec[i].0;
        let right_snp_pos = snp_range_parts_vec[i].1;
        //Populate all_part.txt file
        let mut vec_part: Vec<&Frag> = set.into_iter().cloned().collect();
        vec_part.sort();
        //Non-empty means that we're writing the final partition after path collection
        let gz = if gzip { ".gz" } else { "" };
        let part_fastq_reads = format!("{}/long_reads/{}_part.fastq{}", out_bam_part_dir, i, gz);
        let part_fastq_reads_paired1 = format!(
            "{}/short_reads/{}_part_paired1.fastq{}",
            out_bam_part_dir, i, gz
        );
        let part_fastq_reads_paired2 = format!(
            "{}/short_reads/{}_part_paired2.fastq{}",
            out_bam_part_dir, i, gz
        );

        let fastq_file = File::create(part_fastq_reads).expect("Can't create file");
        let fastq_file1 = File::create(part_fastq_reads_paired1).expect("Can't create file");
        let fastq_file2 = File::create(part_fastq_reads_paired2).expect("Can't create file");
        let mut fastq_writer;
        let mut fastq_writer_paired1;
        let mut fastq_writer_paired2;

        if gzip {
            let gz_encoder: Box<dyn Write> =
                Box::new(GzEncoder::new(fastq_file, Compression::default()));
            let gz_encoder1: Box<dyn Write> =
                Box::new(GzEncoder::new(fastq_file1, Compression::default()));
            let gz_encoder2: Box<dyn Write> =
                Box::new(GzEncoder::new(fastq_file2, Compression::default()));

            fastq_writer = fastq::Writer::new(gz_encoder);
            fastq_writer_paired1 = fastq::Writer::new(gz_encoder1);
            fastq_writer_paired2 = fastq::Writer::new(gz_encoder2);
        } else {
            let enc: Box<dyn Write> = Box::new(fastq_file);
            let enc1: Box<dyn Write> = Box::new(fastq_file1);
            let enc2: Box<dyn Write> = Box::new(fastq_file2);
            fastq_writer = fastq::Writer::new(enc);
            fastq_writer_paired1 = fastq::Writer::new(enc1);
            fastq_writer_paired2 = fastq::Writer::new(enc2);
        }

        //        let mut fastq_writer = fastq::Writer::new(fastq_file);
        //        let mut fastq_writer_paired1 = fastq::Writer::new(fastq_file1);
        //        let mut fastq_writer_paired2 = fastq::Writer::new(fastq_file2);
        let extension = 25;
        for frag in vec_part.iter() {
            let mut found_primary = false;
            for seq in frag.seq_string.iter() {
                if seq.len() != 0 {
                    found_primary = true;
                    break;
                }
            }
            if !found_primary {
                log::trace!(
                    "{} primary not found. Paired: {}",
                    &frag.id,
                    &frag.is_paired
                );
                continue;
            }
            //We can merge haplogroups such that small reads may fall off the
            //snp range intervals
            if frag.first_position > right_snp_pos {
                log::trace!("Read {} past right endpoint; trim until empty.", &frag.id);
                continue;
            }
            if frag.last_position < left_snp_pos {
                log::trace!("Read {} past left endpoint; trim until empty.", &frag.id);
                continue;
            }
            let mut left_seq_pos;
            let mut tmp = left_snp_pos;
            let _left_read_pair;
            if frag.first_position > left_snp_pos && extend_read_clipping {
                left_seq_pos = 0;
                _left_read_pair = 0;
            } else {
                loop {
                    if frag.snp_pos_to_seq_pos.contains_key(&tmp) {
                        let info = frag.snp_pos_to_seq_pos[&tmp];
                        left_seq_pos = info.1;
                        _left_read_pair = info.0;
                        break;
                    }
                    tmp += 1;
                    if tmp - left_snp_pos > 10000000 {
                        dbg!(
                            &frag.first_position,
                            &frag.last_position,
                            left_snp_pos,
                            right_snp_pos,
                            &frag.snp_pos_to_seq_pos,
                        );
                        panic!("left snp position of partition for the read was not found.");
                    }
                }
            }
            if left_seq_pos > extension {
                left_seq_pos -= extension;
            } else {
                left_seq_pos = 0;
            }

            let mut right_seq_pos;
            let mut tmp = right_snp_pos;
            let right_read_pair;
            if frag.last_position < right_snp_pos && extend_read_clipping {
                if frag.is_paired {
                    right_read_pair = 1;
                } else {
                    right_read_pair = 0;
                }
                if frag.seq_string[right_read_pair as usize].len() == 0 {
                    right_seq_pos = 0;
                } else {
                    right_seq_pos = frag.seq_string[right_read_pair as usize].len() - 1;
                }
            } else {
                loop {
                    if frag.snp_pos_to_seq_pos.contains_key(&tmp) {
                        let info = frag.snp_pos_to_seq_pos[&tmp];
                        right_seq_pos = info.1;
                        right_read_pair = info.0;
                        break;
                    }
                    if tmp == 0 {
                        dbg!(left_snp_pos, right_snp_pos);
                    }
                    tmp -= 1;
                }
            }

            if frag.seq_string[right_read_pair as usize].len() == 0 {
                right_seq_pos = 0;
            } else if frag.seq_string[right_read_pair as usize].len() > extension + 1
                && right_seq_pos < frag.seq_string[right_read_pair as usize].len() - extension - 1
            {
                right_seq_pos += extension;
            } else {
                right_seq_pos = frag.seq_string[right_read_pair as usize].len() - 1;
            }

            if frag.is_paired {
                write_paired_reads_no_trim(
                    &mut fastq_writer_paired1,
                    &mut fastq_writer_paired2,
                    &frag,
                );
            } else {
                if left_seq_pos > right_seq_pos {
                    log::trace!(
                            "{} left seq pos > right seq pos at {:?}. Left:{}, Right:{}.
                            This usually happens when a read id is not unique. May happen with suppl. alignments too.",
                            &frag.id, snp_range_parts_vec[i], left_seq_pos, right_seq_pos
                        );
                    continue;
                }
                fastq_writer
                    .write(
                        &frag.id,
                        None,
                        &frag.seq_string[0].to_ascii_vec()[left_seq_pos..right_seq_pos + 1],
                        &frag.qual_string[0].as_slice()[left_seq_pos..right_seq_pos + 1],
                    )
                    .unwrap();
            }
        }
    }
}

//Write a vector of blocks into a file.
pub fn write_blocks_to_file<P>(
    out_part_dir: P,
    blocks: &Vec<HapBlock>,
    lengths: &Vec<usize>,
    snp_to_genome: &Vec<usize>,
    part: &Vec<FxHashSet<&Frag>>,
    _first_iter: bool,
    contig: &String,
    break_positions: &FxHashMap<SnpPosition, FxHashSet<SnpPosition>>,
) where
    P: AsRef<Path>,
{
    let ploidy = blocks[0].blocks.len();
    let filename = out_part_dir
        .as_ref()
        .join(format!("{}_phasing.txt", contig));

    dbg!(&filename);
    let file;
    file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(filename)
        .unwrap();
    //let file = File::create(filename).expect("Can't create file");
    let mut file = LineWriter::new(file);
    let mut length_prev_block = 1;
    let emptydict = FxHashMap::default();
    let unpolished_block = utils_frags::hap_block_from_partition(part, true);
    //dbg!(snp_to_genome.len(),lengths[0] + 1);

    for (i, block) in blocks.iter().enumerate() {
        let title_string = format!("**{}**\n", contig);
        file.write_all(title_string.as_bytes()).unwrap();
        for pos in length_prev_block..length_prev_block + lengths[i] {
            let pos = pos as SnpPosition;
            if break_positions.contains_key(&pos) {
                write!(file, "--------\n").unwrap();
            }
            if snp_to_genome.len() == 0 {
                write!(file, "{}:NA\t", pos).unwrap();
            } else {
                write!(file, "{}:{}\t", pos, snp_to_genome[(pos - 1) as usize]).unwrap();
            }
            //Write haplotypes
            for k in 0..ploidy {
                let allele_map = block.blocks[k].get(&pos).unwrap_or(&emptydict);
                //If a block has no coverage at a position, we write -1.
                if *allele_map == emptydict {
                    file.write_all(b"-1\t").unwrap();
                } else {
                    let best_allele = allele_map.iter().max_by_key(|entry| entry.1).unwrap().0;
                    write!(file, "{}\t", best_allele).unwrap();
                }
            }

            //Write stats
            for k in 0..ploidy {
                let allele_map_unpolish =
                    unpolished_block.blocks[k].get(&pos).unwrap_or(&emptydict);
                if *allele_map_unpolish == emptydict {
                    write!(file, "NA\t").unwrap();
                } else {
                    let mut first = true;
                    for (site, count) in allele_map_unpolish {
                        if !first {
                            write!(file, "|").unwrap();
                        }
                        if first {
                            first = false;
                        }
                        write!(file, "{}:{}", site, count).unwrap();
                    }
                    write!(file, "\t").unwrap();
                }
            }
            write!(file, "\n").unwrap();
        }
        write!(file, "*****\n").unwrap();
        length_prev_block += lengths[i]
    }
}

//Write a vector of sorted fragment files by first position (no guarantees on end position) to a
//file in the same format as H-PoP and other haplotypers.
pub fn write_frags_file(frags: Vec<Frag>, filename: String) {
    let file = File::create(filename).expect("Can't create file");
    let mut file = LineWriter::new(file);
    for frag in frags.into_iter() {
        let frag_id = frag.id.clone();
        let (start_vec, blocks, qual_block) = convert_dict_to_block(frag);
        if start_vec.len() != blocks.len() {
            dbg!(start_vec.len(), blocks.len());
            panic!("Block length diff");
        }

        write!(file, "{}\t", blocks.len()).unwrap();
        write!(file, "{}\t", frag_id).unwrap();
        for i in 0..blocks.len() {
            write!(file, "{}\t", start_vec[i]).unwrap();
            for var in blocks[i].iter() {
                write!(file, "{}", *var).unwrap();
            }
            write!(file, "\t").unwrap();
        }

        for q in qual_block.iter() {
            if *q as usize + 33 > 255 {
                write!(file, "{}", (*q) as char).unwrap();
            } else {
                write!(file, "{}", (*q + 33) as char).unwrap();
            }
        }

        write!(file, "\n").unwrap();
    }
}

#[allow(non_snake_case)]
fn write_haplotypes(
    part: &Vec<FxHashSet<&Frag>>,
    contig: &String,
    snp_range_parts_vec: &Vec<(SnpPosition, SnpPosition)>,
    out_bam_part_dir: &String,
    snp_pos_to_genome_pos: &Vec<usize>,
    hapqs: &Vec<u8>,
    rel_err: &Vec<f64>,
    top_dir: &str,
    avg_err: f64,
) -> FxHashMap<usize, u8> {
    let vartig_file = format!("{}/{}.vartigs", out_bam_part_dir, contig);
    let ploidy_file = format!("{}/ploidy_info.tsv", out_bam_part_dir);
    let top_ploidy_file = format!("{}/contig_ploidy_info.tsv", top_dir);
    let mut longest_vartig_bases = 0;

    let mut snp_covered_count = vec![0.; snp_pos_to_genome_pos.len()];
    let mut coverage_count = vec![0.; snp_pos_to_genome_pos.len()];

    let mut snp_covered_count_g0 = vec![0.; snp_pos_to_genome_pos.len()];
    let mut coverage_count_g0 = vec![0.; snp_pos_to_genome_pos.len()];

    let mut hapQ_scores = FxHashMap::default();
    let mut total_bases_covered = 0;
    //    let mut hapQ_scores = vec![];
    let mut vartig_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(vartig_file)
        .unwrap();

    for (i, set) in part.iter().enumerate() {
        if set.is_empty() {
            continue;
        }

        if !snp_range_parts_vec.is_empty() {
            let left_snp_pos = snp_range_parts_vec[i].0;
            let right_snp_pos = snp_range_parts_vec[i].1;
            if left_snp_pos > right_snp_pos {
                dbg!(&snp_range_parts_vec[i], contig);
                panic!();
            }
            let left_gn_pos = snp_pos_to_genome_pos[(left_snp_pos - 1) as usize];
            let right_gn_pos = snp_pos_to_genome_pos[(right_snp_pos - 1) as usize];

            let bases_covered_vartig = right_gn_pos - left_gn_pos;
            total_bases_covered += bases_covered_vartig;
            if bases_covered_vartig > longest_vartig_bases {
                longest_vartig_bases = bases_covered_vartig;
            }

            let (cov, err, _total_err, _total_cov) =
                utils_frags::get_errors_cov_from_frags(set, left_snp_pos, right_snp_pos);

            let hap_q = hapqs[i];

            for i in left_snp_pos..right_snp_pos + 1 {
                snp_covered_count[(i - 1) as usize] += 1.;
                coverage_count[(i - 1) as usize] += cov
            }
            if hap_q > 0 {
                for i in left_snp_pos..right_snp_pos + 1 {
                    snp_covered_count_g0[(i - 1) as usize] += 1.;
                    coverage_count_g0[(i - 1) as usize] += cov
                }
            }

            hapQ_scores.insert(i, hap_q);

            write!(
                vartig_file,
                ">HAP{}.{}\tCONTIG:{}\tSNPRANGE:{}-{}\tBASERANGE:{}-{}\tCOV:{:.3}\tERR:{:.3}\tHAPQ:{}\tREL_ERR:{:.3}\n",
                i,
                out_bam_part_dir,
                contig,
                left_snp_pos,
                right_snp_pos,
                left_gn_pos + 1,
                right_gn_pos + 1,
                cov,
                err,
                hap_q,
                rel_err[i]
            )
            .unwrap();
            let vec_of_alleles = write_fragset_haplotypes(
                set,
                &format!("{}", i),
                &out_bam_part_dir,
                &snp_pos_to_genome_pos,
                false,
                left_snp_pos,
                right_snp_pos,
            );
            write!(
                vartig_file,
                "{}\n",
                std::str::from_utf8(
                    &vec_of_alleles
                        .into_iter()
                        .map(|x| x + 48)
                        .collect::<Vec<u8>>()
                )
                .unwrap()
            )
            .unwrap();
        }
    }

    let mut ploidy_file = OpenOptions::new()
        .write(true)
        .truncate(true)
        .create(true)
        .open(ploidy_file)
        .unwrap();

    let mut top_ploidy_file = OpenOptions::new()
        .write(true)
        .append(true)
        .create(true)
        .open(top_ploidy_file)
        .unwrap();

    let num_nonzero = snp_covered_count
        .iter()
        .filter(|x| **x > 0.)
        .collect::<Vec<_>>()
        .len();

    let num_nonzero_g0 = snp_covered_count_g0
        .iter()
        .filter(|x| **x > 0.)
        .collect::<Vec<_>>()
        .len();

    let avg_local_ploidy = snp_covered_count.iter().sum::<f64>() / num_nonzero as f64;
    let avg_local_ploidy_g0 = snp_covered_count_g0.iter().sum::<f64>() / num_nonzero_g0 as f64;
    let avg_global_ploidy = snp_covered_count.iter().sum::<f64>() / snp_covered_count.len() as f64;
    let avg_global_ploidy_g0 =
        snp_covered_count_g0.iter().sum::<f64>() / snp_covered_count_g0.len() as f64;
    //let avg_global_ploidy = total_bases_covered /
    let rough_cvg = coverage_count.iter().sum::<f64>() / num_nonzero as f64;
    write!(
        ploidy_file,
        "contig\taverage_local_ploidy\taverage_global_ploidy\tapproximate_coverage_ignoring_indels\ttotal_vartig_bases_covered\taverage_local_ploidy_min1hapq\taverage_global_ploidy_min1hapq\tavg_err\n",
    )
    .unwrap();

    write!(
        ploidy_file,
        "{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{:.3}\t{:.3}\t{:.3}\n",
        contig,
        avg_local_ploidy,
        avg_global_ploidy,
        rough_cvg,
        total_bases_covered,
        avg_local_ploidy_g0,
        avg_global_ploidy_g0,
        avg_err
    )
    .unwrap();

    write!(
        top_ploidy_file,
        "{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{:.3}\t{:.3}\t{:.3}\n",
        contig,
        avg_local_ploidy,
        avg_global_ploidy,
        rough_cvg,
        total_bases_covered,
        avg_local_ploidy_g0,
        avg_global_ploidy_g0,
        avg_err
    )
    .unwrap();

    return hapQ_scores;
}

pub fn write_all_parts_file(
    part: &Vec<FxHashSet<&Frag>>,
    contig: &str,
    snp_range_parts_vec: &Vec<(SnpPosition, SnpPosition)>,
    out_bam_part_dir: &String,
    prefix: &String,
    snp_pos_to_genome_pos: &Vec<usize>,
    hapqs: &Vec<u8>,
    rel_err: &Vec<f64>,
) {
    fs::create_dir_all(&out_bam_part_dir).unwrap();
    let part_path = &format!("{}/{}.haplosets", out_bam_part_dir, prefix);
    let file = File::create(part_path).expect("Can't create file");

    let mut file = LineWriter::new(file);
    let mut total_cov_all = 0.;
    let mut total_err_all = 0.;

    for (i, set) in part.iter().enumerate() {
        if set.is_empty() {
            continue;
        }

        //Populate all_part.txt file
        let mut vec_part: Vec<&Frag> = set.into_iter().cloned().collect();
        vec_part.sort();
        if snp_range_parts_vec.is_empty() {
            write!(file, "#{}\n", i).unwrap();
        } else {
            let left_snp_pos = snp_range_parts_vec[i].0;
            let right_snp_pos = snp_range_parts_vec[i].1;
            let (cov, err, total_err, total_cov) =
                utils_frags::get_errors_cov_from_frags(set, left_snp_pos, right_snp_pos);
            write!(
                file,
                ">HAP{}.{}\tCONTIG:{}\tSNPRANGE:{}-{}\tBASERANGE:{}-{}\tCOV:{:.3}\tERR:{:.3}\tHAPQ:{}\tREL_ERR:{:.3}\n",
                //1-indexed snp poses are output... this is annoying
                i,
                out_bam_part_dir,
                contig,
                left_snp_pos,
                right_snp_pos,
                snp_pos_to_genome_pos[(left_snp_pos - 1) as usize] + 1,
                snp_pos_to_genome_pos[(right_snp_pos - 1) as usize] + 1,
                cov,
                err,
                hapqs[i],
                rel_err[i],
            )
            .unwrap();
            total_cov_all += total_cov;
            total_err_all += total_err;
        }

        for frag in vec_part.iter() {
            //I think this was done because there
            //can be a lot of short reads. I think we should still
            //output it, though.
            write!(
                file,
                "{}\t{}\t{}\n",
                frag.id.clone(),
                frag.first_position,
                frag.last_position
            )
            .unwrap();
        }
    }
    if !snp_range_parts_vec.is_empty() {
        log::info!(
            "Final SNP error rate for all haplogroups is {}",
            total_err_all / total_cov_all
        );
    }
}

//Convert a fragment which stores sequences in a dictionary format to a block format which makes
//writing to frag files easier.
fn convert_dict_to_block(frag: Frag) -> (Vec<SnpPosition>, Vec<Vec<Genotype>>, Vec<u8>) {
    let d = frag.seq_dict;
    let vec_d: BTreeMap<SnpPosition, Genotype> = d.into_iter().collect();
    let vec_q: BTreeMap<SnpPosition, u8> = frag.qual_dict.into_iter().collect();
    let mut prev_pos = 0;
    let mut block_start_pos = Vec::new();
    let mut blocks = Vec::new();
    let mut block = Vec::new();
    let mut qual_block = Vec::new();

    for (pos, var) in &vec_d {
        if prev_pos == 0 {
            prev_pos = *pos;
            block.push(*var);
            block_start_pos.push(*pos);
        } else if pos - prev_pos > 1 {
            blocks.push(block);
            block = vec![*var];
            block_start_pos.push(*pos);
            prev_pos = *pos;
        } else if pos - prev_pos == 1 {
            block.push(*var);
            prev_pos = *pos;
        }
    }

    for (_pos, q) in &vec_q {
        qual_block.push(*q);
    }

    blocks.push(block);
    (block_start_pos, blocks, qual_block)
}

pub fn write_alignment_as_vartig(
    frags: &Vec<Frag>,
    in_file: &str,
    contig: &str,
    snp_pos_to_genome_pos: &Vec<GnPosition>,
    left_snp_pos: SnpPosition,
    right_snp_pos: SnpPosition,
    out: &str
) {
    let set_frag = frags.iter().collect();
    let hap_map = utils_frags::set_to_seq_dict(&set_frag, false);
    let emptydict = FxHashMap::default();
    let mut vec_of_alleles = vec![];
    let rightmost_base = snp_pos_to_genome_pos[(right_snp_pos-1) as usize];
    let leftmost_base = snp_pos_to_genome_pos[(left_snp_pos-1) as usize];
    for pos in left_snp_pos..right_snp_pos + 1 {
        let allele_map = hap_map.get(&pos).unwrap_or(&emptydict);
        if *allele_map == emptydict {
            //This prints ?
            vec_of_alleles.push(15);
        } else {
            let best_allele = allele_map.iter().max_by_key(|entry| entry.1).unwrap().0;
            vec_of_alleles.push(*best_allele as u8);
        }
    }
    let hap_header = format!(">HAP{}\tCONTIG:{}\tSNPRANGE:{}-{}\tBASERANGE:{}-{}\n", in_file, contig, left_snp_pos, right_snp_pos,  leftmost_base, rightmost_base);
    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out)
        .unwrap();
    write!(file, "{}", hap_header).unwrap();
    write!(
        file,
        "{}\n",
        std::str::from_utf8(
            &vec_of_alleles
                .into_iter()
                .map(|x| x + 48)
                .collect::<Vec<u8>>()
        )
        .unwrap()
    )
    .unwrap();

}
