use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn three_ploidy_standard_test() {
    let mut cmd = Command::cargo_bin("flopp").unwrap();
    let assert = cmd
        .arg("-b")
        .arg("./tests/test_bams/pds_ploidy3.bam")
        .arg("-v")
        .arg("./tests/test_vcfs/pds.vcf")
        .arg("-p")
        .arg("3")
        .assert();
    assert
        .success()
        .code(0);

}

#[test]
fn three_ploidy_fragment_test(){
    let mut cmd = Command::cargo_bin("flopp").unwrap();
    let assert = cmd
        .arg("-f")
        .arg("./tests/3xploidy_frags.txt")
        .arg("-p")
        .arg("3")
        .assert();
    assert
        .success()
        .code(0);
}

#[test]
fn frag_dump_test(){
    let mut cmd = Command::cargo_bin("frag-dump").unwrap();
    let assert = cmd
        .arg("-b")
        .arg("./tests/test_bams/pds_ploidy3.bam")
        .arg("-v")
        .arg("./tests/test_vcfs/pds.vcf")
        .arg("-o")
        .arg("./tests/output_frag.txt")
        .assert();
    assert
        .success()
        .code(0);
}

#[test]
fn multiple_ref_test(){
    let mut cmd = Command::cargo_bin("frag-dump").unwrap();
    let assert = cmd
        .arg("-b")
        .arg("./tests/test_bams/sorted_merged_bam_3x.bam")
        .arg("-v")
        .arg("./tests/test_vcfs/merged_vcf.vcf")
        .arg("-o")
        .arg("./tests/output_frag.txt")
        .assert();
    assert
        .success()
        .code(0);

}

