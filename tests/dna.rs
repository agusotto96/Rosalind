use rosalind::monomers::DnaNucleotide;
use rosalind::monomers::Monomer;
use rosalind::monomers::RnaNucleotide;
use rosalind::polymers::Dna;
use rosalind::polymers::Polymer;
use rosalind::polymers::Rna;

use std::collections::HashMap;

#[test]
fn monomer_count() {
    let dna = new_dna("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC");
    let actual = dna.monomer_count();
    let expected = HashMap::from([
        (DnaNucleotide::Adenine, 20),
        (DnaNucleotide::Cytosine, 12),
        (DnaNucleotide::Guanine, 17),
        (DnaNucleotide::Thymine, 21),
    ]);
    assert_eq!(actual, expected);
}

#[test]
fn transcribe() {
    let dna = new_dna("GATGGAACTTGACTACGTAAATT");
    let actual = dna.transcribe();
    let expected = new_rna("GAUGGAACUUGACUACGUAAAUU");
    assert_eq!(actual, expected);
}

#[test]
fn reverse_complement() {
    let dna = new_dna("AAAACCCGGT");
    let actual = dna.reverse_complement();
    let expected = new_dna("ACCGGGTTTT");
    assert_eq!(actual, expected);
}

#[test]
fn gc_content() {
    let dna = new_dna(
        "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT",
    );
    let actual = dna.gc_content();
    let expected = 60.91954;
    assert_eq!(actual, expected);
}

#[test]
fn hamming_distance() {
    let dna = new_dna("GAGCCTACTAACGGGAT");
    let actual = dna.hamming_distance(&new_dna("CATCGTAATGACGGCCT"));
    let expected = 7;
    assert_eq!(actual, expected);
}

#[test]
fn consensus() {
    let dnas = vec![
        new_dna("ATCCAGCT"),
        new_dna("GGGCAACT"),
        new_dna("ATGGATCT"),
        new_dna("AAGCAACC"),
        new_dna("TTGGAACT"),
        new_dna("ATGCCATT"),
        new_dna("ATGGCACT"),
    ];
    let actual = Dna::consensus(&dnas);
    let expected = new_dna("ATGCAACT");
    assert_eq!(actual, expected);
}

#[test]
fn motif_locations() {
    let dna = new_dna("GATATATGCATATACTT");
    let actual = dna.motif_locations(&new_dna("ATAT"));
    let expected = [2, 4, 10];
    assert_eq!(actual, expected);
}

fn new_dna(symbols: &str) -> Dna {
    Polymer::new(symbols, DnaNucleotide::new).unwrap()
}

fn new_rna(symbols: &str) -> Rna {
    Polymer::new(symbols, RnaNucleotide::new).unwrap()
}
