use rosalind::monomers::AminoAcid;
use rosalind::monomers::Monomer;
use rosalind::monomers::RnaNucleotide;
use rosalind::polymers::Polymer;
use rosalind::polymers::Protein;
use rosalind::polymers::Rna;

#[test]
fn translate() {
    let rna = new_rna("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");
    let mut actual = rna.translate();
    let mut expected = [
        new_protein("MAMAPRTEINSTRING"),
        new_protein("MAPRTEINSTRING"),
    ];
    actual.sort();
    expected.sort();
    assert_eq!(actual, expected);
}

fn new_rna(symbols: &str) -> Rna {
    Polymer::new(symbols, RnaNucleotide::new).unwrap()
}

fn new_protein(symbols: &str) -> Protein {
    Polymer::new(symbols, AminoAcid::new).unwrap()
}
