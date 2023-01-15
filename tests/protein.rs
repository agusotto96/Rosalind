use rosalind::monomers::AminoAcid;
use rosalind::monomers::Monomer;
use rosalind::polymers::Polymer;
use rosalind::polymers::Protein;

#[test]
fn mass() {
    let protein = new_protein("SKADYEK");
    let actual = protein.mass();
    let expected = 821.392;
    assert!((actual - expected).abs() <= 0.0001);
}

fn new_protein(symbols: &str) -> Protein {
    Polymer::new(symbols, AminoAcid::new).unwrap()
}
