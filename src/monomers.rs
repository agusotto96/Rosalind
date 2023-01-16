use core::fmt::Debug;
use core::hash::Hash;

//------------------------------------ MONOMER ---------------------------------------//

pub trait Monomer: Sized + PartialEq + Eq + Hash + Copy + Clone + Ord + PartialOrd + Debug {
    fn new(symbol: char) -> Option<Self>;
}

//----------------------------------- NUCLEOTIDE -------------------------------------//

pub trait Nucleotide: Monomer {
    fn complement(&self) -> Self;
    fn is_gc(&self) -> bool;
    fn is_purine(&self) -> bool;
    fn is_pyrimidine(&self) -> bool;
}

//--------------------------------- DNA NUCLEOTIDE -----------------------------------//

#[derive(PartialEq, Eq, Hash, Copy, Clone, Ord, PartialOrd, Debug)]
pub enum DnaNucleotide {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
}

impl Monomer for DnaNucleotide {
    fn new(symbol: char) -> Option<DnaNucleotide> {
        match symbol {
            'A' => Some(DnaNucleotide::Adenine),
            'C' => Some(DnaNucleotide::Cytosine),
            'G' => Some(DnaNucleotide::Guanine),
            'T' => Some(DnaNucleotide::Thymine),
            _ => None,
        }
    }
}

impl Nucleotide for DnaNucleotide {
    fn complement(&self) -> DnaNucleotide {
        match &self {
            DnaNucleotide::Adenine => DnaNucleotide::Thymine,
            DnaNucleotide::Cytosine => DnaNucleotide::Guanine,
            DnaNucleotide::Guanine => DnaNucleotide::Cytosine,
            DnaNucleotide::Thymine => DnaNucleotide::Adenine,
        }
    }
    fn is_gc(&self) -> bool {
        matches!(&self, DnaNucleotide::Guanine | DnaNucleotide::Cytosine)
    }
    fn is_purine(&self) -> bool {
        matches!(&self, DnaNucleotide::Adenine | DnaNucleotide::Guanine)
    }
    fn is_pyrimidine(&self) -> bool {
        matches!(&self, DnaNucleotide::Cytosine | DnaNucleotide::Thymine)
    }
}

impl DnaNucleotide {
    pub fn transcribe(&self) -> RnaNucleotide {
        match &self {
            DnaNucleotide::Adenine => RnaNucleotide::Adenine,
            DnaNucleotide::Cytosine => RnaNucleotide::Cytosine,
            DnaNucleotide::Guanine => RnaNucleotide::Guanine,
            DnaNucleotide::Thymine => RnaNucleotide::Uracil,
        }
    }
}

//--------------------------------- RNA NUCLEOTIDE -----------------------------------//

#[derive(PartialEq, Eq, Hash, Copy, Clone, Ord, PartialOrd, Debug)]
pub enum RnaNucleotide {
    Adenine,
    Cytosine,
    Guanine,
    Uracil,
}

impl Monomer for RnaNucleotide {
    fn new(symbol: char) -> Option<RnaNucleotide> {
        match symbol {
            'A' => Some(RnaNucleotide::Adenine),
            'C' => Some(RnaNucleotide::Cytosine),
            'G' => Some(RnaNucleotide::Guanine),
            'U' => Some(RnaNucleotide::Uracil),
            _ => None,
        }
    }
}

impl Nucleotide for RnaNucleotide {
    fn complement(&self) -> RnaNucleotide {
        match &self {
            RnaNucleotide::Adenine => RnaNucleotide::Uracil,
            RnaNucleotide::Cytosine => RnaNucleotide::Guanine,
            RnaNucleotide::Guanine => RnaNucleotide::Cytosine,
            RnaNucleotide::Uracil => RnaNucleotide::Adenine,
        }
    }
    fn is_gc(&self) -> bool {
        matches!(&self, RnaNucleotide::Guanine | RnaNucleotide::Cytosine)
    }
    fn is_purine(&self) -> bool {
        matches!(&self, RnaNucleotide::Adenine | RnaNucleotide::Guanine)
    }
    fn is_pyrimidine(&self) -> bool {
        matches!(&self, RnaNucleotide::Cytosine | RnaNucleotide::Uracil)
    }
}

impl RnaNucleotide {
    pub fn untranscribe(&self) -> DnaNucleotide {
        match &self {
            RnaNucleotide::Adenine => DnaNucleotide::Adenine,
            RnaNucleotide::Cytosine => DnaNucleotide::Cytosine,
            RnaNucleotide::Guanine => DnaNucleotide::Guanine,
            RnaNucleotide::Uracil => DnaNucleotide::Thymine,
        }
    }
}

//----------------------------------- AMINO ACID -------------------------------------//

#[derive(PartialEq, Eq, Hash, Copy, Clone, Ord, PartialOrd, Debug)]
pub enum AminoAcid {
    Alanine,
    Cysteine,
    AsparticAcid,
    GlutamicAcid,
    Phenylalanine,
    Glycine,
    Histidine,
    Isoleucine,
    Lysine,
    Leucine,
    Methionine,
    Asparagine,
    Proline,
    Glutamine,
    Arginine,
    Serine,
    Threonine,
    Valine,
    Tryptophan,
    Tyrosine,
}

impl Monomer for AminoAcid {
    fn new(symbol: char) -> Option<AminoAcid> {
        match symbol {
            'A' => Some(AminoAcid::Alanine),
            'C' => Some(AminoAcid::Cysteine),
            'D' => Some(AminoAcid::AsparticAcid),
            'E' => Some(AminoAcid::GlutamicAcid),
            'F' => Some(AminoAcid::Phenylalanine),
            'G' => Some(AminoAcid::Glycine),
            'H' => Some(AminoAcid::Histidine),
            'I' => Some(AminoAcid::Isoleucine),
            'K' => Some(AminoAcid::Lysine),
            'L' => Some(AminoAcid::Leucine),
            'M' => Some(AminoAcid::Methionine),
            'N' => Some(AminoAcid::Asparagine),
            'P' => Some(AminoAcid::Proline),
            'Q' => Some(AminoAcid::Glutamine),
            'R' => Some(AminoAcid::Arginine),
            'S' => Some(AminoAcid::Serine),
            'T' => Some(AminoAcid::Threonine),
            'V' => Some(AminoAcid::Valine),
            'W' => Some(AminoAcid::Tryptophan),
            'Y' => Some(AminoAcid::Tyrosine),
            _ => None,
        }
    }
}

impl AminoAcid {
    pub fn mass(&self) -> f64 {
        match &self {
            AminoAcid::Alanine => 71.037_11,
            AminoAcid::Cysteine => 103.009_19,
            AminoAcid::AsparticAcid => 115.026_94,
            AminoAcid::GlutamicAcid => 129.042_59,
            AminoAcid::Phenylalanine => 147.068_41,
            AminoAcid::Glycine => 57.021_46,
            AminoAcid::Histidine => 137.058_91,
            AminoAcid::Isoleucine => 113.084_06,
            AminoAcid::Lysine => 128.094_96,
            AminoAcid::Leucine => 113.084_06,
            AminoAcid::Methionine => 131.040_49,
            AminoAcid::Asparagine => 114.042_93,
            AminoAcid::Proline => 97.052_76,
            AminoAcid::Glutamine => 128.058_58,
            AminoAcid::Arginine => 156.101_11,
            AminoAcid::Serine => 87.032_03,
            AminoAcid::Threonine => 101.047_68,
            AminoAcid::Valine => 99.068_41,
            AminoAcid::Tryptophan => 186.079_31,
            AminoAcid::Tyrosine => 163.063_33,
        }
    }
    pub fn is_start(&self) -> bool {
        matches!(&self, AminoAcid::Methionine)
    }
}

//------------------------------------- CODON ----------------------------------------//

pub struct Codon(pub RnaNucleotide, pub RnaNucleotide, pub RnaNucleotide);

impl Codon {
    pub fn aminoacid(&self) -> Option<AminoAcid> {
        match &self {
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Uracil, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Phenylalanine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Uracil, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Leucine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Uracil, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Isoleucine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Uracil, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Valine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Uracil, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Phenylalanine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Uracil, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Leucine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Uracil, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Isoleucine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Uracil, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Valine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Uracil, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Leucine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Uracil, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Leucine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Uracil, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Isoleucine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Uracil, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Valine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Uracil, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Leucine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Uracil, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Leucine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Uracil, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Valine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Cytosine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Serine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Cytosine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Proline)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Cytosine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Threonine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Cytosine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Alanine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Cytosine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Serine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Cytosine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Proline)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Cytosine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Threonine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Cytosine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Alanine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Cytosine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Serine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Cytosine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Proline)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Cytosine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Threonine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Cytosine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Alanine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Cytosine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Serine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Cytosine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Proline)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Cytosine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Threonine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Cytosine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Alanine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Adenine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Tyrosine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Adenine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Histidine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Adenine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Asparagine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Adenine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::AsparticAcid)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Adenine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Tyrosine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Adenine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Histidine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Adenine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Asparagine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Adenine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::AsparticAcid)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Adenine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Glutamine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Adenine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Lysine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Adenine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::GlutamicAcid)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Adenine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Glutamine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Adenine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Lysine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Adenine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::GlutamicAcid)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Guanine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Cysteine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Guanine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Arginine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Guanine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Serine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Guanine, RnaNucleotide::Uracil) => {
                Some(AminoAcid::Glycine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Guanine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Cysteine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Guanine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Arginine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Guanine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Serine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Guanine, RnaNucleotide::Cytosine) => {
                Some(AminoAcid::Glycine)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Guanine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Arginine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Guanine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Arginine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Guanine, RnaNucleotide::Adenine) => {
                Some(AminoAcid::Glycine)
            }
            Codon(RnaNucleotide::Uracil, RnaNucleotide::Guanine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Tryptophan)
            }
            Codon(RnaNucleotide::Cytosine, RnaNucleotide::Guanine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Arginine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Guanine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Arginine)
            }
            Codon(RnaNucleotide::Guanine, RnaNucleotide::Guanine, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Glycine)
            }
            Codon(RnaNucleotide::Adenine, RnaNucleotide::Uracil, RnaNucleotide::Guanine) => {
                Some(AminoAcid::Methionine)
            }
            _ => None,
        }
    }
}
