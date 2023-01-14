use core::fmt::Debug;
use core::hash::Hash;

//------------------------------------ MONOMER ---------------------------------------//

pub trait Monomer: Sized + PartialEq + Eq + Hash + Copy + Clone + Debug {
    fn new(symbol: char) -> Option<Self>;
}

//----------------------------------- NUCLEOTIDE -------------------------------------//

pub trait Nucleotide: Monomer {
    fn complement(&self) -> Self;
    fn is_gc(&self) -> bool;
}

//--------------------------------- DNA NUCLEOTIDE -----------------------------------//

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
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
        matches!(self, DnaNucleotide::Guanine | DnaNucleotide::Cytosine)
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

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
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
        matches!(self, RnaNucleotide::Guanine | RnaNucleotide::Cytosine)
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
