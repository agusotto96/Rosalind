use crate::monomers::AminoAcid;
use crate::monomers::Codon;
use crate::monomers::DnaNucleotide;
use crate::monomers::Monomer;
use crate::monomers::Nucleotide;
use crate::monomers::RnaNucleotide;
use core::hash::Hash;
use std::collections::HashMap;

//------------------------------------ POLYMER ---------------------------------------//

#[derive(PartialEq, Eq, Hash, Clone, Ord, PartialOrd, Debug)]
pub struct Polymer<M: Monomer> {
    monomers: Vec<M>,
}

impl<M: Monomer> Polymer<M> {
    pub fn new(symbols: &str, monomer: fn(char) -> Option<M>) -> Option<Self> {
        if symbols.is_empty() {
            None
        } else {
            symbols
                .chars()
                .map(monomer)
                .collect::<Option<Vec<M>>>()
                .map(|monomers| Polymer { monomers })
        }
    }
    pub fn monomer_count(&self) -> HashMap<M, usize> {
        let mut count = HashMap::new();
        for monomer in &self.monomers {
            let prev = count.get(monomer);
            let current = prev.unwrap_or(&0) + 1;
            count.insert(*monomer, current);
        }
        count
    }
    pub fn hamming_distance(&self, other: &Self) -> usize {
        self.monomers
            .iter()
            .zip(other.monomers.iter())
            .filter(|p| p.0 != p.1)
            .count()
    }
    pub fn motif_locations(&self, motif: &Self) -> Vec<usize> {
        self.monomers
            .windows(motif.monomers.len())
            .enumerate()
            .filter(|e| e.1 == motif.monomers)
            .map(|e| e.0 + 1)
            .collect()
    }
    pub fn profile(polymers: &[Polymer<M>]) -> HashMap<M, Vec<usize>> {
        let mut profile = HashMap::new();
        let length = polymers[0].monomers.len();
        for i in 0..length {
            for polymer in polymers {
                let monomer = polymer.monomers[i];
                profile.entry(monomer).or_insert_with(|| vec![0; length])[i] += 1;
            }
        }
        profile
    }
    pub fn consensus(polymers: &[Polymer<M>]) -> Polymer<M> {
        let mut count = HashMap::new();
        let mut consensus = Vec::new();
        let length = polymers[0].monomers.len();
        for i in 0..length {
            for polymer in polymers {
                let monomer = polymer.monomers[i];
                *count.entry(monomer).or_insert_with(|| 0) += 1;
            }
            let monomer = *count.iter().max_by_key(|e| e.1).unwrap().0;
            consensus.push(monomer);
            count.clear();
        }
        Polymer {
            monomers: consensus,
        }
    }
    pub fn shared_motif(polymers: &[Polymer<M>]) -> Option<Polymer<M>> {
        let (shortest_i, shortest) = polymers
            .iter()
            .enumerate()
            .min_by_key(|p| p.1.monomers.len())
            .unwrap();
        for size in (1..=shortest.monomers.len()).rev() {
            for motif in shortest.monomers.windows(size) {
                let mut shared = true;
                for (i, polymer) in polymers.iter().enumerate() {
                    if i == shortest_i {
                        continue;
                    }
                    if !polymer.monomers.windows(motif.len()).any(|s| s == motif) {
                        shared = false;
                        break;
                    }
                }
                if shared {
                    let monomers = motif.to_vec();
                    let polymer = Polymer { monomers };
                    return Some(polymer);
                }
            }
        }
        None
    }
}

//---------------------------------- NUCLEIC ACID ------------------------------------//

impl<N: Nucleotide> Polymer<N> {
    pub fn gc_content(&self) -> f64 {
        let count = self.monomers.iter().filter(|n| n.is_gc()).count() as f64;
        let length = self.monomers.len() as f64;
        count * 100.0 / length
    }
    pub fn reverse_complement(&self) -> Polymer<N> {
        let monomers = self
            .monomers
            .iter()
            .rev()
            .map(Nucleotide::complement)
            .collect();
        Polymer { monomers }
    }
    pub fn transition_transversion_ratio(&self, other: &Self) -> f64 {
        let mut transition = 0.0;
        let mut transversion = 0.0;
        for nucleotides in self.monomers.iter().zip(other.monomers.iter()) {
            let (self_nucleotide, other_nucleotide) = nucleotides;
            if self_nucleotide != other_nucleotide {
                if self_nucleotide.is_purine() == other_nucleotide.is_purine() {
                    transition += 1.0;
                } else {
                    transversion += 1.0;
                }
            }
        }
        transition / transversion
    }
}

//-------------------------------------- DNA -----------------------------------------//

pub type Dna = Polymer<DnaNucleotide>;

impl Dna {
    pub fn transcribe(&self) -> Rna {
        let monomers = self
            .monomers
            .iter()
            .map(DnaNucleotide::transcribe)
            .collect();
        Polymer { monomers }
    }
}

//-------------------------------------- RNA -----------------------------------------//

pub type Rna = Polymer<RnaNucleotide>;

impl Rna {
    pub fn untranscribe(&self) -> Dna {
        let monomers = self
            .monomers
            .iter()
            .map(RnaNucleotide::untranscribe)
            .collect();
        Polymer { monomers }
    }
    pub fn translate(&self) -> Vec<Protein> {
        let mut candidates = Vec::new();
        let mut translations = Vec::new();
        let chunks = self.monomers.chunks_exact(3);
        for chunk in chunks {
            let codon = Codon(chunk[0], chunk[1], chunk[2]);
            match codon.aminoacid() {
                Some(aminoacid) => {
                    if aminoacid.is_start() {
                        candidates.push(Vec::new());
                    }
                    for canidate in &mut candidates {
                        canidate.push(aminoacid);
                    }
                }
                None => {
                    translations.append(&mut candidates);
                }
            }
        }
        translations
            .iter()
            .map(|t| t.to_vec())
            .map(|monomers| Polymer { monomers })
            .collect()
    }
    pub fn reading_frames(&self) -> Vec<Rna> {
        let mut frames = Vec::new();
        for i in 0..=2 {
            if self.monomers.len() > i {
                let monomers = self.monomers[i..]
                    .chunks_exact(3)
                    .flat_map(|c| c.to_vec())
                    .collect::<Vec<RnaNucleotide>>();
                let frame = Polymer { monomers };
                frames.push(frame);
            }
        }
        frames
    }
}

//------------------------------------ PROTEIN ---------------------------------------//

pub type Protein = Polymer<AminoAcid>;

impl Protein {
    pub fn mass(&self) -> f64 {
        self.monomers.iter().map(AminoAcid::mass).sum()
    }
}
