use crate::monomers::DnaNucleotide;
use crate::monomers::Monomer;
use crate::monomers::Nucleotide;
use crate::monomers::RnaNucleotide;
use core::hash::Hash;
use std::collections::HashMap;

//-------------------------------------- DNA -----------------------------------------//

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Polymer<M: Monomer> {
    monomers: Vec<M>,
}

pub type Dna = Polymer<DnaNucleotide>;

pub type Rna = Polymer<RnaNucleotide>;

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
        let mut distance = 0;
        for monomers in self.monomers.iter().zip(other.monomers.iter()) {
            let (self_monomer, other_monomer) = monomers;
            if self_monomer != other_monomer {
                distance += 1;
            }
        }
        distance
    }
    pub fn motif_locations(&self, motif: &Self) -> Vec<usize> {
        let mut locations = Vec::new();
        for (i, subpolymer) in self.monomers.windows(motif.monomers.len()).enumerate() {
            if subpolymer == motif.monomers {
                let location = i + 1;
                locations.push(location);
            }
        }
        locations
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
        for size in (1..shortest.monomers.len() + 1).rev() {
            for motif in shortest.monomers.windows(size) {
                let mut shared = true;
                for (i, polymer) in polymers.iter().enumerate() {
                    if i == shortest_i {
                        continue;
                    }
                    if !polymer.contains_motif(&Polymer {
                        monomers: motif.to_vec(),
                    }) {
                        shared = false;
                        break;
                    }
                }
                if shared {
                    return Some(Polymer {
                        monomers: motif.to_vec(),
                    });
                }
            }
        }
        None
    }
    pub fn contains_motif(&self, motif: &Self) -> bool {
        for subpolymer in self.monomers.windows(motif.monomers.len()) {
            if subpolymer == motif.monomers {
                return true;
            }
        }
        false
    }
}

impl<N: Nucleotide> Polymer<N> {
    pub fn gc_content(&self) -> f32 {
        let count = self.monomers.iter().filter(|n| n.is_gc()).count() as f32;
        let length = self.monomers.len() as f32;
        count * 100.0 / length
    }
    pub fn reverse_complement(&self) -> Polymer<N> {
        let monomers = self.monomers.iter().rev().map(|n| n.complement()).collect();
        Polymer { monomers }
    }
}

impl Dna {
    pub fn transcribe(&self) -> Rna {
        let monomers = self.monomers.iter().map(|n| n.transcribe()).collect();
        Polymer { monomers }
    }
}

impl Rna {
    pub fn untranscribe(&self) -> Dna {
        let monomers = self.monomers.iter().map(|n| n.untranscribe()).collect();
        Polymer { monomers }
    }
}
