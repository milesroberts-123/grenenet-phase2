use paraseq::{fastx, ProcessError};
use paraseq::prelude::*;
use std::env::args;
use std::collections::HashMap;

#[derive(Clone, Default)]
struct MyProcessor {
// Your processing state here
}

struct SharedData {
    // ... your data
    let mut kmer_counts: HashMap<String, usize> = HashMap::new();
    let mut qmer_counts: HashMap<String, u32> = HashMap::new();
    let k = 31;
    let q = 35;
}

fn total_qscore(quality: Vec<u8>) -> u32{
    //let result: u32 = quality.iter().map(|&x| x as u32).min();
    //let product: u16 = quality.iter().sum();
    let min_value: u32 = quality.iter().map(|&x| x as u32).min().unwrap();
    return min_value;
}

fn reverse_complement(dna_sequence: &str) -> String {
    let mut complemented_sequence = String::new();

    // Step 1: Complement each base
    for base in dna_sequence.chars() {
        let complemented_base = match base {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            // Handle other characters (e.g., 'N' for unknown bases) if needed
            _ => base, // Keep non-DNA bases as they are
        };
        complemented_sequence.push(complemented_base);
    }

    // Step 2: Reverse the complemented sequence
    complemented_sequence.chars().rev().collect()
}

fn get_can(forward: String, reverse: String) -> String {
    if forward <= reverse {
        return forward;
    } else {
        return reverse;
    }
}

impl ParallelProcessor for MyProcessor {
    fn process_record<R: Record>(&mut self, record: R) -> Result<(), ProcessError> {
        // Process record in parallel
        let sequence = record.seq_str();
        let quality = record.qual();

        if k == 0 || k > sequence.len() {
            continue;
            //return (kmer_counts, qmer_counts); // Handle invalid k or sequence length
        }

        for i in 0..=(sequence.len() - k) {

            let qmer = &quality[i..i + k];
            let min_score = total_qscore(qmer.to_vec());

            if min_score < q {
                continue;
            }

            let kmer = sequence[i..i + k].to_string();

            // find cannonical k-mer
            let rev_kmer = reverse_complement(&kmer);

            let can_kmer = get_can(kmer, rev_kmer);

            *kmer_counts.entry(can_kmer.clone()).or_insert(0) += 1;
            *qmer_counts.entry(can_kmer).or_insert(0) += min_score;
        }


        Ok(())
    }
}

fn main() -> Result<(), ProcessError> {
    //let path = "./data/sample.fastq";
    let path: String = args().nth(1).unwrap();
    let reader = fastx::Reader::from_path(path)?;
    let processor = MyProcessor::default();
    let num_threads = 8;

    reader.process_parallel(processor, num_threads)?;
    Ok(())

}
