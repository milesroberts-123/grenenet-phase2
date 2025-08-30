use paraseq::ParallelFastq;
use dashmap::DashMap;
use std::sync::Arc;
use std::fs::File;
use std::io::{self, BufReader};
use std::env::args;

// A helper function to generate canonical k-mers from a sequence
// using a sliding window.
fn canonical_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    let mut kmers = Vec::new();
    if seq.len() < k {
        return kmers;
    }

    for window in seq.windows(k) {
        let mut kmer = 0;
        let mut rc_kmer = 0;
        for &base in window {
            let nuc = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => continue, // Skip non-ACGT characters
            };
            kmer = (kmer << 2) | nuc;
            rc_kmer = (rc_kmer >> 2) | ((3 - nuc) << (2 * (k - 1)));
        }
        kmers.push(u64::min(kmer, rc_kmer));
    }
    kmers
}

fn main() -> io::Result<()> {
    let k = 31; // K-mer size
    //let fastq_path = "your_data.fastq"; // Replace with your input file
    let path: String = args().nth(1).unwrap();
    
    // Use a buffered reader for efficient I/O.
    let file = BufReader::new(File::open(fastq_path)?);

    // Create a shared concurrent hash map wrapped in an Arc.
    let kmer_counts = Arc::new(DashMap::<u64, usize>::new());

    // Create the parallel fastq iterator from `paraseq`.
    let iterator = ParallelFastq::new(file, None);
    
    // Process the records in parallel.
    iterator.for_each(|record_set| {
        if let Some(record_set) = record_set {
            // Get the DashMap reference
            let counts = kmer_counts.clone();
            // Process each record within the set
            for record in record_set.records() {
                if let Some(seq) = record.seq() {
                    for kmer in canonical_kmers(seq, k) {
                        *counts.entry(kmer).or_insert(0) += 1;
                    }
                }
            }
        }
    })?;

    // Print the final counts
    println!("Total unique k-mers counted: {}", kmer_counts.len());
    for entry in kmer_counts.iter() {
        println!("k-mer hash: {} -> count: {}", entry.key(), entry.value());
    }

    Ok(())
}

