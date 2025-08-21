use std::env::args;
use kseq::parse_path;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

fn qual_str_to_int(quality_string: &str) -> Vec<u8>{
        let phred_scores: Vec<u8> = quality_string
        .chars()
        .map(|c| c as u8 - 33) // Subtract 33 for Phred+33 encoding
        .collect();

        return phred_scores;
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

//fn count_kmers(sequence: &str, quality: Vec<u8>, k: usize) -> (HashMap<String, usize>,  HashMap<String, u32>) {
//    let mut kmer_counts: HashMap<String, usize> = HashMap::new();
//    let mut qmer_counts: HashMap<String, u32> = HashMap::new();

  //  if k == 0 || k > sequence.len() {
    //    return (kmer_counts, qmer_counts); // Handle invalid k or sequence length
    //}

    //for i in 0..=(sequence.len() - k) {

      //  let qmer = &quality[i..i + k];
        //let min_score = total_qscore(qmer.to_vec());

        //if min_score < 35{
          //  continue;
        //}

        //let kmer = sequence[i..i + k].to_string();
        //println!("{}", kmer);
        

        //println!("{:?}", total_qscore(qmer.to_vec()));
        
        //*kmer_counts.entry(kmer.clone()).or_insert(0) += 1;
        //*qmer_counts.entry(kmer).or_insert(0) += min_score;
    //}

    //(kmer_counts, qmer_counts)
//}

// write vector to file
fn write_file(output: &Vec<String>, prefix: &String, suffix: &str) {
    println!("{}", format!("Saving results to {}_{}...", prefix, suffix));
    let mut file = File::create(format!("{}_{}", prefix, suffix)).expect("Unable to create file");
    writeln!(file, "{}", output.join("\n")).expect("Unable to write to file");
}

fn main() -> Result<(), Box<dyn std::error::Error>>{
	let path: String = args().nth(1).unwrap();
	let mut records = parse_path(path).unwrap();
    let mut kmer_counts: HashMap<String, usize> = HashMap::new();
    let mut qmer_counts: HashMap<String, u32> = HashMap::new();
    let k = 31;
    let q = 35;
    let mut w = 0;
    println!("Counting k-mers...");
    // let mut records = parse_reader(File::open(path).unwrap()).unwrap();
	while let Some(record) = records.iter_record().unwrap() {
		//println!("head:{} des:{} seq:{} qual:{} len:{}", 
		//	record.head(), record.des(), record.seq(), 
		//	record.qual(), record.len());
        //println!("{:?}", qual_str_to_int(record.qual()));
        //count_kmers(record.seq(), qual_str_to_int(record.qual()), 31);
        if (w + 1) % 10000 == 0 {
            println!("Iteration {}", w + 1);
        }

        w = w + 1;

        let sequence = record.seq();
        let quality = qual_str_to_int(record.qual());

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
    }

    // Iterate over the keys of the HashMap

    println!("Writing counts...");
    // Open a file for writing, creating it if it doesn't exist, and truncating it if it does
    let mut file1 = File::create("kmers.txt")?;
    let mut file2 = File::create("qmers.txt")?;

    // Iterate over the HashMap and write each key-value pair to the file
    for (key, value) in &kmer_counts {
        writeln!(file1, "{} {}", key, value)?;
    }

    for (key, value) in &qmer_counts {
        writeln!(file2, "{} {}", key, value)?;
    }

    println!("HashMap written to output.txt successfully!");
    Ok(())
}
