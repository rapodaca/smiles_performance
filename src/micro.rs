use std::{ fs, io, time };
use io::prelude::*;
use chemcore::daylight;

fn main() -> io::Result<()> {
    let in_path = "./data/rdkit_2017.03.3.smi";
    let in_file = fs::File::open(in_path)?;
    let reader = io::BufReader::new(in_file);
    let mut inputs = Vec::new();

    for result in reader.lines() {
        let line = result?;
        let parts = line.split_whitespace().collect::<Vec<_>>();

        if parts.len() == 1 {
            continue;
        }

        let smiles = parts.first().expect("no smiles").to_string();

        inputs.push(smiles);
    }

    let start = time::Instant::now();

    for smiles in inputs {
        let _ = daylight::read(&smiles);
    }

    println!("elapsed: {}", start.elapsed().as_secs_f32());

    Ok(())
}