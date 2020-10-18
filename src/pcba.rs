use std::{ fs, io };
use io::prelude::*;

fn main() -> io::Result<()> {
    let _ = fs::create_dir("./results");
    let in_path = "./data/pcba.csv";
    let out_path = "./data/pcba.smi";
    let in_file = fs::File::open(in_path)?;
    let out_file = fs::File::create(out_path)?;
    let reader = io::BufReader::new(in_file);
    let mut writer = io::LineWriter::new(out_file);
    let mut lines = reader.lines();

    let _header = lines.next().unwrap();

    for (id, result) in lines.enumerate() {
        let line = result?;
        let parts = line.split(",").collect::<Vec<_>>();
        let smiles = parts.last().expect("no smiles");

        writeln!(writer, "{} {}", smiles, id);
    }

    Ok(())
}