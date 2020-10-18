use std::{ fs, io };
use io::prelude::*;
use chemcore::daylight;
use chemcore::molecule::{ Molecule, Error };
use gamma::graph::Graph;

fn main() -> io::Result<()> {
    let _ = fs::create_dir("./results");
    let in_path = "./data/rdkit_2017.03.3.smi";
    let out_path = "./results/chemcore-out.txt";
    // let in_path = "./data/pcba.smi";
    // let out_path = "./results/chemcore-pcba-out.txt";
    let in_file = fs::File::open(in_path)?;
    let out_file = fs::File::create(out_path)?;
    let reader = io::BufReader::new(in_file);
    let mut writer = io::LineWriter::new(out_file);

    for result in reader.lines() {
        let line = result?;
        let parts = line.split_whitespace().collect::<Vec<_>>();
        let id = parts.last().expect("no parts");

        if parts.len() == 1 {
            writeln!(&mut writer, "# {} No_input", id)?;

            continue;
        }

        let smiles = parts.first().expect("no smiles");

        match daylight::read(smiles) {
            Ok(molecule) => {
                let hcounts = molecule.nodes().iter().map(|id| {
                    molecule.hydrogens(*id).unwrap().to_string()
                }).collect::<Vec<_>>();

                writeln!(&mut writer, "{} {}", id, hcounts.join(" "))?;
            },
            Err(Error::CanNotKekulize) => {
                writeln!(writer, "# {} Kekulization_failure", id)?;
            },
            Err(Error::Hypervalent(_)) => {
                writeln!(writer, "# {} Bad_valence", id)?;
            },
            Err(error) => {
                writeln!(writer, "# {} ERROR: {:?}", id, error)?;
            }
        }
    }

    Ok(())
}