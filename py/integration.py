import sys
import pathlib
from io import StringIO
from rdkit import Chem

in_path = "./data/rdkit_2017.03.3.smi"
out_path = "./results/rdkit-out.txt"
# in_path = "./data/pcba.smi"
# out_path = "./results/rdkit-pcba-out.txt"

pathlib.Path('./results').mkdir(exist_ok=True)
Chem.WrapLogs()

sio = sys.stderr = StringIO()

with open(in_path) as in_file, open(out_path, 'w') as out_file:
    for line in in_file:
        parts = line.strip().split()
        id = parts[-1]

        if len(parts) == 1:
            out_file.write(f'# {id} No_input\n')
            continue

        smiles = parts[0]

        mol = Chem.MolFromSmiles(smiles)
        err = sio.getvalue()
        if err:
            sio = sys.stderr = StringIO()
            if "Can't kekulize" in err:
                out_file.write(f'# {id} Kekulization_failure\n')
                continue
            elif "Explicit valence" in err:
                out_file.write(f'# {id} Bad_valence\n')
                continue
            elif "SMILES Parse Error" in err:
                out_file.write(f'# {id} SMILES_parse_error\n')
                continue
            elif "Aromatic bonds on non aromatic atom" in err:
                out_file.write(f'# {id} Aromatic_bonds_on_non_aromatic_atom\n')
                continue
            elif "non-ring" in err and "marked aromatic" in err:
                out_file.write(f'# {id} Non_ring_atom_marked_aromatic\n')
                continue
            elif "WARNING" not in err:
                out_file.write(f'# {id} ERROR: {err}\n')
                continue

        if mol is None:
            out_file.write(f'# {id} No_output\n')
            continue
    
        counts = ' '.join(str(atom.GetTotalNumHs()) for atom in mol.GetAtoms())

        out_file.write(f'{id} {counts}\n')