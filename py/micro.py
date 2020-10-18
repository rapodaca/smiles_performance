import sys
import time
from rdkit import Chem

in_path = "./data/rdkit_2017.03.3.smi"

with open(in_path) as in_file:
    inputs = [ ]
    
    for line in in_file:
        parts = line.strip().split()
        
        if len(parts) == 1:
            continue

        inputs.append(parts[0])

    start = time.time()

    for smiles in inputs:
        Chem.MolFromSmiles(smiles)

    print("elapsed: " + str(time.time() - start))
