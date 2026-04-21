import sys
import os
from meeko import MoleculePreparation, Polymer, ResidueChemTemplates
from rdkit import Chem

def fix_pdb_and_prep(pdb_in, out_pdbqt):
    with open(pdb_in, "r") as f:
        lines = f.readlines()
        
    filtered = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            filtered.append(line)
        elif line.startswith("TER") or line.startswith("END"):
            filtered.append(line)
            
    pdb_str = "".join(filtered)
    try:
        # PDBQTWriterLegacy uses Meeko's internal PDB writer
        mk_prep = MoleculePreparation.from_config({"charge_model": "gasteiger"})
        templates = ResidueChemTemplates.create_from_defaults()
        polymer = Polymer.from_pdb_string(pdb_str, templates, mk_prep, {}, [], True)
        from meeko.writer import PDBQTWriterLegacy
        rigid_pdbqt, flex_pdbqt_dict = PDBQTWriterLegacy.write_from_polymer(polymer)
        with open(out_pdbqt, "w") as f:
            f.write(rigid_pdbqt)
        print("Success!")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    fix_pdb_and_prep(sys.argv[1], sys.argv[2])
