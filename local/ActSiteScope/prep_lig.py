import sys
import argparse
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTMolecule

def prep_ligand(sdf_in, pdbqt_out):
    mol = Chem.MolFromMolFile(sdf_in, removeHs=False, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_FINDRADICALS | Chem.SanitizeFlags.SANITIZE_KEKULIZE | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION | Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
    mol = Chem.AddHs(mol, addCoords=True)
    
    prep = MoleculePreparation()
    prep.prepare(mol)
    prep.write_pdbqt_file(pdbqt_out)

if __name__ == '__main__':
    prep_ligand(sys.argv[1], sys.argv[2])
