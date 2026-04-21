import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

def is_2d_molecule(mol):
    """
    Check if a molecule is 2D by examining Z-coordinates of all atoms.
    """
    if not mol.GetNumConformers():
        return True
    conf = mol.GetConformer()
    if not conf.Is3D():
        return True
    
    # Check if all Z coordinates are approximately zero
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        if abs(pos.z) > 0.05: # Threshold for 3D detection
            return False
    return True

def prep_ligand(sdf_in, pdbqt_out):
    try:
        # Load molecule with sanitized properties
        mol = Chem.MolFromMolFile(sdf_in, removeHs=False, sanitize=True)
        if not mol:
            # Fallback for sanitization errors
            mol = Chem.MolFromMolFile(sdf_in, removeHs=False, sanitize=False)
            mol.UpdatePropertyCache(strict=False)

        # 1. Automatic 2D -> 3D Conversion & Optimization
        if is_2d_molecule(mol):
            # Add Hydrogens needed for 3D geometry
            mol = Chem.AddHs(mol)
            # Generate 3D Conformation (ETKDGv3 is standard for accuracy)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            success = AllChem.EmbedMolecule(mol, params)
            
            if success != 0:
                # Fallback to random coordinates if ETKDG fails
                AllChem.EmbedMolecule(mol, AllChem.ETKDG(), randomSeed=42)
            
            # 2. Energy Minimization (MMFF94)
            # This ensures realistic bond lengths and angles for Vina
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except:
                AllChem.UFFOptimizeMolecule(mol)
        else:
            # For existing 3D structures, still add Hs and perform a quick cleanup
            mol = Chem.AddHs(mol, addCoords=True)
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=100)
            except:
                pass

        # 3. PDBQT Conversion using Meeko
        prep = MoleculePreparation()
        prep.prepare(mol)
        prep.write_pdbqt_file(pdbqt_out)

        # Post-process ligand PDBQT to replace unsupported Vina atom types (e.g. Boron -> Carbon)
        with open(pdbqt_out, "r") as f:
            lines = f.readlines()
        with open(pdbqt_out, "w") as f:
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    etype = line[77:].strip()
                    if etype == "B" or etype == "Si":
                        # Morph unsupported elements to Carbon for empirical scoring
                        line = line[:77] + "C\n"
                f.write(line)
        return True
    except Exception as e:
        print(f"Meeko Ligand Prep Error: {e}", file=sys.stderr)
        return False

def prep_receptor(pdb_in, pdbqt_out):
    try:
        # Simplistic and highly fault-tolerant PDB to PDBQT converter
        # AutoDock Vina's empirical scoring function requires receptor coordinates and AD4 types.
        # It ignores receptor partial charges completely.
        with open(pdb_in, "r") as f:
            lines = f.readlines()
        
        with open(pdbqt_out, "w") as f:
            for line in lines:
                if line.startswith("ATOM  "):
                    element = line[76:78].strip()
                    if not element:
                        name = line[12:16].strip()
                        element = ''.join([c for c in name if c.isalpha()])[0]
                    
                    ad4_type = element
                    if element == "C": ad4_type = "C"
                    elif element == "N": ad4_type = "N"
                    elif element == "O": ad4_type = "OA"
                    elif element == "S": ad4_type = "SA"
                    elif element == "H": ad4_type = "HD"
                    elif element == "P": ad4_type = "P"
                    else: ad4_type = "C"

                    part1 = line[:66].ljust(66)
                    charge = " 0.000"
                    part2 = ad4_type.ljust(2)
                    new_line = f"{part1}    {charge} {part2}\n"
                    f.write(new_line)
                elif line.startswith("TER"):
                    f.write("TER\n")
        return True
    except Exception as e:
        print(f"Receptor Prep Error: {e}", file=sys.stderr)
        return False

if __name__ == '__main__':
    mode = sys.argv[1]
    fin = sys.argv[2]
    fout = sys.argv[3]
    if mode == "ligand":
        success = prep_ligand(fin, fout)
    elif mode == "receptor":
        success = prep_receptor(fin, fout)
    else:
        success = False
    sys.exit(0 if success else 1)
