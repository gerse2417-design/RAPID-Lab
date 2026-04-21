import os, sys, glob
sys.path.append('.')
from scripts.validate_binding import VinaValidator

curr_dir = os.path.dirname(os.path.abspath(__file__))
val = VinaValidator(curr_dir)
pdb = "inputs/clean_1m40.pdb1"
lig = "inputs/1m40_I_CB4.sdf"
center = [15.0, -2.0, 41.0]
size = [20.0, 20.0, 20.0]

print("Preparing Receptor...")
cmd = f"{val.venv_python} scripts/prepare_vina.py receptor {pdb} results/site_TEST_receptor.pdbqt"
os.system(cmd)

print("Preparing Ligand...")
cmd = f"{val.venv_python} scripts/prepare_vina.py ligand {lig} results/site_TEST_ligand.pdbqt"
os.system(cmd)

print("Running Vina...")
cmd = f"scripts/bin/vina --receptor results/site_TEST_receptor.pdbqt --ligand results/site_TEST_ligand.pdbqt --center_x {center[0]} --center_y {center[1]} --center_z {center[2]} --size_x {size[0]} --size_y {size[1]} --size_z {size[2]} --out results/site_TEST_out.pdbqt > results/site_TEST_docking.log 2>&1"
os.system(cmd)
with open("results/site_TEST_docking.log") as f:
    print("Vina log:")
    print(f.read())
