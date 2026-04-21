import os
import sys
import glob
sys.path.append('.')
from scripts.validate_binding import VinaValidator

curr_dir = os.path.dirname(os.path.abspath(__file__))
val = VinaValidator(curr_dir)
pdb = "inputs/clean_1m40.pdb1"
lig = "inputs/1m40_I_CB4.sdf"
center = [15.0, -2.0, 41.0]
size = [20.0, 20.0, 20.0]

print("Starting docking test...")
res = val.run_docking(os.path.abspath(pdb), os.path.abspath(lig), center, size, site_id="TEST")
print("Docking complete.")
print("Result:", res)
