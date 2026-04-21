import sys, os
sys.path.append('/home/kcpak/projects/hotspot_prediction')
from scripts.pocket_utils import run_p2rank
print("Testing p2rank on 1m40.pdb1")
res, err = run_p2rank('/home/kcpak/projects/hotspot_prediction/inputs/1m40.pdb1', '/home/kcpak/projects/hotspot_prediction/results')
print(f"Results len: {len(res)}")
print(f"Err: {err}")
