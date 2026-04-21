import sys
import os

curr_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(curr_dir)

from scripts.apbs_utils import run_apbs_pipeline
test_pdb = os.path.join(curr_dir, "inputs", "orig_1m40.pdb1")
out_dir = os.path.join(curr_dir, "results")
res = run_apbs_pipeline(test_pdb, out_dir)
print("APBS DX Output Path:", res)
