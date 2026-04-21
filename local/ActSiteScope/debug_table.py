import os
import sys

from app.streamlit_app import get_pdb_residue_map, format_residue_info_columns, generate_styled_table
from scripts.predict_site import get_real_mpbind_predictions

pdb_path = "results/tmp_input/1m40.pdb1"
with open(pdb_path, "r") as f:
    pdb_data = f.read()

rm = get_pdb_residue_map(pdb_data)

# Real MPBind Output
res = get_real_mpbind_predictions(pdb_path, "deploy/actsitescope-full")

print("Real Hotspots:")
for h in res:
    print(h)

print("\nTable:")
html = generate_styled_table(res, ["site","prob","residues"], rm, widths=[15,20,35,30])
print(html)
