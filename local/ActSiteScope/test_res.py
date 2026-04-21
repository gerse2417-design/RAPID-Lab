import sys
from app.streamlit_app import get_pdb_residue_map, format_residue_info_columns
from scripts.predict_site import MPBindWrapper

with open("results/tmp_input/1m40.pdb1", "r") as f:
    pdb_str = f.read()

res_map = get_pdb_residue_map(pdb_str)
print("Keys in res_map sample:", list(res_map.keys())[:10])

# MPBind simulation
wrapper = MPBindWrapper("deploy/actsitescope-full")
hotspots = [{"residues": "UNK 26, UNK 71, UNK 98", "key_residues": "UNK 26, UNK 71"}]
for h in hotspots:
    a, b = format_residue_info_columns(h["residues"], res_map)
    print("Formatted:", b)

