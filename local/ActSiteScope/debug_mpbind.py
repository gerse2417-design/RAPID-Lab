import os
import sys

# load streamlit_app as a module and test format_residue_info_columns
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app.streamlit_app import get_pdb_residue_map, format_residue_info_columns

with open("results/tmp_input/1m40.pdb1", "r") as f:
    pdb_data = f.read()

rm = get_pdb_residue_map(pdb_data)

# simulate an MPBind output
test_str = "UNK 26, UNK 71, UNK 98"
ident, names = format_residue_info_columns(test_str, rm)
print("IDENT:", ident)
print("NAMES:", names)

