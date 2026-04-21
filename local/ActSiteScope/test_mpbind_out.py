import sys
import json
from scripts.predict_site import get_real_mpbind_predictions

res = get_real_mpbind_predictions("results/tmp_input/1m40.pdb1", ".")
for hotspot in res:
    print(hotspot["site"], hotspot["residues"])

