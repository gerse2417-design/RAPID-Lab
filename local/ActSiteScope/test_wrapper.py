import os
import sys
# Add scripts to path
sys.path.append(os.path.abspath("scripts"))
from predict_site import MPBindWrapper

def test_mpbind():
    print("Testing MPBindWrapper...")
    wrapper = MPBindWrapper(project_root=os.getcwd())
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(curr_dir, "inputs", "1m40.pdb1")
    output_dir = os.path.join(curr_dir, "results", "test_wrapper")
    os.makedirs(output_dir, exist_ok=True)
    
    results = wrapper.predict(pdb_path, output_dir)
    if results:
        print(f"Success! Generated files: {results}")
    else:
        print("Failed to generate results.")
        if os.path.exists("mpbind_error.log"):
            with open("mpbind_error.log", "r") as f:
                print("Error Log Content:")
                print(f.read())

if __name__ == "__main__":
    test_mpbind()
