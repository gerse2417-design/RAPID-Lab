import sys
import os

# scripts 폴더를 시스템 경로에 추가
sys.path.append(os.path.join(os.getcwd(), "scripts"))

from pocket_utils import run_p2rank, run_fpocket

def test_pockets():
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(curr_dir, "inputs", "1m40.pdb1")
    output_dir = "results/test_pockets"
    os.makedirs(output_dir, exist_ok=True)
    
    print("Testing P2Rank...")
    p2rank_res = run_p2rank(pdb_path, output_dir)
    print(f"P2Rank results: Found {len(p2rank_res)} pockets")
    if p2rank_res:
        print(f"First pocket residue count: {len(p2rank_res[0].get('residues', []))}")
        
    print("\nTesting Fpocket...")
    fpocket_res = run_fpocket(pdb_path, output_dir)
    print(f"Fpocket results: Found {len(fpocket_res)} pockets")
    if fpocket_res:
        print(f"First pocket residue count: {len(fpocket_res[0].get('residues', []))}")

if __name__ == "__main__":
    test_pockets()
