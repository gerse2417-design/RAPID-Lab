import os
import subprocess
import logging
import pandas as pd
import numpy as np

def run_fpocket(pdb_path, output_dir):
    """
    Execute fpocket and parse results.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # 환경 변수 지원 (AWS 배포 환경 대응)
    fpocket_bin = os.environ.get('FPOCKET_PATH') or os.path.join(script_dir, "bin", "fpocket", "fpocket")
    
    if not (os.path.exists(pdb_path) and os.path.exists(fpocket_bin)):
        logging.error("PDB or Fpocket binary not found.")
        return []
    
    try:
        # fpocket creates a directory with _out suffix in the SAME folder as PDB
        # To avoid cluttering inputs, we copy PDB to output_dir first (optional but safer)
        cmd = [fpocket_bin, "-f", pdb_path]
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Output is in [pdb_basename]_out folder
        pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
        pdb_dir = os.path.dirname(pdb_path)
        out_folder = os.path.join(pdb_dir, f"{pdb_base}_out")
        info_file = os.path.join(out_folder, f"{pdb_base}_info.txt")
        
        if not os.path.exists(info_file):
            logging.warning(f"Fpocket info file not found: {info_file}")
            return []
            
        results = []
        with open(info_file, "r") as f:
            lines = f.readlines()
            
        current_pocket = None
        for line in lines:
            if line.startswith("Pocket"):
                parts = line.split()
                if len(parts) >= 2:
                    site_id = parts[1].replace(":", "")
                    current_pocket = {
                        "id": int(site_id),
                        "site": f"Site {site_id}",
                        "score": 0.0,
                        "prob": 0.0,
                        "center": None,
                        "residues": []
                    }
                    results.append(current_pocket)
            elif "Score" in line and current_pocket:
                try: 
                    current_pocket["score"] = float(line.split(":")[1].strip())
                    current_pocket["prob"] = current_pocket["score"]
                except: pass
            elif "Center" in line and current_pocket:
                try:
                    vals = line.split(":")[1].strip().split()
                    current_pocket["center"] = [float(vals[0]), float(vals[1]), float(vals[2])]
                except: pass
        
        # For residues, we need to parse the pocketX_atm.pdb files in [out_folder]/pockets
        for p in results:
            p_file = os.path.join(out_folder, "pockets", f"pocket{p['id']}_atm.pdb")
            if os.path.exists(p_file):
                res_set = set()
                with open(p_file, "r") as pf:
                    for pline in pf:
                        if pline.startswith("ATOM"):
                            res_name = pline[17:20].strip()
                            res_seq = pline[22:26].strip()
                            res_set.add(f"{res_name} {res_seq}")
                p["residues"] = ", ".join(sorted(list(res_set)))
                
        return results
    except Exception as e:
        logging.error(f"Error running fpocket: {e}")
        return []

def run_p2rank(pdb_path, output_dir):
    """
    Execute p2rank using a local JRE and parse results.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_bin = os.path.join(script_dir, "bin")
    
    # 환경 변수 지원 (AWS 배포 환경 대응)
    p2rank_bin = os.environ.get('P2RANK_PATH') or os.path.join(base_bin, "p2rank", "prank")
    # Docker 등 외부 JRE가 설치된 환경 지원
    jre_exec = os.environ.get('JRE_PATH') or os.path.join(base_bin, "jre", "bin", "java")
    jre_bin = os.path.dirname(jre_exec)
    
    import shutil
    if not os.path.exists(p2rank_bin):
        logging.error(f"P2Rank binary not found: {p2rank_bin}")
        return [], {"error": f"P2Rank 실행 파일을 찾을 수 없습니다: {p2rank_bin}"}
    
    if not (os.path.exists(jre_exec) or shutil.which(jre_exec)):
        logging.error(f"JRE not found: {jre_exec}")
        return [], {"error": f"Java 실행 환경(JRE)을 찾을 수 없습니다: {jre_exec}"}
    
    try:
        # Create a unique output subfolder
        pdb_basename = os.path.basename(pdb_path)
        p2rank_out_dir = os.path.join(output_dir, "p2rank_out")
        os.makedirs(p2rank_out_dir, exist_ok=True)
        
        p2rank_out = os.path.join(output_dir, "p2rank_out")
        os.makedirs(p2rank_out, exist_ok=True)
        
        import uuid
        unique_id = str(uuid.uuid4())[:8]
        tmp_pdb = os.path.join(output_dir, f"temp_input_{unique_id}.pdb")
        if os.path.lexists(tmp_pdb): os.remove(tmp_pdb)
        try:
            import shutil
            shutil.copy(pdb_path, tmp_pdb)
        except Exception as e:
            logging.error(f"Failed to copy PDB for P2Rank: {e}")
            return [], {"error": f"파일 복사 실패: {e}"}
        
        # Use 'prank' script
        cmd = [p2rank_bin, "predict", "-f", tmp_pdb, "-o", p2rank_out]
        
        env = os.environ.copy()
        env["PATH"] = f"{jre_bin}:" + env.get("PATH", "")
        env["JAVA_HOME"] = os.path.dirname(jre_bin)
        
        try:
            result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                logging.error(f"P2Rank execution failed with code {result.returncode}. Stderr: {result.stderr}")
                if os.path.lexists(tmp_pdb): os.remove(tmp_pdb)
                return [], {"error": f"P2Rank 실행 실패 (Code {result.returncode}): {result.stderr[:200]}"}
        except Exception as e:
            logging.error(f"P2Rank subprocess exception: {str(e)}")
            if os.path.lexists(tmp_pdb): os.remove(tmp_pdb)
            return [], {"error": f"P2Rank 실행 중 예외 발생: {str(e)}"}

        # Clean up temp file
        if os.path.lexists(tmp_pdb): os.remove(tmp_pdb)

        # In newer P2Rank, the output CSV is named after the input file
        expected_csv = os.path.join(p2rank_out, f"temp_input_{unique_id}.pdb_predictions.csv")
        if not os.path.exists(expected_csv):
            # Try alternate naming
            alt_csv = os.path.join(p2rank_out, "temp_input_predictions.csv")
            if os.path.exists(alt_csv):
                expected_csv = alt_csv
            else:
                logging.warning(f"P2Rank predictions file not found in {p2rank_out}")
                return [], {"error": "P2Rank 분석 결과 파일(CSV)이 생성되지 않았습니다."}
                
        # 1. Pockets CSV
        df = pd.read_csv(expected_csv)
        df.columns = [c.strip() for c in df.columns]
        
        results = []
        for i, row in df.iterrows():
            results.append({
                "id": i + 1,
                "site": f"Site {i+1}",
                "score": float(row['score']),
                "prob": float(row.get('probability', row['score'])),
                "center": [float(row['center_x']), float(row['center_y']), float(row['center_z'])],
                "residues": row['residue_ids']
            })

        # 2. Residues CSV (for detailed tie-break)
        res_scores = {}
        res_csv = expected_csv.replace("_predictions.csv", "_residues.csv")
        if not os.path.exists(res_csv):
            res_csv = expected_csv.replace(".pdb_predictions.csv", ".pdb_residues.csv")
            
        if os.path.exists(res_csv):
            try:
                res_df = pd.read_csv(res_csv)
                res_df.columns = [c.strip() for c in res_df.columns]
                for _, r in res_df.iterrows():
                    key = (r['chain'], int(r['residue']))
                    res_scores[key] = float(r['score'])
            except: pass
        
        # Fallback: if no residues.csv, map pocket score to its residues
        if not res_scores:
            for p in results:
                for rid in p['residues'].split():
                    try:
                        # P2Rank residues.ids format: Chain.ResNo
                        parts = rid.split('.')
                        if len(parts) == 2:
                            res_scores[(parts[0], int(parts[1]))] = p['score']
                    except: continue

        return results, res_scores
    except Exception as e:
        logging.error(f"Error running P2Rank: {e}")
        return [], {}

def simulate_fpocket_results(pdb_path):
    """Return mock pocket data for development/demo."""
    # Simulated pockets based on some realistic defaults
    return [
        {"id": 1, "site": "Site 1", "score": 0.542, "prob": 0.542, "center": [10.0, 20.0, 30.0], "residues": "SER 70, LYS 73, SER 130"},
        {"id": 2, "site": "Site 2", "score": 0.421, "prob": 0.421, "center": [15.0, 25.0, 35.0], "residues": "ASP 101, HIS 105"}
    ]

def simulate_p2rank_results(pdb_path):
    """Return mock p2rank pocket data."""
    return [
        {"id": 1, "site": "Site 1", "score": 0.85, "prob": 0.85, "center": [10.2, 20.1, 30.3], "residues": "SER 70, LYS 73, SER 130, GLU 166"},
        {"id": 2, "site": "Site 2", "score": 0.72, "prob": 0.72, "center": [15.5, 24.8, 35.2], "residues": "ASP 101, GLY 102"}
    ]

def get_pocket_overlap(pocket1_res, pocket2_res):
    """Calculate overlap between two residue lists."""
    res1 = set([r.strip() for r in pocket1_res.split(",")])
    res2 = set([r.strip() for r in pocket2_res.split(",")])
    intersection = res1.intersection(res2)
    return len(intersection) / max(len(res1), len(res2))

if __name__ == "__main__":
    print("Testing Pocket Simulation...")
    pockets = simulate_p2rank_results("dummy.pdb")
    print(pockets)
