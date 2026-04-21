import os
import sys
import argparse
import logging
import shutil
import glob
import numpy as np
from pathlib import Path
from dotenv import load_dotenv

# Add project root to sys.path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)

from scripts.mcsa_utils import get_mcsa_residues_by_uniprot, get_mcsa_residues_by_pdb_id, get_mcsa_residues_by_keyword
from scripts.pdb_utils import extract_pdb_info, clean_pdb_file
from scripts.pocket_utils import run_p2rank, run_fpocket
from scripts.apbs_utils import run_apbs_pipeline, score_pocket_electrostatics
from scripts.predict_site import get_real_mpbind_predictions, get_residue_center
from scripts.consensus_analysis import ConsensusEngine
from scripts.validate_binding import VinaValidator

# Setup Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler()]
)

def download_s3_path(s3_uri, local_dir):
    """
    Detect S3 URI and download to local directory.
    Requires boto3 and active AWS credentials.
    """
    if not s3_uri.startswith("s3://"):
        return s3_uri

    try:
        import boto3
        from urllib.parse import urlparse
        
        logging.info(f"[*] S3 Download Started: {s3_uri}")
        parsed = urlparse(s3_uri)
        bucket_name = parsed.netloc
        s3_key = parsed.path.lstrip('/')
        
        s3 = boto3.client('s3')
        local_path = os.path.join(local_dir, os.path.basename(s3_key))
        os.makedirs(local_dir, exist_ok=True)
        
        s3.download_file(bucket_name, s3_key, local_path)
        return local_path
    except ImportError:
        logging.error("boto3 is not installed. S3 download failed.")
        return s3_uri
    except Exception as e:
        logging.error(f"S3 Download Error: {e}")
        return s3_uri

def run_pipeline(args):
    """
    Core Analytical Pipeline. 
    Can be called from CLI (main) or UI (streamlit_app).
    """
    # 1. Environment Settings
    env_path = getattr(args, 'env', None) or os.path.join(BASE_DIR, "config", ".env.aws")
    if os.path.exists(env_path):
        load_dotenv(env_path)
        logging.info(f"[*] Environment Config Loaded: {env_path}")
    
    # Defaults from Env
    INPUT_DIR = os.environ.get("INPUT_DIR", os.path.join(BASE_DIR, "inputs"))
    RESULT_DIR = getattr(args, 'output', None) or os.environ.get("RESULT_DIR", os.path.join(BASE_DIR, "results"))
    
    os.makedirs(INPUT_DIR, exist_ok=True)
    os.makedirs(RESULT_DIR, exist_ok=True)

    # 2. File Preparation (Handles S3 or Local)
    pdb_path = download_s3_path(args.pdb, INPUT_DIR)
    ligand_paths = []
    if getattr(args, 'ligands', None):
        for lp in args.ligands:
            ligand_paths.append(download_s3_path(lp, INPUT_DIR))

    logging.info(f"[*] Starting Analysis for: {pdb_path}")
    
    # 3. PDB Preprocessing
    clean_path = os.path.join(INPUT_DIR, f"clean_{os.path.basename(pdb_path)}")
    if not clean_pdb_file(pdb_path, clean_path):
        logging.warning("PDB cleaning failed, using original file.")
        clean_path = pdb_path
    
    with open(clean_path, "r") as f:
        pdb_data = f.read()
    pdb_info = extract_pdb_info(pdb_data)
    
    # 4. M-CSA Standard Catalytic Motif Search
    logging.info("Step 1: M-CSA Catalytic Motif Search...")
    mcsa_residues = []
    if getattr(args, 'uniprot', None):
        mcsa_residues = get_mcsa_residues_by_uniprot(args.uniprot)
    if not mcsa_residues and pdb_info.get("pdb_id"):
        mcsa_residues = get_mcsa_residues_by_pdb_id(pdb_info.get("pdb_id"))
    if not mcsa_residues and pdb_info.get("title_keyword"):
        mcsa_residues = get_mcsa_residues_by_keyword(pdb_info.get("title_keyword"))
    
    enriched_mcsa = []
    for res in mcsa_residues:
        c = get_residue_center(pdb_data, str(res.get("res_num", "")))
        if c: res["center"] = c
        enriched_mcsa.append(res)
    logging.info(f"   -> {len(enriched_mcsa)} Catalytic residues found.")

    # 5. AI-Based Binding Prediction (MPBind)
    logging.info("Step 2: AI-Based Binding Site Prediction (MPBind)...")
    mp_hotspots = [h for h in get_real_mpbind_predictions(clean_path) if h.get("prob", 0) >= args.th_mp]
    logging.info(f"   -> {len(mp_hotspots)} AI hotspots detected (threshold {args.th_mp}).")

    # 6. Pocket Search (P2Rank & fPocket)
    logging.info("Step 3: Geometry-Based Pocket Search (P2Rank, fPocket)...")
    p2r_all = run_p2rank(clean_path, RESULT_DIR)
    p2rank_filtered = []
    if p2r_all:
        mx = max(p.get("score", 0) for p in p2r_all)
        p2r_th = 1.0 - (args.th_p2r / 100.0)
        p2rank_filtered = [p for p in p2r_all if p.get("score", 0) >= mx * p2r_th]
    
    fpocket_results = [f for f in run_fpocket(clean_path, RESULT_DIR) if f.get("score", 0) >= args.th_fpc]
    logging.info(f"   -> P2Rank: {len(p2rank_filtered)}, fPocket: {len(fpocket_results)} pockets detected.")

    # 7. APBS Electrostatic Potential Generation
    logging.info("Step 4: APBS Electrostatic Potential Calculation...")
    dx_path = run_apbs_pipeline(clean_path, RESULT_DIR)
    if dx_path:
        for key_list in [mp_hotspots, p2rank_filtered, fpocket_results]:
            for p in key_list:
                ns, rs = score_pocket_electrostatics(p.get("residues", ""), dx_path, clean_path)
                p["apbs_score"], p["raw_apbs_score"] = ns, rs
        logging.info("   -> Electrostatic mapping completed.")

    # 8. Consensus Model Integration
    logging.info("Step 5: Consensus Model Integration (Consensus Hotspots)...")
    engine = ConsensusEngine(
        clean_path, mp_hotspots, p2rank_filtered, fpocket_results, enriched_mcsa,
        mp_th=args.th_mp, p2r_pct=args.th_p2r, fpc_th=args.th_fpc, apbs_pct=args.th_apbs
    )
    final_res = engine.run()
    
    enriched_final = []
    for res in final_res:
        c = get_residue_center(pdb_data, str(res.get("ResNo", "")))
        if c: res["center"] = c
        enriched_final.append(res)
    
    report_path = os.path.join(RESULT_DIR, "Final_Hotspot_Report.txt")
    engine.generate_report(
        enriched_final, report_path, 
        all_mcsa=enriched_mcsa, all_mp=mp_hotspots, 
        all_p2r=p2rank_filtered, all_fpc=fpocket_results
    )
    logging.info(f"   -> Consensus analysis completed. Report: {report_path}")

    # 9. AutoDock Vina Docking Validation
    vina_summaries = {}
    if ligand_paths and enriched_final:
        logging.info("Step 6: AutoDock Vina Docking Validation...")
        validator = VinaValidator(BASE_DIR)
        target_coords = [h["center"] for h in enriched_final if "center" in h]
        if target_coords:
            center_point = np.mean(target_coords, axis=0).tolist()
            for l_path in ligand_paths:
                l_name = os.path.basename(l_path)
                logging.info(f"   - Docking: {l_name} ...")
                res = validator.run_docking(clean_path, l_path, center_point, [20, 20, 20], output_dir=RESULT_DIR, site_id=0)
                if "affinity" in res:
                    vina_summaries[l_name] = res["affinity"]
                    logging.info(f"     -> Result: {res['affinity']:.3f} kcal/mol")
        
        # Final Report Update with Vina
        engine.generate_report(
            enriched_final, report_path, vina_results=vina_summaries,
            all_mcsa=enriched_mcsa, all_mp=mp_hotspots, 
            all_p2r=p2rank_filtered, all_fpc=fpocket_results
        )

    logging.info(f"[*] Pipeline Complete! Results saved in: {RESULT_DIR}")
    
    # Return results for UI consumption
    return {
        "enriched_final": enriched_final,
        "enriched_mcsa": enriched_mcsa,
        "mp_hotspots": mp_hotspots,
        "p2rank_results": p2rank_filtered,
        "fpocket_results": fpocket_results,
        "vina_results": vina_summaries,
        "report_path": report_path
    }

def main():
    parser = argparse.ArgumentParser(description="ActSiteScope Hotspot Prediction Pipeline (CLI Version)")
    parser.add_argument("--pdb", required=True, help="Input Protein PDB path (Local or S3)")
    parser.add_argument("--uniprot", help="UniProt ID for M-CSA lookup")
    parser.add_argument("--ligands", nargs="*", help="List of Ligand files (SDF/MOL2/PDB) for batch docking")
    parser.add_argument("--env", help="Path to .env file (Default: config/.env.aws)")
    parser.add_argument("--output", help="Override RESULT_DIR from env")
    
    # Threshold Overrides
    parser.add_argument("--th_mp", type=float, default=0.5, help="MPBind probability threshold")
    parser.add_argument("--th_p2r", type=int, default=30, help="P2Rank top percentage threshold")
    parser.add_argument("--th_fpc", type=float, default=0.5, help="fPocket score threshold")
    parser.add_argument("--th_apbs", type=int, default=30, help="APBS top percentage threshold")
    
    args = parser.parse_args()
    
    run_pipeline(args)

if __name__ == "__main__":
    main()
