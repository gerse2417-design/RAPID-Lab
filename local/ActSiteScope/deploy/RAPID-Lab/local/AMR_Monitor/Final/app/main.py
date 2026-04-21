import os
import argparse
import sys
from dotenv import load_dotenv

# Try to load environments specifically for AWS CLI context
env_paths = [os.path.join(os.path.dirname(__file__), "..", "config", ".env.aws"), ".env"]
for env_path in env_paths:
    if os.path.exists(env_path):
        load_dotenv(env_path)
        print(f"Loaded environment variables from {env_path}")
        break

# Add parent directory to path to import pipeline
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline import FinalWorkflowPipeline
from database import init_db

def main():
    parser = argparse.ArgumentParser(description="Bioinformatics Command Line Pipeline")
    
    INPUT_DIR = os.environ.get("INPUT_DIR")
    OUTPUT_DIR = os.environ.get("OUTPUT_DIR", "results")

    # Input/Output Arguments
    parser.add_argument("--input", "-i", default=INPUT_DIR, help="Input FASTA file path or S3 mapped directory")
    parser.add_argument("--output", "-o", default=OUTPUT_DIR, help="Output directory")
    
    # RADAR arguments
    parser.add_argument("--cutoff", "-c", type=float, default=float(os.environ.get("RADAR_CUTOFF", 0.95)), help="Identity cutoff for RADAR BOARDS")
    parser.add_argument("--sample-name", "-s", default=os.environ.get("RADAR_SAMPLE_NAME", ""), help="Sample name")
    parser.add_argument("--strain-tag", "-t", default=os.environ.get("RADAR_STRAIN_TAG", ""), help="Strain tag")
    
    # AlphaFold arguments
    parser.add_argument("--run-alphafold", action="store_true", help="Execute AlphaFold2 after RADAR engine")
    parser.add_argument("--af-mode", default=os.environ.get("AF2_MODE", "local"), choices=["local", "docker"], help="AlphaFold2 runtime mode")
    parser.add_argument("--msa-mode", default=os.environ.get("AF2_MSA_MODE", "mmseqs2_uniref_env"), help="MSA Mode")
    parser.add_argument("--num-models", type=int, default=int(os.environ.get("AF2_NUM_MODELS", 5)), help="Number of models to generate")
    parser.add_argument("--num-recycles", type=int, default=int(os.environ.get("AF2_NUM_RECYCLES", 3)), help="Number of recycles")
    
    args = parser.parse_args()

    if not args.input:
        print("[Error] No input specified. Provide --input or set INPUT_DIR environment variable.")
        sys.exit(1)
        
    input_fasta = os.path.abspath(args.input)
    if not os.path.exists(input_fasta):
        print(f"[Error] Input file not found: {input_fasta}")
        sys.exit(1)
        
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)

    print(f"========== Bioinformatics Pipeline CLI ==========")
    print(f"Input: {input_fasta}")
    print(f"Output Directory: {output_dir}")

    session = init_db(db_path="sqlite:////tmp/final_workflow.db")
    custom_config = {
        "radar": {"db_name": "BOARDS", "cutoff": args.cutoff, "sample_name": args.sample_name, "strain_tag": args.strain_tag},
        "alphafold": {"mode": args.af_mode, "msa_mode": args.msa_mode, "num_models": args.num_models, "num_recycle": args.num_recycles, "output_dir": output_dir}
    }
    pipeline = FinalWorkflowPipeline(session, custom_config=custom_config)
    
    print("\n--- [Stage 1] Running RADAR Engine ---")
    try:
        is_hit, radar_data = pipeline.run_actual_radar_engine(input_fasta)
        if is_hit:
            print("=> RADAR Analysis Complete: Hit Found.")
            # Copy generated reports from internal radar path to the user's output directory
            import shutil
            import glob
            
            # Extract sample prefix from args
            prefix = args.sample_name.lower().replace(" ", "_") if args.sample_name else "unknown"
            
            # Find TSV reports
            internal_report_dir = os.path.join(os.path.dirname(__file__), "..", "radar", "pipeline", "result", "wgs", "table", prefix)
            if os.path.exists(internal_report_dir):
                for tsv_file in glob.glob(os.path.join(internal_report_dir, "*.tsv")):
                    shutil.copy(tsv_file, output_dir)
                    print(f"   [Copy] Saved report: {os.path.basename(tsv_file)} -> {output_dir}")
            
            # Find target protein
            if radar_data and radar_data.get("query_faa") and os.path.exists(radar_data["query_faa"]):
                shutil.copy(radar_data["query_faa"], output_dir)
                print(f"   [Copy] Saved protein: {os.path.basename(radar_data['query_faa'])} -> {output_dir}")
                
        else:
            print("=> RADAR Analysis Complete: No Hit Found.")
    except Exception as e:
        print(f"[Error] RADAR Engine failed: {e}")
        sys.exit(1)

    if args.run_alphafold:
        print("\n--- [Stage 2] Running AlphaFold2 Engine ---")
        try:
            af2_input_seq = radar_data.get("query_faa") if is_hit and radar_data else None
            if not af2_input_seq or not os.path.exists(af2_input_seq):
                af2_input_seq = input_fasta
                
            af_out_dir = os.path.join(output_dir, "af2_results")
            os.makedirs(af_out_dir, exist_ok=True)
            
            try:
                from af2_module.af2_runner import run_alphafold_local, run_alphafold_docker
            except ImportError:
                print("[Warning] af2_module not found. Calling CLI tool_paths fallback.")
                pipeline.run_alphafold("CLI_RUN", af2_input_seq, af_out_dir)
                run_alphafold_local = None

            if run_alphafold_local:
                if args.af_mode == "docker":
                    docker_img = os.environ.get("AF2_DOCKER_IMAGE", "ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2")
                    gen = run_alphafold_docker(af2_input_seq, af_out_dir, msa_mode=args.msa_mode, docker_image=docker_img, gpu_mem_fraction="0.9", num_models=args.num_models, num_recycles=args.num_recycles)
                else:
                    gen = run_alphafold_local(af2_input_seq, af_out_dir, msa_mode=args.msa_mode, num_models=args.num_models, num_recycles=args.num_recycles, gpu_mem_fraction="0.9")
                
                for log_line in gen:
                    print(log_line, end="")
            print(f"\n=> AlphaFold2 Execution Complete. Results in: {af_out_dir}")
        except Exception as e:
            print(f"[Error] AlphaFold2 Engine failed: {e}")
            sys.exit(1)

    print("\n========== Pipeline Execution Finished ==========")

if __name__ == "__main__":
    main()
