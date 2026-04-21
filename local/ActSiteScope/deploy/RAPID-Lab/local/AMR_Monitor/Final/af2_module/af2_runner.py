import os
import subprocess
import argparse
import sys
import time
from datetime import datetime

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def _run_generator(cmd, env=None):
    """Helper to run a command and yield stdout lines in real-time."""
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, env=env
    )
    for line in process.stdout:
        yield line
    process.stdout.close()
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd)

def run_alphafold_docker(input_file, output_dir, docker_image="ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2", 
                         msa_mode="mmseqs2_uniref_env", gpu_mem_fraction="0.7",
                         num_models=5, num_recycles=3):
    """
    Runs LocalColabFold using the official Docker image.
    Yields stdout lines for real-time progress tracking.
    """
    input_abs = os.path.abspath(input_file)
    output_abs = os.path.abspath(output_dir)
    os.makedirs(output_abs, exist_ok=True)
    
    cache_abs = os.path.join(SCRIPT_DIR, "af2_cache")
    jax_cache_abs = os.path.join(SCRIPT_DIR, "af2_jax_cache")
    os.makedirs(cache_abs, exist_ok=True)
    os.makedirs(jax_cache_abs, exist_ok=True)
    
    cmd = [
        "docker", "run", "--rm", "--gpus", "all",
        "-e", "XLA_PYTHON_CLIENT_PREALLOCATE=false",
        "-e", f"XLA_PYTHON_CLIENT_MEM_FRACTION={gpu_mem_fraction}",
        "-e", "TF_FORCE_UNIFIED_MEMORY=1",
        "-e", "JAX_COMPILATION_CACHE_DIR=/jax_cache",
        "-e", "JAX_COMPILATION_CACHE_MIN_COMPILE_TIME_SECS=0",
        "-e", "COLABFOLD_CACHE_DIR=/cache",
        "-v", f"{input_abs}:/input.fasta:ro",
        "-v", f"{output_abs}:/output",
        "-v", f"{cache_abs}:/cache",
        "-v", f"{jax_cache_abs}:/jax_cache",
        docker_image,
        "colabfold_batch",
        "/input.fasta",
        "/output",
        "--msa-mode", msa_mode,
        "--num-models", str(num_models),
        "--num-recycle", str(num_recycles)
    ]
    
    yield f"🚀 Starting AlphaFold2 Prediction via Docker ({docker_image})...\n"
    yield f"Input: {input_abs}\nOutput: {output_abs}\n"
    yield "-"*50 + "\n"
    
    try:
        yield from _run_generator(cmd)
    except Exception as e:
        yield f"\n❌ Error during AlphaFold2 execution: {e}\n"
        raise

def run_alphafold_local(input_file, output_dir, msa_mode="mmseqs2_uniref_env", 
                        num_models=5, num_recycles=3, gpu_mem_fraction="0.90"):
    """
    Runs LocalColabFold using the 'af_pipeline' Conda environment.
    Yields stdout lines for real-time progress tracking.
    
    환경 변수 CONDA_ENV_LIB_PATH 가 설정되어 있으면 해당 경로를 사용합니다.
    (로컬: 미설정 시 하드코딩 경로 사용 / Docker: Dockerfile에서 ENV로 주입)
    """
    input_abs = os.path.abspath(input_file)
    output_abs = os.path.abspath(output_dir)
    os.makedirs(output_abs, exist_ok=True)
    
    import shutil
    has_conda = shutil.which("conda") is not None

    env = os.environ.copy()
    
    # [Peak Performance & Stability] JAX/CUDA 환경 최적화
    # 2. XLA/GPU 성능 및 안정성 플래그
    env["HOME"] = "/tmp"  # non-root user를 위한 홈 디렉토리 우회
    env["MPLCONFIGDIR"] = "/tmp/matplotlib"
    env["XDG_CACHE_HOME"] = "/tmp/.cache"
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false" 
    env["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
    env["TF_FORCE_GPU_ALLOW_GROWTH"] = "true" 
    env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.85"
    env["TF_FORCE_UNIFIED_MEMORY"] = "1"
    env["JAX_COMPILATION_CACHE_DIR"] = "/tmp/jax_cache"
    env["JAX_COMPILATION_CACHE_MIN_COMPILE_TIME_SECS"] = "0"
    
    os.makedirs("/tmp/jax_cache", exist_ok=True)
    os.makedirs("/tmp/matplotlib", exist_ok=True)
    os.makedirs("/tmp/.cache", exist_ok=True)
    
    if has_conda:
        # 1. 라이브러리 경로: 환경 변수 우선 → 없으면 로컬 하드코딩 경로 폴백
        _LOCAL_FALLBACK = "/home/chlee/miniconda3/envs/af_pipeline_v2/lib"
        conda_lib = os.environ.get("CONDA_ENV_LIB_PATH", _LOCAL_FALLBACK)
        env["LD_LIBRARY_PATH"] = f"{conda_lib}:{env.get('LD_LIBRARY_PATH', '')}"
        env["CUDA_DIR"] = conda_lib
        env["XLA_FLAGS"] = f"--xla_gpu_cuda_data_dir={conda_lib} --xla_gpu_graph_level=0 --xla_gpu_enable_triton_gemm=false"
    else:
        conda_lib = "/usr/local/cuda"
        env["XLA_FLAGS"] = "--xla_gpu_graph_level=0 --xla_gpu_enable_triton_gemm=false"
    
    conda_env_name = os.environ.get("CONDA_ENV_NAME", "af_pipeline_v2")
    
    import shutil
    has_conda = shutil.which("conda") is not None

    if has_conda:
        cmd = [
            "conda", "run", "-n", conda_env_name, "--no-capture-output",
            "colabfold_batch",
            input_abs,
            output_abs,
            "--msa-mode", msa_mode,
            "--num-models", str(num_models),
            "--num-recycle", str(num_recycles)
        ]
        runner_msg = f"🚀 Starting AlphaFold2 Prediction via Local (Conda: {conda_env_name})...\n"
    else:
        cmd = [
            "colabfold_batch",
            input_abs,
            output_abs,
            "--msa-mode", msa_mode,
            "--num-models", str(num_models),
            "--num-recycle", str(num_recycles)
        ]
        runner_msg = f"🚀 Starting AlphaFold2 Prediction via Local (Native colabfold_batch)...\n"

    yield runner_msg
    yield f"   conda_lib: {conda_lib}\n"
    yield f"Input: {input_abs}\nOutput: {output_abs}\n"
    yield "-"*50 + "\n"
    
    try:
        yield from _run_generator(cmd, env=env)
    except Exception as e:
        yield f"\n❌ Error during AlphaFold2 execution: {e}\n"
        raise

if __name__ == "__main__":
    # Standard CLI fallback for run_pipeline.sh compatibility
    parser = argparse.ArgumentParser(description="AlphaFold2 Local/Docker Runner")
    parser.add_argument("input", help="Path to input .faa (FASTA) file")
    parser.add_argument("--output", default="af2_results", help="Output directory name")
    parser.add_argument("--mode", choices=["docker", "local"], default="docker")
    
    args = parser.parse_args()
    
    if args.mode == "docker":
        for line in run_alphafold_docker(args.input, args.output):
            print(line, end="")
    else:
        for line in run_alphafold_local(args.input, args.output):
            print(line, end="")
