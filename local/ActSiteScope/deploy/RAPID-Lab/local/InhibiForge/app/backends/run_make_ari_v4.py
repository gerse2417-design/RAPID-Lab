"""
run_make_ari_v4.py
──────────────────────────────────────────────────────────────
RFdiffusion Pipeline v4 — Docker / AWS CLI 배포용.

v4 변경 사항 (vs v3):
  - 모든 하드코딩 경로를 환경변수로 전환 (python-dotenv 지원)
  - venv_activate를 VENV_CMD 환경변수로 대체 (Docker에선 빈 문자열)
  - GPU 메모리 비율을 GPU_MEMORY_FRACTION 환경변수로 설정
  - argparse 인터페이스 완전 보존

사용 예시 (CLI):
    python run_make_ari_v4.py --name my_job --contigs "100" --num_designs 5 --step all

환경변수 (config/.env.aws 또는 직접 설정):
    RFDIFF_INPUT_DIR, RFDIFF_BASE_DIR, RFDIFF_CODE_DIR, COLABDESIGN_DIR,
    AF2_CACHE_DIR, AF2_PARAM_DIR, RFDIFF_OUTPUT_DIR, VENV_CMD,
    GPU_MEMORY_FRACTION, CUSTOM_LIB_DIR
"""

import os
import sys
import subprocess
import argparse
import time
import datetime
import glob
import math
import re
import json

from dotenv import load_dotenv

# .env 파일 자동 로드 (존재할 경우)
_env_file = os.environ.get("ENV_FILE", "")
if _env_file and os.path.isfile(_env_file):
    load_dotenv(_env_file)
elif os.path.isfile("/config/.env.aws"):
    load_dotenv("/config/.env.aws")

import torch

# ── 환경변수 기반 경로 상수 ────────────────────────────────────────────────
RFDIFF_BASE_DIR = os.environ.get("RFDIFF_BASE_DIR", "/home/sooyeon/amr/RFdiffusion")
RFDIFF_INPUT_DIR = os.environ.get("RFDIFF_INPUT_DIR", os.path.join(RFDIFF_BASE_DIR, "input"))
RFDIFF_CODE_DIR = os.environ.get("RFDIFF_CODE_DIR", os.path.join(RFDIFF_BASE_DIR, "RFdiffusion"))
COLABDESIGN_DIR = os.environ.get("COLABDESIGN_DIR", os.path.join(RFDIFF_BASE_DIR, "ColabDesign_repo"))
AF2_CACHE_DIR = os.environ.get("AF2_CACHE_DIR", os.path.join(RFDIFF_BASE_DIR, "af2_cache"))
AF2_PARAM_DIR = os.environ.get("AF2_PARAM_DIR", os.path.join(AF2_CACHE_DIR, "colabfold"))
RFDIFF_OUTPUT_DIR = os.environ.get("RFDIFF_OUTPUT_DIR", os.path.join(RFDIFF_BASE_DIR, "outputs"))
CUSTOM_LIB_DIR = os.environ.get("CUSTOM_LIB_DIR", "/home/sooyeon/amr/lib")

GPU_MEMORY_FRACTION = float(os.environ.get("GPU_MEMORY_FRACTION", "0.7"))


def _build_venv_activate():
    """
    Docker 환경(VENV_CMD="") → "true" (shell no-op) 반환.
    로컬 환경(VENV_CMD 미설정) → 기존 venv 활성화 + CUDA 경로 설정.
    """
    venv_cmd_env = os.environ.get("VENV_CMD")
    if venv_cmd_env is not None:
        return venv_cmd_env if venv_cmd_env else "true"

    venv_path = os.environ.get("VENV_ACTIVATE", "/home/sooyeon/amr/bin/activate")
    cuda_lib = os.environ.get("CUDA_LIB_DIR", "/usr/local/cuda-12.6/lib64")
    cuda_bin = os.environ.get("CUDA_BIN_DIR", "/usr/local/cuda-12.6/bin")
    return (
        f"source {venv_path} && "
        f"export LD_LIBRARY_PATH={cuda_lib}:{CUSTOM_LIB_DIR}:$LD_LIBRARY_PATH && "
        f"export PATH={cuda_bin}:$PATH && "
        "export E3NN_JIT_COMPILE=0"
    )


def run_command(command, log_file=None, cwd=None, env=None):
    """실시간으로 출력을 확인하며 쉘 커맨드를 실행하는 유틸리티. log_file이 있으면 동시에 파일에도 저장."""
    print(f"Executing: {command}")
    if log_file:
        log_file.write(f"Executing: {command}\n")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        text=True,
        executable="/bin/bash",
        cwd=cwd,
        env=env
    )

    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
        if line:
            print(line.strip())
            if log_file:
                log_file.write(line)
                log_file.flush()

    return process.returncode


def trim_pdb_by_hotspot(pdb_path, hotspot_str, radius, output_path):
    """
    핫스팟 잔기 중심에서 radius(Å) 이내의 잔기만 남긴 PDB를 output_path에 저장.
    hotspot_str 예: "A425,A426,A430"
    반환값: (성공 여부, 포함된 잔기 수)
    """
    if not hotspot_str:
        return False, 0

    hotspot_keys = set()
    for token in hotspot_str.split(","):
        token = token.strip()
        if len(token) >= 2:
            chain = token[0]
            try:
                resnum = int(token[1:])
                hotspot_keys.add((chain, resnum))
            except ValueError:
                pass

    if not hotspot_keys:
        return False, 0

    atom_lines = []
    other_lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                chain = line[21]
                try:
                    resnum = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atom_lines.append((chain, resnum, x, y, z, line))
                except ValueError:
                    other_lines.append(line)
            elif rec != "HETATM":
                other_lines.append(line)

    hotspot_coords = [
        (x, y, z)
        for chain, resnum, x, y, z, _ in atom_lines
        if (chain, resnum) in hotspot_keys
    ]

    if not hotspot_coords:
        found_keys = {(chain, resnum) for chain, resnum, *_ in atom_lines}
        missing = hotspot_keys - found_keys
        print(f"[Warning] trim_pdb_by_hotspot: 요청한 hotspot 잔기가 PDB에 존재하지 않습니다.\n"
              f"  요청: {sorted(hotspot_keys)}\n  누락: {sorted(missing)}")
        return False, 0

    def min_dist_to_hotspot(x, y, z):
        return min(
            math.sqrt((x - hx)**2 + (y - hy)**2 + (z - hz)**2)
            for hx, hy, hz in hotspot_coords
        )

    kept_residues = set()
    for chain, resnum, x, y, z, _ in atom_lines:
        if min_dist_to_hotspot(x, y, z) <= radius:
            kept_residues.add((chain, resnum))

    with open(output_path, "w") as out:
        for line in other_lines:
            if line[:6].strip() not in ("ATOM", "HETATM", "END"):
                out.write(line)
        for chain, resnum, x, y, z, line in atom_lines:
            if (chain, resnum) in kept_residues:
                out.write(line)
        out.write("END\n")

    return True, len(kept_residues)


def renumber_pdb(input_path, output_path):
    """
    PDB 잔기를 체인별로 1부터 순차 재번호 매김 (gap 제거).
    반환: {(chain, old_resnum): new_resnum} 매핑
    """
    residue_mapping = {}
    chain_counters = {}
    output_lines = []

    with open(input_path, "r") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                chain = line[21]
                try:
                    old_resnum = int(line[22:26].strip())
                except ValueError:
                    output_lines.append(line)
                    continue
                key = (chain, old_resnum)
                if key not in residue_mapping:
                    chain_counters[chain] = chain_counters.get(chain, 0) + 1
                    residue_mapping[key] = chain_counters[chain]
                new_resnum = residue_mapping[key]
                new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                output_lines.append(new_line)
            elif rec != "HETATM":
                output_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(output_lines)

    return residue_mapping


def remap_hotspot(hotspot_str, residue_mapping):
    """
    재번호 매핑을 이용해 핫스팟 잔기 번호 업데이트.
    예: 'A425,A426' → 'A3,A4'
    """
    if not hotspot_str:
        return hotspot_str
    new_tokens = []
    for token in hotspot_str.split(","):
        token = token.strip()
        if len(token) >= 2:
            chain = token[0]
            try:
                resnum = int(token[1:])
                new_resnum = residue_mapping.get((chain, resnum), resnum)
                new_tokens.append(f"{chain}{new_resnum}")
            except ValueError:
                new_tokens.append(token)
        else:
            new_tokens.append(token)
    return ",".join(new_tokens)


def update_contig_after_renumber(contigs_str, residue_mapping):
    """
    재번호 매핑 기반으로 contig의 체인 범위를 업데이트.
    재번호 이후 잔기는 1..N으로 연속적이므로 정확한 범위 사용 가능.
    숫자만 있는 부분(바인더 길이 등)은 그대로 유지.
    """
    chain_ranges = {}
    for (chain, _old), new_res in residue_mapping.items():
        chain_ranges.setdefault(chain, []).append(new_res)

    def replace_chain_range(match):
        chain = match.group(1)
        if chain in chain_ranges:
            resnums = sorted(chain_ranges[chain])
            return f"{chain}{resnums[0]}-{resnums[-1]}"
        return match.group(0)

    return re.sub(r'([A-Z])(\d+)-(\d+)', replace_chain_range, contigs_str)


def _run_rfdiff_batched(num_designs, batch_size, rfdiff_opts_base,
                         output_dir, name, venv_activate, rfdiff_dir,
                         log_file):
    """num_designs를 batch_size씩 순차 배치, 기존 PDB가 있으면 자동 재개."""
    # 기존 PDB 파일 스캔 → 이미 완료된 디자인 수 계산
    existing = glob.glob(os.path.join(output_dir, f"{name}_*.pdb"))
    done_indices = []
    for e in existing:
        m = re.match(r".*_(\d+)\.pdb$", e)
        if m:
            done_indices.append(int(m.group(1)))

    if done_indices:
        resume_from = max(done_indices) + 1
    else:
        resume_from = 0

    if resume_from >= num_designs:
        skip_msg = (f"[Step 1] 이미 {resume_from}개 디자인 완료 "
                    f"(목표: {num_designs}개). Step 1 건너뜀.\n")
        print(skip_msg)
        log_file.write(skip_msg)
        return 0

    if resume_from > 0:
        resume_msg = (f"[Step 1] 기존 디자인 {resume_from}개 감지. "
                      f"design {resume_from}부터 재개합니다.\n")
        print(resume_msg)
        log_file.write(resume_msg)

    remaining = num_designs - resume_from
    num_batches = math.ceil(remaining / batch_size)

    for i in range(num_batches):
        start_idx = resume_from + i * batch_size
        batch_n = min(batch_size, num_designs - start_idx)

        opts = list(rfdiff_opts_base)
        opts.append(f"inference.num_designs={batch_n}")
        opts.append(f"inference.design_startnum={start_idx}")

        batch_msg = (f"\n[Step 1] Batch {i+1}/{num_batches}: "
                     f"designs {start_idx}~{start_idx + batch_n - 1} "
                     f"({batch_n}개)\n")
        print(batch_msg)
        log_file.write(batch_msg)

        cmd = (f"{venv_activate} && "
               f"export PYTHONPATH={rfdiff_dir}:$PYTHONPATH && "
               f"python {rfdiff_dir}/scripts/run_inference.py "
               f"{' '.join(opts)}")

        rc = run_command(cmd, log_file=log_file, cwd=rfdiff_dir)
        if rc != 0:
            fail_msg = (f"[Step 1] Batch {i+1} 실패 (design {start_idx}~). "
                        f"생성 완료된 디자인은 보존됩니다. "
                        f"같은 프로젝트명으로 재실행하면 자동 재개됩니다.\n")
            print(fail_msg)
            log_file.write(fail_msg)
            return rc

    return 0


def main():
    # GPU 메모리 사용량 제한
    if torch.cuda.is_available():
        torch.cuda.set_per_process_memory_fraction(GPU_MEMORY_FRACTION, device=0)
        print(f"[GPU] VRAM 사용량 {GPU_MEMORY_FRACTION*100:.0f}%로 제한 설정 완료")

    parser = argparse.ArgumentParser(description="RFdiffusion Pipeline v4 (Docker CLI)")
    parser.add_argument("--name", type=str, default="test", help="Job folder name")
    parser.add_argument("--contigs", type=str, default="100", help="Contig map (e.g. '100' or 'A:50')")
    parser.add_argument("--num_designs", type=int, default=5, help="Number of designs to generate (5~10 권장)")
    parser.add_argument("--iterations", type=int, default=50, help="Number of diffusion iterations (T)")
    parser.add_argument("--num_seqs", type=int, default=50, help="Number of sequences per design for ProteinMPNN")
    parser.add_argument("--num_recycles", type=int, default=3, help="Number of AlphaFold recycles")
    parser.add_argument("--input_pdb", type=str, default="", help="Input PDB file path or filename (for binders/motifs). Filename-only is resolved under RFDIFF_INPUT_DIR.")
    parser.add_argument("--hotspot", type=str, default="", help="Hotspot residues (e.g. 'A30,A33')")
    parser.add_argument("--hotspot_radius", type=float, default=10.0, help="Trimming radius around hotspot (Å)")
    parser.add_argument("--use_multimer", action="store_true", help="Use AlphaFold-multimer for validation")
    parser.add_argument("--initial_guess", action="store_true", help="mk_af_model의 initial_guess 플래그 (Step 2/3)")
    parser.add_argument("--copies", type=int, default=1, help="Homooligomer copies (Step 2/3, partial/fixbb 프로토콜)")
    parser.add_argument("--param_dir", type=str, default=AF2_PARAM_DIR, help="AlphaFold parameter files directory")
    parser.add_argument("--step", type=str, default="all",
                        choices=["0", "1", "2", "3", "all"],
                        help="실행할 단계 지정: '0'=트리밍, '1'=RFdiffusion, "
                             "'2'=ProteinMPNN, '3'=AF2(-multimer), 'all'=전체")
    parser.add_argument("--backends_dir", type=str,
                        default=os.environ.get("BACKENDS_DIR", "/app/backends"),
                        help="Directory containing run_mpnn.py / run_af2.py")

    args = parser.parse_args()

    # --input_pdb 경로 해석: 파일명만 입력 시 RFDIFF_INPUT_DIR에서 찾기
    if args.input_pdb:
        if os.path.isabs(args.input_pdb):
            resolved_input_pdb = args.input_pdb
        else:
            resolved_input_pdb = os.path.join(RFDIFF_INPUT_DIR, args.input_pdb)
    else:
        resolved_input_pdb = ""

    base_dir = RFDIFF_BASE_DIR
    rfdiff_dir = RFDIFF_CODE_DIR
    colabdesign_dir = COLABDESIGN_DIR
    af2_cache_dir = AF2_CACHE_DIR
    output_dir = os.path.join(RFDIFF_OUTPUT_DIR, args.name)
    os.makedirs(output_dir, exist_ok=True)

    log_path = os.path.join(output_dir, "pipeline.log")
    run_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    venv_activate = _build_venv_activate()

    pipeline_start = time.time()
    with open(log_path, "a", encoding="utf-8") as log_file:
        log_header = (
            f"\n{'='*60}\n"
            f"Pipeline Run: {run_time}\n"
            f"Job: {args.name} | PDB: {resolved_input_pdb} | Contigs: {args.contigs}\n"
            f"Hotspot: {args.hotspot} | Radius: {args.hotspot_radius}Å\n"
            f"Designs: {args.num_designs} | Iters(T): {args.iterations} | "
            f"Seqs: {args.num_seqs} | Recycles: {args.num_recycles}\n"
            f"{'='*60}\n"
        )
        print(log_header)
        log_file.write(log_header)

        # 파일 존재 여부 확인
        if resolved_input_pdb:
            if not os.path.exists(resolved_input_pdb):
                err_msg = f"[Error] 파일 없음: {resolved_input_pdb}\n"
                print(err_msg)
                log_file.write(err_msg)
                sys.exit(1)
            else:
                ok_msg = f"파일 확인됨: {resolved_input_pdb}\n"
                print(ok_msg)
                log_file.write(ok_msg)

        state_file = os.path.join(output_dir, "step_state.json")

        # ── Step 0: 트리밍 ──────────────────────────────────────────────────
        input_pdb = resolved_input_pdb
        rfdiff_contigs = args.contigs
        rfdiff_hotspot = args.hotspot
        step0_elapsed = 0.0

        if args.step in ("0", "all"):
            if resolved_input_pdb and args.hotspot:
                step0_header = f"\n--- [Step 0] Trimming PDB to hotspot ±{args.hotspot_radius}Å ---\n"
                print(step0_header)
                log_file.write(step0_header)
                step0_start = time.time()

                trimmed_pdb = os.path.join(output_dir, f"{args.name}_trimmed.pdb")
                ok, n_res = trim_pdb_by_hotspot(resolved_input_pdb, args.hotspot, args.hotspot_radius, trimmed_pdb)

                if ok:
                    renum_pdb = os.path.join(output_dir, f"{args.name}_trimmed_renum.pdb")
                    res_mapping = renumber_pdb(trimmed_pdb, renum_pdb)

                    input_pdb = renum_pdb
                    rfdiff_contigs = update_contig_after_renumber(args.contigs, res_mapping)
                    rfdiff_hotspot = remap_hotspot(args.hotspot, res_mapping)

                    step0_elapsed = time.time() - step0_start
                    msg = (
                        f"[Step 0] Trimming 완료: {n_res}개 잔기 유지 "
                        f"(반경 {args.hotspot_radius}Å) → {trimmed_pdb}\n"
                        f"         재번호 매김 완료 → {renum_pdb}\n"
                        f"         Contig:  {args.contigs} → {rfdiff_contigs}\n"
                        f"         Hotspot: {args.hotspot} → {rfdiff_hotspot}\n"
                        f"         소요: {step0_elapsed:.1f}초\n"
                    )
                else:
                    step0_elapsed = time.time() - step0_start
                    msg = f"[Step 0] Trimming 건너뜀 (핫스팟 파싱 실패 또는 매칭 없음) - 원본 PDB 사용\n"
                print(msg)
                log_file.write(msg)

            with open(state_file, "w", encoding="utf-8") as sf:
                json.dump({
                    "input_pdb": input_pdb,
                    "rfdiff_contigs": rfdiff_contigs,
                    "rfdiff_hotspot": rfdiff_hotspot,
                }, sf, ensure_ascii=False)

            if args.step == "0":
                done_msg = (
                    f"\n--- [Step 0] 완료 ---\n"
                    f"[Step 0] Trimming: {step0_elapsed:.1f}초\n"
                    f"Results are in: {output_dir}\n"
                )
                print(done_msg)
                log_file.write(done_msg)
                return

        # ── Step 1: RFdiffusion ─────────────────────────────────────────────
        step1_elapsed = 0.0

        if args.step in ("1", "all"):
            if args.step == "1" and os.path.exists(state_file):
                with open(state_file, "r", encoding="utf-8") as sf:
                    state = json.load(sf)
                input_pdb = state.get("input_pdb", resolved_input_pdb)
                rfdiff_contigs = state.get("rfdiff_contigs", args.contigs)
                rfdiff_hotspot = state.get("rfdiff_hotspot", args.hotspot)
                load_msg = f"[Step 1] step_state.json 로드: input_pdb={input_pdb}, contigs={rfdiff_contigs}, hotspot={rfdiff_hotspot}\n"
                print(load_msg)
                log_file.write(load_msg)

            step1_header = "\n--- [Step 1] Running RFdiffusion ---\n"
            print(step1_header)
            log_file.write(step1_header)
            step1_start = time.time()

            BATCH_THRESHOLD = 100
            BATCH_SIZE = 50

            rfdiff_opts = [
                f"inference.output_prefix={output_dir}/{args.name}",
                f"diffuser.T={args.iterations}",
                f"'contigmap.contigs=[{rfdiff_contigs}]'"
            ]
            if input_pdb:
                rfdiff_opts.append(f"inference.input_pdb={input_pdb}")
            if rfdiff_hotspot:
                rfdiff_opts.append(f"'ppi.hotspot_res=[{rfdiff_hotspot}]'")

            if args.num_designs > BATCH_THRESHOLD:
                batch_msg = (f"[Step 1] 배치 모드: {args.num_designs}개를 "
                             f"{BATCH_SIZE}개씩 "
                             f"{math.ceil(args.num_designs/BATCH_SIZE)}배치로 분할\n")
                print(batch_msg)
                log_file.write(batch_msg)

                rc = _run_rfdiff_batched(
                    args.num_designs, BATCH_SIZE, rfdiff_opts,
                    output_dir, args.name,
                    venv_activate, rfdiff_dir, log_file
                )
            else:
                rfdiff_opts.append(f"inference.num_designs={args.num_designs}")
                cmd_rfdiff = (f"{venv_activate} && "
                              f"export PYTHONPATH={rfdiff_dir}:$PYTHONPATH && "
                              f"python {rfdiff_dir}/scripts/run_inference.py "
                              f"{' '.join(rfdiff_opts)}")
                rc = run_command(cmd_rfdiff, log_file=log_file, cwd=rfdiff_dir)
            step1_elapsed = time.time() - step1_start
            step1_msg = f"[Step 1] RFdiffusion 완료: {step1_elapsed:.1f}초 ({step1_elapsed/60:.1f}분)\n"
            print(step1_msg)
            log_file.write(step1_msg)
            if rc != 0:
                msg = f"Error in RFdiffusion step. Return code: {rc}\n"
                print(msg)
                log_file.write(msg)
                sys.exit(1)

            if args.step == "1":
                done_msg = (
                    f"\n--- [Step 1] 완료 ---\n"
                    f"[Step 1] RFdiffusion: {step1_elapsed:.1f}초 ({step1_elapsed/60:.1f}분)\n"
                    f"Results are in: {output_dir}\n"
                )
                print(done_msg)
                log_file.write(done_msg)
                return

        # state dict: loaded-or-synthesized, shared by Step 2 & Step 3.
        def _load_state():
            if os.path.exists(state_file):
                with open(state_file, "r", encoding="utf-8") as sf:
                    return json.load(sf)
            return {
                "input_pdb": input_pdb,
                "rfdiff_contigs": rfdiff_contigs,
                "rfdiff_hotspot": rfdiff_hotspot,
            }

        def _save_state(st):
            with open(state_file, "w", encoding="utf-8") as sf:
                json.dump(st, sf, ensure_ascii=False)

        # ── Step 2: ProteinMPNN ────────────────────────────────────────────
        step2_elapsed = 0.0

        if args.step in ("2", "all"):
            state = _load_state()
            if args.step == "2":
                load_msg = f"[Step 2] step_state.json 로드: {state}\n"
                print(load_msg)
                log_file.write(load_msg)

            step2_header = "\n--- [Step 2] Running ProteinMPNN ---\n"
            print(step2_header)
            log_file.write(step2_header)
            step2_start = time.time()

            mpnn_dir = os.path.join(output_dir, "mpnn")
            mpnn_opts = [
                f"--pdb={output_dir}/{args.name}_0.pdb",
                f"--loc={mpnn_dir}",
                f"--num_seqs={args.num_seqs}",
                f"--num_designs={args.num_designs}",
                f'--contigs={state["rfdiff_contigs"]}',
                f"--param_dir={args.param_dir}",
                f"--copies={args.copies}",
            ]
            if args.use_multimer:
                mpnn_opts.append("--use_multimer")
            if args.initial_guess:
                mpnn_opts.append("--initial_guess")

            cmd_mpnn = (
                f"{venv_activate} && "
                f"export PYTHONPATH=$PYTHONPATH:{colabdesign_dir} && "
                f"python {args.backends_dir}/run_mpnn.py {' '.join(mpnn_opts)}"
            )

            rc = run_command(cmd_mpnn, log_file=log_file, cwd=base_dir)
            step2_elapsed = time.time() - step2_start
            step2_msg = f"[Step 2] ProteinMPNN 완료: {step2_elapsed:.1f}초\n"
            print(step2_msg)
            log_file.write(step2_msg)
            if rc != 0:
                msg = f"Error in Step 2 (ProteinMPNN). Return code: {rc}\n"
                print(msg)
                log_file.write(msg)
                sys.exit(1)

            state["mpnn_done"] = True
            state["mpnn_fasta"] = "mpnn/design.fasta"
            _save_state(state)

            if args.step == "2":
                done_msg = (
                    f"\n--- [Step 2] 완료 ---\n"
                    f"[Step 2] ProteinMPNN: {step2_elapsed:.1f}초\n"
                    f"Results are in: {mpnn_dir}\n"
                )
                print(done_msg)
                log_file.write(done_msg)
                return

        # ── Step 3: AlphaFold2 / AF2-multimer ──────────────────────────────
        step3_elapsed = 0.0

        if args.step in ("3", "all"):
            state = _load_state()
            if args.step == "3":
                load_msg = f"[Step 3] step_state.json 로드: {state}\n"
                print(load_msg)
                log_file.write(load_msg)
                if not state.get("mpnn_done"):
                    warn = ("[Step 3] 경고: step_state.json에 mpnn_done 플래그가 "
                            "없음. Step 2를 먼저 실행했는지 확인하세요.\n")
                    print(warn)
                    log_file.write(warn)

            step3_header = "\n--- [Step 3] Running AlphaFold2 ---\n"
            print(step3_header)
            log_file.write(step3_header)
            step3_start = time.time()

            af2_dir = os.path.join(output_dir, "af2")
            fasta_rel = state.get("mpnn_fasta", "mpnn/design.fasta")
            fasta_abs = os.path.join(output_dir, fasta_rel)

            af_opts = [
                f"--pdb={output_dir}/{args.name}_0.pdb",
                f"--fasta={fasta_abs}",
                f"--loc={af2_dir}",
                f"--num_designs={args.num_designs}",
                f"--num_recycles={args.num_recycles}",
                f"--param_dir={args.param_dir}",
                f'--contigs={state["rfdiff_contigs"]}',
                f"--copies={args.copies}",
            ]
            if args.use_multimer:
                af_opts.append("--use_multimer")
            if args.initial_guess:
                af_opts.append("--initial_guess")

            cmd_af = (
                f"{venv_activate} && "
                f"export PYTHONPATH=$PYTHONPATH:{colabdesign_dir} && "
                f"export JAX_COMPILATION_CACHE_DIR={af2_cache_dir} && "
                f"python {args.backends_dir}/run_af2.py {' '.join(af_opts)}"
            )

            rc = run_command(cmd_af, log_file=log_file, cwd=base_dir)
            step3_elapsed = time.time() - step3_start
            step3_msg = f"[Step 3] AF2 완료: {step3_elapsed:.1f}초 ({step3_elapsed/60:.1f}분)\n"
            print(step3_msg)
            log_file.write(step3_msg)
            if rc != 0:
                msg = f"Error in Step 3 (AF2). Return code: {rc}\n"
                print(msg)
                log_file.write(msg)
                sys.exit(1)

        total_elapsed = time.time() - pipeline_start
        done_msg = (
            f"\n--- Pipeline Completed Successfully ---\n"
            f"[Step 0] Trimming:     {step0_elapsed:.1f}초\n"
            f"[Step 1] RFdiffusion:  {step1_elapsed:.1f}초 ({step1_elapsed/60:.1f}분)\n"
            f"[Step 2] ProteinMPNN:  {step2_elapsed:.1f}초 ({step2_elapsed/60:.1f}분)\n"
            f"[Step 3] AlphaFold2:   {step3_elapsed:.1f}초 ({step3_elapsed/60:.1f}분)\n"
            f"[Total]  전체 소요:    {total_elapsed:.1f}초 ({total_elapsed/60:.1f}분)\n"
            f"Results are in: {output_dir}\n"
            f"Log saved to: {log_path}\n"
        )
        print(done_msg)
        log_file.write(done_msg)

if __name__ == "__main__":
    main()
