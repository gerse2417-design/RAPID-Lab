"""
run_make_ari_v5.py (호스트 통합 러너 — 프로젝트 단일 폴더 + KST)
────────────────────────────────────────
v4 베이스 + 다음 사항:
  - `--project-dir` CLI 옵션 및 `INHIBIFORGE_PROJECT_DIR` 환경변수 지원.
    지정되면 output_dir 가 이 경로로 고정되고, 상대 경로 input_pdb 도 이 폴더
    기준으로 먼저 해석한다 (없으면 기존 INPUT_DIR fallback).
  - 모든 로그/상태 파일 타임스탬프를 KST(UTC+9) 로 명시.

v4 원본은 건드리지 않고 이 파일을 새로 추가.
"""
import os
import sys
import subprocess
import argparse
import time
import datetime
import math
import re
import json
import torch

# 공통 inputs/outputs 루트 헬퍼
#   컨테이너(Docker)에서는 AMR_ROOT=/app 로 주입. 호스트에선 기본값(/home/sooyeon/amr) 사용.
_AMR_ROOT = os.environ.get("AMR_ROOT", "/home/sooyeon/amr")
if _AMR_ROOT not in sys.path:
    sys.path.insert(0, _AMR_ROOT)
from common.amr_paths import input_root, output_root, ensure_roots
ensure_roots()

INPUT_DIR = str(input_root())

# Docker 이미지와 공용으로 쓰는 분리 러너 스크립트들이 있는 디렉토리.
# 기본값: amr CLI 배포 디렉토리 내부. BACKENDS_DIR 환경변수로 오버라이드 가능.
BACKENDS_DIR = os.environ.get(
    "BACKENDS_DIR",
    "/home/sooyeon/amr/deploy/inhibiforge-cli/app/backends",
)

# ── KST 타임존 헬퍼 ────────────────────────────────────────────────────────
KST = datetime.timezone(datetime.timedelta(hours=9))


def _kst_now_str() -> str:
    """KST 기준 현재 시각을 '%Y-%m-%d %H:%M:%S%z' 형식 문자열로 반환."""
    return datetime.datetime.now(KST).strftime("%Y-%m-%d %H:%M:%S%z")


def run_command(command, log_file=None, cwd=None, env=None):
    """실시간으로 출력을 확인하며 쉘 커맨드를 실행하는 유틸리티."""
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
        env=env,
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
    """핫스팟 잔기 중심 radius(Å) 이내 잔기만 남긴 PDB 저장. (v3/v4와 동일)"""
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
            math.sqrt((x - hx) ** 2 + (y - hy) ** 2 + (z - hz) ** 2)
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


def _resolve_project_override(args) -> str | None:
    """--project-dir CLI 또는 INHIBIFORGE_PROJECT_DIR env. 없으면 None."""
    return args.project_dir or os.environ.get("INHIBIFORGE_PROJECT_DIR")


def _select_best_pdb(output_dir: str, job_name: str, num_designs: int, log_file=None) -> None:
    """mpnn_results.csv 에서 최저 rmsd 행의 design index 를 찾아
    해당 best_design{m}.pdb 를 output_dir/best.pdb 로 복사한다.

    실패해도 예외를 던지지 않고 경고 메시지만 남긴다 (디버깅 편의).
    """
    import csv as _csv
    import shutil as _shutil

    def _log(msg: str) -> None:
        print(msg)
        if log_file is not None:
            try:
                log_file.write(msg if msg.endswith("\n") else msg + "\n")
            except Exception:
                pass

    csv_path = os.path.join(output_dir, "mpnn_results.csv")
    if not os.path.exists(csv_path):
        _log(f"[best.pdb] mpnn_results.csv 없음 → 선정 건너뜀 ({csv_path})")
        return

    best_design = None
    best_rmsd = float("inf")
    try:
        with open(csv_path, newline="") as f:
            reader = _csv.DictReader(f)
            for row in reader:
                try:
                    rmsd = float(row.get("rmsd", "nan"))
                    m    = int(row.get("design", -1))
                except (TypeError, ValueError):
                    continue
                if rmsd == rmsd and rmsd < best_rmsd:  # NaN-safe
                    best_rmsd = rmsd
                    best_design = m
    except Exception as e:
        _log(f"[best.pdb] CSV 파싱 실패: {e}")
        return

    # 후보 경로: 먼저 best_design{m}.pdb, fallback 은 {job_name}_{m}.pdb (RFdiffusion raw)
    candidates = []
    if best_design is not None:
        candidates.append(os.path.join(output_dir, f"best_design{best_design}.pdb"))
        candidates.append(os.path.join(output_dir, f"{job_name}_{best_design}.pdb"))
    # 최후 fallback — 0번 디자인
    candidates.append(os.path.join(output_dir, "best_design0.pdb"))
    candidates.append(os.path.join(output_dir, f"{job_name}_0.pdb"))

    src = next((p for p in candidates if os.path.exists(p)), None)
    if src is None:
        _log(f"[best.pdb] 후보 PDB 없음 — 선정 실패. 후보: {candidates}")
        return

    dst = os.path.join(output_dir, "best.pdb")
    try:
        _shutil.copy2(src, dst)
        _log(f"[best.pdb] 선정: design={best_design} rmsd={best_rmsd:.3f} "
             f"src={os.path.basename(src)} → best.pdb")
    except Exception as e:
        _log(f"[best.pdb] 복사 실패: {e}")


def main():
    if torch.cuda.is_available():
        torch.cuda.set_per_process_memory_fraction(0.7, device=0)
        print("[GPU] VRAM 사용량 70%로 제한 설정 완료")

    parser = argparse.ArgumentParser(description="RFdiffusion Pipeline v5 (host, Step 2/3 split, project-dir, KST)")
    parser.add_argument("--name", type=str, default="test", help="Job folder name")
    parser.add_argument("--contigs", type=str, default="100", help="Contig map (e.g. '100' or 'A:50')")
    parser.add_argument("--num_designs", type=int, default=5, help="Number of designs to generate (5~10 권장)")
    parser.add_argument("--iterations", type=int, default=50, help="Number of diffusion iterations (T)")
    parser.add_argument("--num_seqs", type=int, default=50, help="Number of sequences per design for ProteinMPNN")
    parser.add_argument("--num_recycles", type=int, default=3, help="Number of AlphaFold recycles")
    parser.add_argument("--input_pdb", type=str, default="", help="Input PDB file path or filename")
    parser.add_argument("--hotspot", type=str, default="", help="Hotspot residues (e.g. 'A30,A33')")
    parser.add_argument("--hotspot_radius", type=float, default=10.0, help="Trimming radius around hotspot (Å)")
    parser.add_argument("--use_multimer", action="store_true", help="Use AlphaFold-multimer for validation")
    parser.add_argument("--param_dir", type=str,
                        default=os.environ.get(
                            "AF2_PARAM_DIR",
                            "/home/sooyeon/amr/RFdiffusion/af2_cache/colabfold",
                        ),
                        help="AlphaFold parameter files directory (env AF2_PARAM_DIR)")
    parser.add_argument("--step", type=str, default="all",
                        choices=["0", "1", "2", "3", "all"],
                        help="실행할 단계: '0'=트리밍, '1'=RFdiffusion, "
                             "'2'=ProteinMPNN, '3'=AF2(-multimer), 'all'=전체")
    parser.add_argument("--backends_dir", type=str, default=BACKENDS_DIR,
                        help="Directory containing run_mpnn.py / run_af2.py")
    parser.add_argument("--project-dir", dest="project_dir", type=str, default=None,
                        help="프로젝트 폴더 경로. 지정하면 output_dir 가 이 경로로 고정되고, "
                             "상대 경로의 input PDB 도 이 폴더 기준으로 먼저 해석한다. "
                             "(INHIBIFORGE_PROJECT_DIR env 로도 지정 가능)")

    args = parser.parse_args()

    project_override = _resolve_project_override(args)

    if args.input_pdb:
        if os.path.isabs(args.input_pdb):
            resolved_input_pdb = args.input_pdb
        else:
            # project_override 가 있으면 프로젝트 폴더에서 먼저 찾고, 없으면 기존 INPUT_DIR fallback
            candidate = None
            if project_override:
                candidate = os.path.join(project_override, args.input_pdb)
                if not os.path.exists(candidate):
                    candidate = None
            if candidate:
                resolved_input_pdb = candidate
            else:
                resolved_input_pdb = os.path.join(INPUT_DIR, args.input_pdb)
    else:
        resolved_input_pdb = ""

    # 컨테이너(Docker)에서는 RFDIFF_BASE_DIR=/mnt/ebs/rfdiffusion 로 주입.
    base_dir = os.environ.get("RFDIFF_BASE_DIR", "/home/sooyeon/amr/RFdiffusion")
    rfdiff_dir = os.path.join(base_dir, "RFdiffusion")
    colabdesign_dir = os.path.join(base_dir, "ColabDesign_repo")
    af2_cache_dir = os.path.join(base_dir, "af2_cache")

    if project_override:
        output_dir = os.path.abspath(project_override)
    else:
        output_dir = str(output_root() / args.name)
    os.makedirs(output_dir, exist_ok=True)

    log_path = os.path.join(output_dir, "pipeline.log")
    run_time = _kst_now_str()

    # venv 활성화 + CUDA 경로 세팅.
    #   - 호스트: 기본값(아래 _default_venv) 사용 — 기존 v3/v4 호환
    #   - 컨테이너: ENV `VENV_CMD=""` 로 주입 → 아래 `or ":"` 로 no-op 쉘 명령(`:`)으로 치환
    #     (뒤에 `&& python ...` 이 붙으므로 공백이면 쉘 문법 에러, 반드시 유효한 명령이어야 함)
    _default_venv = (
        "source /home/sooyeon/amr/bin/activate && "
        "export LD_LIBRARY_PATH=/usr/local/cuda-12.6/lib64:/home/sooyeon/amr/lib:$LD_LIBRARY_PATH && "
        "export PATH=/usr/local/cuda-12.6/bin:$PATH && "
        "export E3NN_JIT_COMPILE=0"
    )
    venv_activate = os.environ.get("VENV_CMD", _default_venv) or ":"

    pipeline_start = time.time()
    with open(log_path, "a", encoding="utf-8") as log_file:
        log_header = (
            f"\n{'='*60}\n"
            f"Pipeline Run (v5 host, KST): {run_time}\n"
            f"Job: {args.name} | PDB: {resolved_input_pdb} | Contigs: {args.contigs}\n"
            f"Hotspot: {args.hotspot} | Radius: {args.hotspot_radius}Å\n"
            f"Designs: {args.num_designs} | Iters(T): {args.iterations} | "
            f"Seqs: {args.num_seqs} | Recycles: {args.num_recycles}\n"
            f"Output dir: {output_dir}\n"
            f"{'='*60}\n"
        )
        print(log_header)
        log_file.write(log_header)

        if resolved_input_pdb:
            if not os.path.exists(resolved_input_pdb):
                err_msg = f"[Error] 파일 없음: {resolved_input_pdb}\n"
                print(err_msg); log_file.write(err_msg); sys.exit(1)
            else:
                ok_msg = f"파일 확인됨: {resolved_input_pdb}\n"
                print(ok_msg); log_file.write(ok_msg)

        state_file = os.path.join(output_dir, "step_state.json")

        # ── Step 0 ─────────────────────────────────────────────────────────
        input_pdb = resolved_input_pdb
        rfdiff_contigs = args.contigs
        rfdiff_hotspot = args.hotspot
        step0_elapsed = 0.0

        if args.step in ("0", "all"):
            if resolved_input_pdb and args.hotspot:
                step0_header = f"\n--- [Step 0] Trimming PDB to hotspot ±{args.hotspot_radius}Å ---\n"
                print(step0_header); log_file.write(step0_header)
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
                print(msg); log_file.write(msg)

            with open(state_file, "w", encoding="utf-8") as sf:
                json.dump({
                    "input_pdb": input_pdb,
                    "rfdiff_contigs": rfdiff_contigs,
                    "rfdiff_hotspot": rfdiff_hotspot,
                    "updated_at": _kst_now_str(),
                }, sf, ensure_ascii=False)

            if args.step == "0":
                done_msg = (
                    f"\n--- [Step 0] 완료 ---\n"
                    f"[Step 0] Trimming: {step0_elapsed:.1f}초\n"
                    f"Results are in: {output_dir}\n"
                )
                print(done_msg); log_file.write(done_msg); return

        # ── Step 1: RFdiffusion ────────────────────────────────────────────
        step1_elapsed = 0.0

        if args.step in ("1", "all"):
            if args.step == "1" and os.path.exists(state_file):
                with open(state_file, "r", encoding="utf-8") as sf:
                    state = json.load(sf)
                input_pdb = state.get("input_pdb", resolved_input_pdb)
                rfdiff_contigs = state.get("rfdiff_contigs", args.contigs)
                rfdiff_hotspot = state.get("rfdiff_hotspot", args.hotspot)
                load_msg = f"[Step 1] step_state.json 로드: {state}\n"
                print(load_msg); log_file.write(load_msg)

            step1_header = "\n--- [Step 1] Running RFdiffusion ---\n"
            print(step1_header); log_file.write(step1_header)
            step1_start = time.time()

            rfdiff_opts = [
                f"inference.output_prefix={output_dir}/{args.name}",
                f"inference.num_designs={args.num_designs}",
                f"diffuser.T={args.iterations}",
                f"'contigmap.contigs=[{rfdiff_contigs}]'",
                f"hydra.run.dir={output_dir}/hydra",
            ]
            if input_pdb:
                rfdiff_opts.append(f"inference.input_pdb={input_pdb}")
            if rfdiff_hotspot:
                rfdiff_opts.append(f"'ppi.hotspot_res=[{rfdiff_hotspot}]'")

            cmd_rfdiff = f"{venv_activate} && python {rfdiff_dir}/scripts/run_inference.py {' '.join(rfdiff_opts)}"

            rc = run_command(cmd_rfdiff, log_file=log_file, cwd=rfdiff_dir)
            step1_elapsed = time.time() - step1_start
            step1_msg = f"[Step 1] RFdiffusion 완료: {step1_elapsed:.1f}초 ({step1_elapsed/60:.1f}분)\n"
            print(step1_msg); log_file.write(step1_msg)
            if rc != 0:
                msg = f"Error in RFdiffusion step. Return code: {rc}\n"
                print(msg); log_file.write(msg); sys.exit(1)

            if args.step == "1":
                done_msg = (
                    f"\n--- [Step 1] 완료 ---\n"
                    f"[Step 1] RFdiffusion: {step1_elapsed:.1f}초 ({step1_elapsed/60:.1f}분)\n"
                    f"Results are in: {output_dir}\n"
                )
                print(done_msg); log_file.write(done_msg); return

        # state 로딩/저장 헬퍼 (Step 2/3 공용)
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
            st["updated_at"] = _kst_now_str()
            with open(state_file, "w", encoding="utf-8") as sf:
                json.dump(st, sf, ensure_ascii=False)

        # ── Step 2: ProteinMPNN ────────────────────────────────────────────
        step2_elapsed = 0.0

        if args.step in ("2", "all"):
            state = _load_state()
            if args.step == "2":
                load_msg = f"[Step 2] step_state.json 로드: {state}\n"
                print(load_msg); log_file.write(load_msg)

            step2_header = "\n--- [Step 2] Running ProteinMPNN ---\n"
            print(step2_header); log_file.write(step2_header)
            step2_start = time.time()

            # v5: design.fasta / mpnn_results.csv / best_design*.pdb 를 프로젝트 루트로.
            mpnn_dir = output_dir
            # run_mpnn.py 는 --contigs / --param_dir 을 필수로 요구한다.
            mpnn_contigs = state.get("rfdiff_contigs", rfdiff_contigs)
            mpnn_opts = [
                f"--pdb={output_dir}/{args.name}_0.pdb",
                f"--loc={mpnn_dir}",
                f"--contigs='{mpnn_contigs}'",
                f"--param_dir={args.param_dir}",
                f"--num_seqs={args.num_seqs}",
                f"--num_designs={args.num_designs}",
            ]
            if args.use_multimer:
                mpnn_opts.append("--use_multimer")

            cmd_mpnn = (
                f"{venv_activate} && "
                f"export PYTHONPATH=$PYTHONPATH:{colabdesign_dir} && "
                f"python {args.backends_dir}/run_mpnn.py {' '.join(mpnn_opts)}"
            )

            rc = run_command(cmd_mpnn, log_file=log_file, cwd=base_dir)
            step2_elapsed = time.time() - step2_start
            step2_msg = f"[Step 2] ProteinMPNN 완료: {step2_elapsed:.1f}초\n"
            print(step2_msg); log_file.write(step2_msg)
            if rc != 0:
                msg = f"Error in Step 2 (ProteinMPNN). Return code: {rc}\n"
                print(msg); log_file.write(msg); sys.exit(1)

            state["mpnn_done"] = True
            state["mpnn_fasta"] = "design.fasta"
            _save_state(state)

            if args.step == "2":
                done_msg = (
                    f"\n--- [Step 2] 완료 ---\n"
                    f"[Step 2] ProteinMPNN: {step2_elapsed:.1f}초\n"
                    f"Results are in: {mpnn_dir}\n"
                )
                print(done_msg); log_file.write(done_msg); return

        # ── Step 3: AlphaFold2 ─────────────────────────────────────────────
        step3_elapsed = 0.0

        if args.step in ("3", "all"):
            state = _load_state()
            if args.step == "3":
                load_msg = f"[Step 3] step_state.json 로드: {state}\n"
                print(load_msg); log_file.write(load_msg)
                if not state.get("mpnn_done"):
                    warn = ("[Step 3] 경고: step_state.json에 mpnn_done 플래그가 "
                            "없음. Step 2를 먼저 실행했는지 확인하세요.\n")
                    print(warn); log_file.write(warn)

            step3_header = "\n--- [Step 3] Running AlphaFold2 ---\n"
            print(step3_header); log_file.write(step3_header)
            step3_start = time.time()

            # v5: AF2 결과도 프로젝트 루트에. all_pdb/ 서브디렉터리는 AF2 내부 생성.
            af2_dir = output_dir
            fasta_rel = state.get("mpnn_fasta", "design.fasta")
            fasta_abs = os.path.join(output_dir, fasta_rel)

            # run_af2.py 도 --contigs / --param_dir 을 필수로 요구한다.
            af_contigs = state.get("rfdiff_contigs", rfdiff_contigs)
            af_opts = [
                f"--pdb={output_dir}/{args.name}_0.pdb",
                f"--fasta={fasta_abs}",
                f"--loc={af2_dir}",
                f"--contigs='{af_contigs}'",
                f"--num_designs={args.num_designs}",
                f"--num_recycles={args.num_recycles}",
                f"--param_dir={args.param_dir}",
            ]
            if args.use_multimer:
                af_opts.append("--use_multimer")

            # JAX compilation cache 경로. 기본은 af2_cache_dir (= 호스트 공용) 이지만,
            # 컨테이너처럼 CUDA/cuDNN 빌드가 다른 환경에서는 호스트 cache 를 재사용하면
            # CUDA_ERROR_INVALID_VALUE 같은 graph 에러가 발생하므로 AF2_JAX_CACHE_DIR
            # env 로 컨테이너 전용 경로(/tmp/jax_cache 등)를 지정 가능.
            jax_cache = os.environ.get("AF2_JAX_CACHE_DIR", af2_cache_dir)
            cmd_af = (
                f"{venv_activate} && "
                f"export PYTHONPATH=$PYTHONPATH:{colabdesign_dir} && "
                f"export JAX_COMPILATION_CACHE_DIR={jax_cache} && "
                f"python {args.backends_dir}/run_af2.py {' '.join(af_opts)}"
            )

            rc = run_command(cmd_af, log_file=log_file, cwd=base_dir)
            step3_elapsed = time.time() - step3_start
            step3_msg = f"[Step 3] AF2 완료: {step3_elapsed:.1f}초 ({step3_elapsed/60:.1f}분)\n"
            print(step3_msg); log_file.write(step3_msg)
            if rc != 0:
                msg = f"Error in Step 3 (AF2). Return code: {rc}\n"
                print(msg); log_file.write(msg); sys.exit(1)

            # ── best.pdb 선정: mpnn_results.csv 에서 최저 rmsd 행의 design 으로부터 복사 ──
            # 이후 단계(체인 분리, DNAWorks 등)가 <output_dir>/best.pdb 를 읽는다.
            _select_best_pdb(output_dir, args.name, args.num_designs, log_file)

        total_elapsed = time.time() - pipeline_start
        done_msg = (
            f"\n--- Pipeline Completed Successfully ---\n"
            f"Finished at: {_kst_now_str()}\n"
            f"[Step 0] Trimming:     {step0_elapsed:.1f}초\n"
            f"[Step 1] RFdiffusion:  {step1_elapsed:.1f}초 ({step1_elapsed/60:.1f}분)\n"
            f"[Step 2] ProteinMPNN:  {step2_elapsed:.1f}초 ({step2_elapsed/60:.1f}분)\n"
            f"[Step 3] AlphaFold2:   {step3_elapsed:.1f}초 ({step3_elapsed/60:.1f}분)\n"
            f"[Total]  전체 소요:    {total_elapsed:.1f}초 ({total_elapsed/60:.1f}분)\n"
            f"Results are in: {output_dir}\n"
            f"Log saved to: {log_path}\n"
        )
        print(done_msg); log_file.write(done_msg)


if __name__ == "__main__":
    main()
