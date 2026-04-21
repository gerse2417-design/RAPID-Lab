"""
amr_dock_backend/pipeline.py
──────────────────────────────────────────────────────────────
LightDock 도킹 파이프라인 오케스트레이터 (Streamlit 의존 없음).

책임 범위:
  - BSAS 클러스터링 실행 (lgd_cluster_bsas.py, 스왐 루프)
  - LightDock 전체 파이프라인 순차 실행
      Step 2: Setup (lightdock3_setup.py)
      Step 3: GSO 시뮬레이션 (lightdock3.py)
      Step 4: 포즈 추출 (lgd_generate_conformations.py) + BSAS 클러스터링
      Step 5: 랭킹 생성 (lgd_rank.py)
      Step 5.5: fnat 필터링 (lgd_filter_restraints.py, restraints 있을 때)
  - Resume 모드 (기존 산출물 감지 → 완료된 단계 건너뜀)

콜백 인터페이스 (Streamlit 없이도 동작하도록 의존성 역전):
  - progress_cb(step_id, status, duration)  → 단계 진행 상황 알림
  - log_cb(line)                            → subprocess 실시간 출력
  - warn_cb(msg)                            → 비치명적 경고 메시지

의존성:
  - 외부: subprocess, json, shutil, datetime, re
  - 내부: amr_dock_backend.pdb_utils, restraints, swarm_utils, subprocess_utils, env
"""

import datetime
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from .env import get_optimal_cores
from .pdb_utils import preprocess_pdb_for_docking
from .restraints import generate_restraints_list
from .subprocess_utils import format_duration_min, now_kst, run_subprocess, save_timing
from .swarm_utils import iter_swarm_dirs


# ============================================================================
# Resume 모드 판정 헬퍼 (private)
# ============================================================================

def _has_pose_outputs(job_dir: Path) -> bool:
    """적어도 하나의 스왐에 lightdock_*.pdb 파일이 있으면 True."""
    for sd in job_dir.glob("swarm_*"):
        if any(sd.glob("lightdock_*.pdb")):
            return True
    return False


def _step3_done(job_dir: Path) -> bool:
    """모든 swarm_N 디렉토리에 gso_100.out이 있으면 True."""
    sds = list(job_dir.glob("swarm_*"))
    if not sds:
        return False
    return all((sd / "gso_100.out").exists() for sd in sds)


def _step4_done(job_dir: Path) -> bool:
    """모든 swarm_N 디렉토리에 lightdock_*.pdb가 있으면 True."""
    sds = list(job_dir.glob("swarm_*"))
    if not sds:
        return False
    return all(any(sd.glob("lightdock_*.pdb")) for sd in sds)


def _clustering_done(job_dir: Path) -> bool:
    """모든 swarm_N 디렉토리에 cluster.repr이 있으면 True."""
    sds = list(job_dir.glob("swarm_*"))
    if not sds:
        return False
    return all((sd / "cluster.repr").exists() for sd in sds)


def _parse_setup_json(job_dir: Path) -> tuple:
    """setup.json에서 (receptor_filename, ligand_filename)을 추출한다."""
    setup_file = job_dir / "setup.json"
    with open(setup_file, "r") as f:
        data = json.load(f)
    if "receptor_pdb" in data and "ligand_pdb" in data:
        return data["receptor_pdb"], data["ligand_pdb"]
    return data["receptor"][0], data["ligands"][0][0]


# ============================================================================
# BSAS 클러스터링
# ============================================================================

def run_bsas_clustering(
    job_dir: Path,
    rmsd_cutoff: float = 4.0,
    progress_cb=None,
) -> dict:
    """
    모든 swarm_N 디렉토리에서 lgd_cluster_bsas.py를 실행한다.

    Parameters
    ----------
    job_dir      : LightDock 작업 디렉토리
    rmsd_cutoff  : 클러스터링 RMSD 기준 (Å). lgd_cluster_bsas.py가 -c 플래그를
                   지원하는 경우에만 적용된다.
    progress_cb  : (idx: int, total: int, swarm_name: str) -> None

    Returns
    -------
    dict: {"total": int, "failed": list[str], "skipped": list[str]}
    """
    job_dir = Path(job_dir)
    swarm_dirs = iter_swarm_dirs(job_dir)
    total = len(swarm_dirs)
    failed = []
    skipped = []

    # -c 플래그 지원 여부를 루프 밖에서 1회만 확인
    probe = subprocess.run(
        ["lgd_cluster_bsas.py", "--help"],
        capture_output=True, text=True,
    )
    use_cutoff_flag = "-c" in (probe.stdout + probe.stderr)

    for idx, swarm_dir in enumerate(swarm_dirs):
        gso_files = sorted(
            swarm_dir.glob("gso_*.out"),
            key=lambda f: int(re.search(r'gso_(\d+)', f.name).group(1))
            if re.search(r'gso_(\d+)', f.name) else 0,
        )
        if not gso_files:
            skipped.append(swarm_dir.name)
            if progress_cb:
                progress_cb(idx, total, swarm_dir.name)
            continue

        gso_filename = gso_files[-1].name
        cmd = ["lgd_cluster_bsas.py", gso_filename]
        if use_cutoff_flag:
            cmd += ["-c", str(rmsd_cutoff)]

        log_path = job_dir / f"cluster_{swarm_dir.name}.log"
        rc = run_subprocess(cmd, cwd=swarm_dir, log_path=log_path)
        if rc != 0:
            failed.append(swarm_dir.name)

        if progress_cb:
            progress_cb(idx, total, swarm_dir.name)

    return {"total": total, "failed": failed, "skipped": skipped}


# ============================================================================
# 전체 파이프라인
# ============================================================================

def run_full_docking_pipeline(
    job_dir: Path,
    inputs_dir: Path,
    rec_name: str,
    lig_name: str,
    *,
    pdb_cleaning: bool = True,
    anm_option: bool = False,
    substitute_nonstandard: bool = False,
    use_pdbfixer: bool = False,
    use_restraints: bool = False,
    hotspot_input: str = "",
    hotspot_role: str = "R",
    is_aws: bool = False,
    progress_cb=None,
    log_cb=None,
    warn_cb=None,
) -> dict:
    """
    LightDock 전체 파이프라인을 순차 실행한다.

    Step 2: Setup
    Step 3: GSO 시뮬레이션
    Step 4: 포즈 추출 + BSAS 클러스터링
    Step 5: 랭킹 생성
    Step 5.5: fnat 필터링 (restraints 있을 때)

    Parameters
    ----------
    job_dir      : 작업 결과 디렉토리
    inputs_dir   : 입력 PDB 파일이 있는 디렉토리
    rec_name     : 수용체 파일명 (inputs_dir 기준)
    lig_name     : 리간드 파일명 (inputs_dir 기준)
    pdb_cleaning : True이면 HETATM 제거 클리닝 수행
    anm_option   : True이면 -anm 플래그 추가
    substitute_nonstandard: True이면 비표준 아미노산 치환 수행
    use_pdbfixer : True이면 PDBFixer 누락 원자 보완 수행
    use_restraints: True이면 hotspot_input 기반 restraints.list 생성
    hotspot_input: 핫스팟 원자 번호 문자열 (예: "A425, A426")
    hotspot_role : "R" (수용체) 또는 "L" (리간드)
    is_aws       : True이면 AWS 최적 코어 수 사용
    progress_cb  : (step_id: str, status: str, duration: str|None) -> None
    log_cb       : (line: str) -> None  (subprocess 실시간 출력)
    warn_cb      : (msg: str) -> None   (비치명적 경고)

    Returns
    -------
    dict: {
      "success": bool,
      "steps":   {"step2": "done"|"skip"|"error", ...},
      "warnings": list[str],
      "errors":   list[str],
    }
    """
    job_dir   = Path(job_dir)
    inputs_dir = Path(inputs_dir)
    job_dir.mkdir(parents=True, exist_ok=True)

    result = {"success": False, "steps": {}, "warnings": [], "errors": []}

    def _progress(step_id, status, duration=None):
        result["steps"][step_id] = status
        if progress_cb:
            progress_cb(step_id, status, duration)

    def _warn(msg):
        result["warnings"].append(msg)
        if warn_cb:
            warn_cb(msg)

    resume_mode = (job_dir / "setup.json").exists() and _has_pose_outputs(job_dir)

    # ── STEP 2: Setup ────────────────────────────────────────────────────────
    _start2 = now_kst()

    if resume_mode and (job_dir / "setup.json").exists():
        _progress("step2", "skip")
    else:
        _progress("step2", "running")

        # 기존 산출물 정리 (fresh run)
        for f in job_dir.glob("lightdock_*.pdb"):
            f.unlink()
        for d in [job_dir / "init", *job_dir.glob("swarm_*")]:
            if d.exists():
                shutil.rmtree(d)
        old_setup = job_dir / "setup.json"
        if old_setup.exists():
            old_setup.unlink()

        # PDB 파일 준비 (클리닝 or 단순 복사)
        actual_rec, actual_lig = rec_name, lig_name
        try:
            if pdb_cleaning:
                for label, src_name, dst_name in [
                    ("receptor", rec_name, f"cleaned_{Path(rec_name).stem}.pdb"),
                    ("ligand",   lig_name, f"cleaned_{Path(lig_name).stem}.pdb"),
                ]:
                    summary = preprocess_pdb_for_docking(
                        input_path=inputs_dir / src_name,
                        output_path=job_dir / dst_name,
                        substitute_nonstandard=substitute_nonstandard,
                        run_pdbfixer=use_pdbfixer,
                        strict_hetatm_removal=True,
                        logger=None,
                    )
                    for w in summary.get("warnings", []):
                        _warn(f"{label}: {w}")
                actual_rec = f"cleaned_{Path(rec_name).stem}.pdb"
                actual_lig = f"cleaned_{Path(lig_name).stem}.pdb"
            else:
                rec_dst = Path(rec_name).stem + ".pdb"
                lig_dst = Path(lig_name).stem + ".pdb"
                shutil.copy2(inputs_dir / rec_name, job_dir / rec_dst)
                shutil.copy2(inputs_dir / lig_name, job_dir / lig_dst)
                actual_rec, actual_lig = rec_dst, lig_dst
        except Exception as e:
            result["errors"].append(f"PDB 파일 준비 실패: {e}")
            _progress("step2", "error")
            return result

        # Restraints 생성
        restraints_arg = []
        if use_restraints and hotspot_input.strip():
            restraints_list_path = job_dir / "restraints.list"
            rec_pdb_path = job_dir / actual_rec
            try:
                lines_written, not_found = generate_restraints_list(
                    receptor_pdb=rec_pdb_path,
                    raw_input=hotspot_input,
                    output_path=restraints_list_path,
                    role=hotspot_role,
                )
                if not_found:
                    _warn("PDB에서 찾지 못한 원자: " +
                          ", ".join(f"{c}{n}" for c, n in not_found))
                if lines_written:
                    restraints_arg = ["-rst", "restraints.list"]
                else:
                    _warn("restraints.list에 기록할 잔기가 없습니다. 원자 번호를 확인하세요.")
            except ValueError as e:
                result["errors"].append(f"원자 번호 파싱 실패: {e}")
                _progress("step2", "error")
                return result

        # ANM ProDy 사전 체크
        if anm_option:
            try:
                import prody  # noqa: F401
            except ImportError:
                result["errors"].append(
                    "ANM 모드에 필요한 ProDy가 설치되지 않았습니다. "
                    "`conda install -c conda-forge prody`"
                )
                _progress("step2", "error")
                return result

        # lightdock3_setup.py 실행
        cmd2 = ["lightdock3_setup.py", actual_rec, actual_lig, "--noxt", "--noh", "--now"]
        if anm_option:
            cmd2.append("-anm")
        cmd2.extend(restraints_arg)

        rc2 = run_subprocess(cmd2, cwd=job_dir,
                             log_path=job_dir / "step2_setup.log",
                             on_line=log_cb)
        _end2 = now_kst()
        save_timing(job_dir, "step2", _start2, _end2)
        dur2 = format_duration_min(f"{(_end2 - _start2).total_seconds()}s")

        if not (job_dir / "setup.json").exists():
            result["errors"].append("Setup 실패. step2_setup.log를 확인하세요.")
            _progress("step2", "error", dur2)
            return result
        _progress("step2", "done", dur2)

    # ── STEP 3: GSO 시뮬레이션 ───────────────────────────────────────────────
    _start3 = now_kst()
    if resume_mode and _step3_done(job_dir):
        _progress("step3", "skip")
    else:
        _progress("step3", "running")
        cores = str(get_optimal_cores(is_aws=is_aws))
        run_subprocess(
            ["lightdock3.py", "setup.json", "100", "-c", cores],
            cwd=job_dir,
            log_path=job_dir / "step3_simulation.log",
            on_line=log_cb,
        )
        _end3 = now_kst()
        save_timing(job_dir, "step3", _start3, _end3)
        _progress("step3", "done",
                  format_duration_min(f"{(_end3 - _start3).total_seconds()}s"))

    # ── STEP 4: 포즈 추출 + BSAS ─────────────────────────────────────────────
    _start4 = now_kst()
    if resume_mode and _step4_done(job_dir):
        _progress("step4", "skip")
    else:
        _progress("step4", "running")

        # setup.json에서 PDB 파일명 복원 (폴백 포함)
        try:
            target_rec, target_lig = _parse_setup_json(job_dir)
        except Exception as e:
            _warn(f"setup.json 파싱 예외({e}). 파일명 자동 탐색 중...")
            pdb_files = [f for f in list(job_dir.glob("*.pdb")) + list(job_dir.glob("*.pdb1"))
                         if not f.name.startswith("lightdock_")]
            if len(pdb_files) >= 2:
                r_files = [f for f in pdb_files if "rec" in f.name.lower()]
                l_files = [f for f in pdb_files if "lig" in f.name.lower()]
                if r_files and l_files:
                    target_rec, target_lig = r_files[0].name, l_files[0].name
                else:
                    target_rec, target_lig = pdb_files[0].name, pdb_files[1].name
            else:
                target_rec, target_lig = "receptor.pdb", "ligand.pdb"

        gen_swarm_dirs = iter_swarm_dirs(job_dir)
        total_sw = len(gen_swarm_dirs)
        for gi, sd in enumerate(gen_swarm_dirs):
            gso_files = sorted(
                sd.glob("gso_*.out"),
                key=lambda f: int(re.search(r"gso_(\d+)", f.name).group(1))
                if re.search(r"gso_(\d+)", f.name) else 0,
            )
            if not gso_files:
                continue
            run_subprocess(
                ["lgd_generate_conformations.py",
                 f"../{target_rec}", f"../{target_lig}",
                 gso_files[-1].name, "200"],
                cwd=sd,
                log_path=job_dir / f"step4_generate_{sd.name}.log",
                on_line=log_cb,
            )
            if progress_cb:
                progress_cb("step4_progress", f"{gi+1}/{total_sw}", None)

        _end4 = now_kst()
        save_timing(job_dir, "step4", _start4, _end4)
        _progress("step4", "done",
                  format_duration_min(f"{(_end4 - _start4).total_seconds()}s"))

    # BSAS 클러스터링 (skip 조건 독립 판정)
    if not (resume_mode and _clustering_done(job_dir)):
        run_bsas_clustering(job_dir)

    # ── STEP 5: 랭킹 ────────────────────────────────────────────────────────
    _start5 = now_kst()
    if resume_mode and (job_dir / "rank_by_scoring.list").exists():
        _progress("step5", "skip")
    else:
        _progress("step5", "running")
        sw_count = len(list(job_dir.glob("swarm_*")))
        gso_r5 = sorted(
            (job_dir / "swarm_0").glob("gso_*.out"),
            key=lambda f: int(re.search(r"gso_(\d+)", f.name).group(1))
            if re.search(r"gso_(\d+)", f.name) else 0,
        )
        sim_steps = int(re.search(r"gso_(\d+)", gso_r5[-1].name).group(1)) if gso_r5 else 100
        rc5 = run_subprocess(
            ["lgd_rank.py", str(sw_count), str(sim_steps)],
            cwd=job_dir,
            log_path=job_dir / "step5_rank.log",
            on_line=log_cb,
        )
        _end5 = now_kst()
        save_timing(job_dir, "step5", _start5, _end5)
        dur5 = format_duration_min(f"{(_end5 - _start5).total_seconds()}s")
        _progress("step5", "done" if rc5 == 0 else "error", dur5)
        if rc5 != 0:
            result["errors"].append("lgd_rank.py 실행 실패. step5_rank.log를 확인하세요.")
            return result

    # ── STEP 5.5: fnat 필터링 ────────────────────────────────────────────────
    restraints_file = job_dir / "restraints.list"
    rank_scoring    = job_dir / "rank_by_scoring.list"
    if restraints_file.exists() and rank_scoring.exists():
        filtered_dir = job_dir / "filtered"
        if filtered_dir.exists():
            shutil.rmtree(filtered_dir)
        try:
            # sys.executable로 직접 호출해 실행 권한 문제를 우회
            _filter_script = Path(sys.executable).parent / "lgd_filter_restraints.py"
            filter_rc = run_subprocess(
                [sys.executable, str(_filter_script),
                 "rank_by_scoring.list", "restraints.list", "A", "B", "--fnat", "0.0"],
                cwd=job_dir,
                log_path=job_dir / "step5_filter_restraints.log",
                on_line=log_cb,
            )
            if filter_rc != 0:
                _warn(
                    "lgd_filter_restraints.py 실행 실패 — fnat 결과 미생성. "
                    "체인 ID(A/B)가 다른 경우 수동 재실행이 필요합니다."
                )
        except Exception as fe:
            _warn(f"restraints 필터링 단계 예외: {fe}")

    result["success"] = True
    return result
