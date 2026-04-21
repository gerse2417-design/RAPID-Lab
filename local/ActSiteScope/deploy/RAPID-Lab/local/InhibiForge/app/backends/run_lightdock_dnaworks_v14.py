"""
run_lightdock_dnaworks_v14.py
──────────────────────────────────────────────────────────────
LightDock + DNAWorks 통합 도킹 파이프라인 — 단일 파일 버전 (v14).

v13 베이스 + 다음 변경:
  - `save_timing()` 이 KST(UTC+9) 타임존을 명시하여 step_timing.json 에 기록.
    (포맷: 'YYYY-MM-DD HH:MM:SS+0900')
  - 시각 기록용 헬퍼 `_now_kst()` 추가. datetime.datetime.now() 호출을 대체.

v13 원본은 건드리지 않고 이 파일을 새로 추가.

amr_dock_backend/ 패키지의 모든 모듈을 하나의 파일로 통합한 버전.
app_amr_dock_interactive_v13.py 가 필요로 하는 모든 함수/상수를 포함한다.

v13 변경 사항 (vs v12):
  - 신규 파서: parse_rank_filtered(), parse_dnaworks_logfile()
  - 신규 스크립트 생성기: build_pymol_heatmap_script(), build_dnaworks_input()
  - 패키지(amr_dock_backend/) 의 최신 pipeline.py 반영

포함 섹션:
  1. env          — 외부 도구 설치 감지, 경로 상수, CPU 코어 계산
  2. subprocess   — run_subprocess, save_timing, format_duration_min
  3. pdb_utils    — PDB 전처리(비표준 치환, PDBFixer, HETATM 제거), 서열 추출
  4. parsers      — rank_list / cluster.repr / rank_filtered / dnaworks logfile 파서
  5. swarm_utils  — swarm_N 순회·cluster 수집·PDB 목록
  6. restraints   — 핫스팟 원자 번호 파싱, restraints.list 생성
  7. scripts      — PyMOL Heatmap / ChimeraX / DNAWorks 입력 스크립트 빌더
  8. analysis     — 인터페이스 잔기 탐지, 에너지 지형, swarm 좌표
  9. job          — get_job_summary (단계 완료 상태 분석)
 10. pipeline     — run_full_docking_pipeline, run_bsas_clustering
 11. CLI          — argparse 기반 __main__ 엔트리포인트

사용 예시 (CLI):
    python run_lightdock_dnaworks_v14.py \\
        --job-dir ./Projects/2UUY \\
        --inputs-dir ./Projects/2UUY \\
        --rec-name 2UUY_rec.pdb \\
        --lig-name 2UUY_lig.pdb

사용 예시 (import):
    from run_lightdock_dnaworks_v14 import (
        run_full_docking_pipeline, parse_rank_list, get_tool_status
    )
"""

# ── 통합 import 섹션 ────────────────────────────────────────────────────────
import datetime
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

# ── KST 타임존 헬퍼 (v14 추가) ─────────────────────────────────────────────
KST = datetime.timezone(datetime.timedelta(hours=9))


def _now_kst() -> datetime.datetime:
    """KST(UTC+9) 기준의 timezone-aware 현재 시각."""
    return datetime.datetime.now(KST)


def _to_kst(dt: datetime.datetime) -> datetime.datetime:
    """naive datetime 은 KST 로 간주해 tz 를 부여하고, aware 는 KST 로 변환."""
    if dt.tzinfo is None:
        return dt.replace(tzinfo=KST)
    return dt.astimezone(KST)

import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# PDBFixer 선택적 의존성 가드 (미설치 시에도 모듈 정상 로드)
try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    PDBFIXER_AVAILABLE = True
except ImportError:
    PDBFIXER_AVAILABLE = False


# ============================================================================
# 1) ENV — 외부 도구 설치 감지 및 실행 환경 상수
# ============================================================================

# ── 경로 상수 (환경변수 override 지원) ──────────────────────────────────────

DNAWORKS_EXECUTABLE = Path(
    os.environ.get("DNAWORKS_BIN", "/home/sooyeon/amr/DNAWorks/dnaworks")
)

CHIMERAX_BIN = os.environ.get("CHIMERAX_BIN", "/usr/bin/chimerax-daily")

_PYMOL_FIXED = Path(os.environ.get("PYMOL_BIN", "/usr/bin/pymol"))


def get_optimal_cores(is_aws: bool = False) -> int:
    """
    운영 환경에 맞춰 안정적으로 사용할 CPU 코어 수를 계산한다.

    Parameters
    ----------
    is_aws : True이면 AWS 배포 환경으로 간주 (코어-1 반환).
             False이면 로컬 환경으로 간주 (전체의 약 70% 반환).
    """
    total_cores = os.cpu_count()
    if not total_cores:
        return 1
    if is_aws:
        return max(1, total_cores - 1)
    return max(1, math.floor(total_cores * 0.7))


def _check_lightdock() -> bool:
    try:
        result = subprocess.run(
            ["lightdock3_setup.py", "-h"],
            capture_output=True, timeout=5,
        )
        return result.returncode == 0
    except Exception:
        return False


def _resolve_pymol_bin() -> str | None:
    """PyMOL 실행 경로를 반환한다. 고정 경로 우선, 없으면 PATH 탐색."""
    if _PYMOL_FIXED.exists():
        return str(_PYMOL_FIXED)
    return shutil.which("pymol")


def _check_pymol(pymol_bin: str | None) -> bool:
    if pymol_bin is None:
        return False
    try:
        _clean_env = {k: v for k, v in os.environ.items()
                      if k not in ("PYTHONHOME", "PYTHONPATH", "VIRTUAL_ENV")}
        result = subprocess.run(
            [pymol_bin, "--version"],
            capture_output=True, text=True, timeout=10,
            env=_clean_env,
        )
        return result.returncode == 0
    except Exception:
        return False


def get_tool_status() -> dict:
    """
    현재 환경에서 외부 도구의 설치 여부를 확인하여 반환한다.

    Returns
    -------
    dict: {
      "lightdock":  bool,
      "dnaworks":   bool,
      "pymol":      bool,
      "pymol_bin":  str | None,
    }
    """
    pymol_bin = _resolve_pymol_bin()
    return {
        "lightdock": _check_lightdock(),
        "dnaworks":  DNAWORKS_EXECUTABLE.exists(),
        "pymol":     _check_pymol(pymol_bin),
        "pymol_bin": pymol_bin,
    }


# ============================================================================
# 2) SUBPROCESS — 외부 명령어 실행 래퍼 및 타이밍 유틸
# ============================================================================

def run_subprocess(
    cmd: list,
    cwd: Path,
    log_path: Path | None = None,
    on_line=None,
) -> int:
    """
    외부 명령어를 실행하고 실시간으로 출력 라인을 콜백으로 전달한다.

    Parameters
    ----------
    cmd      : 실행할 명령어 리스트
    cwd      : 작업 디렉토리
    log_path : 로그 파일 저장 경로 (None이면 저장 안 함)
    on_line  : 실시간 출력 라인 콜백 (line: str) -> None

    Returns
    -------
    int : 프로세스 종료 코드
    """
    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    output_lines = []

    try:
        process = subprocess.Popen(
            cmd, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        for line in process.stdout:
            stripped = line.rstrip()
            output_lines.append(stripped)
            if on_line:
                on_line(stripped)
        returncode = process.wait()
    except Exception as e:
        output_lines.append(f"[ERROR] {e}")
        returncode = -1

    if log_path:
        try:
            Path(log_path).write_text("\n".join(output_lines), encoding="utf-8")
        except Exception:
            pass

    return returncode


def save_timing(
    job_dir,
    step_id: str,
    start_time: datetime.datetime,
    end_time: datetime.datetime,
) -> None:
    """단계별 실행 시간을 step_timing.json에 누적 저장한다."""
    timing_file = Path(job_dir) / "step_timing.json"
    data = {}
    if timing_file.exists():
        try:
            with open(timing_file, "r") as f:
                data = json.load(f)
        except Exception:
            pass

    duration = (end_time - start_time).total_seconds()
    start_kst = _to_kst(start_time)
    end_kst   = _to_kst(end_time)
    data[step_id] = {
        "start":    start_kst.strftime("%Y-%m-%d %H:%M:%S%z"),
        "end":      end_kst.strftime("%Y-%m-%d %H:%M:%S%z"),
        "duration": f"{duration:.2f}s",
    }
    with open(timing_file, "w") as f:
        json.dump(data, f)


def format_duration_min(dur_str: str) -> str:
    """
    '1442.02s' 형식의 문자열을 '24.03min' 형식으로 변환한다.
    파싱 실패 시 원본 문자열을 그대로 반환한다.
    """
    try:
        sec = float(str(dur_str).rstrip("s"))
        return f"{sec / 60:.2f}min"
    except Exception:
        return dur_str


# ============================================================================
# 3) PDB_UTILS — PDB 파일 전처리 유틸리티
# ============================================================================

# 비표준 아미노산 → 표준 아미노산 매핑 테이블
NONSTANDARD_RESIDUE_MAP = {
    "MSE": ("MET", "ATOM"),   # Selenomethionine
    "HYP": ("PRO", "ATOM"),   # Hydroxyproline
    "CSE": ("CYS", "ATOM"),   # Selenocysteine
    "TPO": ("THR", "ATOM"),   # Phosphothreonine
    "SEP": ("SER", "ATOM"),   # Phosphoserine
    "PTR": ("TYR", "ATOM"),   # Phosphotyrosine
    "MLY": ("LYS", "ATOM"),   # N-methyl lysine
    "CME": ("CYS", "ATOM"),   # S,S-(2-hydroxyethyl)thiocysteine
    "OCS": ("CYS", "ATOM"),   # Cysteinesulfonic acid
    "KCX": ("LYS", "ATOM"),   # Lysine NZ-carboxylic acid
}


def get_sequence_from_pdb(pdb_path) -> str:
    """PDB 파일에서 아미노산 서열(1-letter code)을 추출."""
    if not pdb_path or not Path(pdb_path).exists():
        return ""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        full_seq = []
        for model in structure:
            for chain in model:
                chain_seq = ""
                for residue in chain:
                    if residue.get_id()[0] == " ":  # standard amino acids only
                        try:
                            aa = seq1(residue.get_resname())
                            if aa and aa != "X" and aa != "*":
                                chain_seq += aa
                        except Exception:
                            continue
                if chain_seq:
                    full_seq.append(chain_seq)
        return "".join(full_seq)
    except Exception as e:
        return f"Error: {str(e)}"


def substitute_nonstandard_residues(input_path, output_path, residue_map=None) -> dict:
    """
    PDB 파일의 HETATM 레코드 중 알려진 비표준 아미노산을
    표준 아미노산 ATOM 레코드로 치환한다.

    Returns
    -------
    dict: {"substitutions": {원래잔기명: 치환횟수}}
    """
    if residue_map is None:
        residue_map = NONSTANDARD_RESIDUE_MAP

    substitutions = {}
    input_path = Path(input_path)
    if not input_path.exists():
        return {"substitutions": substitutions}

    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            if line[:6] == "HETATM":
                resname = line[17:20].strip()
                if resname in residue_map:
                    canonical, new_rec = residue_map[resname]
                    line = f"{new_rec:<6}" + line[6:17] + f"{canonical:<3}" + line[20:]
                    substitutions[resname] = substitutions.get(resname, 0) + 1
            fout.write(line)

    return {"substitutions": substitutions}


def clean_pdb_for_lightdock(input_path, output_path) -> bool:
    """
    LightDock ANM 계산 시 원자 수 불일치 오류 방지:
    표준 아미노산(ATOM)만 남기고 HETATM(물, 이온 등)을 제거.
    """
    input_path = Path(input_path)
    if not input_path.exists():
        return False
    try:
        with open(input_path, "r") as fin, open(output_path, "w") as fout:
            for line in fin:
                if line.startswith("ATOM  "):
                    fout.write(line)
                elif line.startswith("TER") or line.startswith("END"):
                    fout.write(line)
        return True
    except Exception as e:
        raise RuntimeError(f"PDB 클리닝 중 오류: {e}") from e


def run_pdbfixer_pipeline(input_path, output_path,
                          add_missing_atoms=True, remove_heterogens=True) -> dict:
    """
    OpenMM PDBFixer로 PDB 파일의 누락 원자(백본/측쇄)를 보완하고
    잔여 HETATM을 제거한다.

    Returns
    -------
    dict: {"success": bool, "missing_atoms": int, "skipped_reason": str or None}
    """
    if not PDBFIXER_AVAILABLE:
        return {
            "success": False,
            "missing_atoms": 0,
            "skipped_reason": "pdbfixer 미설치 — 이 단계를 건너뜁니다.",
        }
    try:
        fixer = PDBFixer(filename=str(input_path))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        n_atoms = sum(len(v) for v in fixer.missingAtoms.values())
        if add_missing_atoms:
            fixer.addMissingAtoms()
        if remove_heterogens:
            fixer.removeHeterogens(keepWater=False)
        with open(output_path, "w") as fout:
            PDBFile.writeFile(fixer.topology, fixer.positions, fout)
        return {"success": True, "missing_atoms": n_atoms, "skipped_reason": None}
    except Exception as exc:
        return {
            "success": False,
            "missing_atoms": 0,
            "skipped_reason": f"pdbfixer 오류: {exc}",
        }


def preprocess_pdb_for_docking(input_path, output_path,
                                substitute_nonstandard=False,
                                run_pdbfixer=False,
                                strict_hetatm_removal=True,
                                logger=None) -> dict:
    """
    LightDock 도킹 전 PDB 전처리 파이프라인 오케스트레이터.

    실행 순서:
      1. 비표준 아미노산 치환  (substitute_nonstandard=True 시)
      2. PDBFixer 누락 원자 보완 (run_pdbfixer=True 시)
      3. 엄격 HETATM 제거       (strict_hetatm_removal=True 시, 기본값)

    Returns
    -------
    dict: {"substitutions": dict, "pdbfixer_result": dict or None, "warnings": list}
    """
    summary = {"substitutions": {}, "pdbfixer_result": None, "warnings": []}

    with tempfile.TemporaryDirectory() as tmpdir:
        current = str(input_path)

        if substitute_nonstandard:
            step1_out = f"{tmpdir}/step1_subst.pdb"
            r = substitute_nonstandard_residues(current, step1_out)
            summary["substitutions"] = r["substitutions"]
            if logger and r["substitutions"]:
                logger("치환: " + ", ".join(
                    f"{k}→{NONSTANDARD_RESIDUE_MAP[k][0]}({v}건)"
                    for k, v in r["substitutions"].items()
                ))
            current = step1_out

        if run_pdbfixer:
            step2_out = f"{tmpdir}/step2_fixed.pdb"
            pf = run_pdbfixer_pipeline(current, step2_out)
            summary["pdbfixer_result"] = pf
            if pf["success"]:
                current = step2_out
                if logger:
                    logger(f"PDBFixer: {pf['missing_atoms']}개 원자 보완 완료")
            else:
                summary["warnings"].append(pf["skipped_reason"])
                if logger:
                    logger(f"WARNING: {pf['skipped_reason']}")

        if strict_hetatm_removal:
            step3_out = f"{tmpdir}/step3_clean.pdb"
            clean_pdb_for_lightdock(current, step3_out)
            current = step3_out

        shutil.copy2(current, output_path)

    return summary


# ============================================================================
# 4) PARSERS — LightDock / DNAWorks 결과 파일 파싱
# ============================================================================

def parse_rank_list(rank_file: Path) -> pd.DataFrame:
    """
    LightDock 랭킹 파일(rank_by_scoring.list 등)을 DataFrame으로 변환한다.

    컬럼: Swarm, Glowworm, Coordinates, RecID, LigID, Luciferin,
          Neigh, VR, RMSD, PDB, Clashes, Scoring
    """
    rank_file = Path(rank_file)
    if not rank_file.exists():
        return pd.DataFrame()

    COLUMNS = [
        "Swarm", "Glowworm", "Coordinates", "RecID", "LigID",
        "Luciferin", "Neigh", "VR", "RMSD", "PDB", "Clashes", "Scoring"
    ]

    rows = []
    try:
        with open(rank_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("Swarm"):
                    continue
                collapsed = re.sub(r'\(.*?\)', 'COORDS', line)
                tokens = collapsed.split()
                if len(tokens) < len(COLUMNS):
                    continue
                row = tokens[:len(COLUMNS)]
                rows.append(row)
    except Exception:
        return pd.DataFrame()

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows, columns=COLUMNS)
    for col in ["Swarm", "Glowworm", "RecID", "LigID", "Neigh", "Clashes"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["Luciferin", "VR", "RMSD", "Scoring"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def parse_cluster_repr(swarm_dir: Path) -> pd.DataFrame:
    """
    swarm_N/cluster.repr 파일을 파싱하여 클러스터 정보 DataFrame을 반환한다.

    반환 컬럼: Cluster ID, Population, Best Score, Best Glowworm, Representative PDB
    """
    cluster_file = Path(swarm_dir) / "cluster.repr"
    if not cluster_file.exists():
        return pd.DataFrame()

    rows = []
    try:
        with open(cluster_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                tokens = line.split(":")
                if len(tokens) < 5:
                    continue
                try:
                    rows.append({
                        "Cluster ID":        int(tokens[0]),
                        "Population":        int(tokens[1]),
                        "Best Score":        float(tokens[2]),
                        "Best Glowworm":     int(tokens[3]),
                        "Representative PDB": tokens[4],
                    })
                except (ValueError, IndexError):
                    continue
    except Exception:
        return pd.DataFrame()

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def parse_rank_filtered(rank_filtered_path: Path) -> pd.DataFrame:
    """
    filtered/rank_filtered.list 파일을 파싱하여 fnat 결과 DataFrame을 반환한다.

    반환 컬럼: Swarm, Glowworm, Scoring, fnat
    """
    rank_filtered_path = Path(rank_filtered_path)
    if not rank_filtered_path.exists():
        return pd.DataFrame()

    pat = re.compile(r"swarm_(\d+)_(\d+)\.pdb")
    rows = []
    try:
        for raw in rank_filtered_path.read_text(encoding="utf-8", errors="replace").splitlines():
            raw = raw.strip()
            if not raw:
                continue
            parts = raw.split()
            if len(parts) < 3:
                continue
            m = pat.search(parts[0])
            if not m:
                continue
            try:
                rows.append({
                    "Swarm":    int(m.group(1)),
                    "Glowworm": int(m.group(2)),
                    "Scoring":  float(parts[1]),
                    "fnat":     float(parts[2]),
                })
            except ValueError:
                continue
    except Exception:
        return pd.DataFrame()

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def parse_dnaworks_logfile(log_path: Path) -> tuple:
    """
    DNAWorks LOGFILE.txt를 파싱하여 Trial별 올리고 목록과 최종 요약을 반환한다.

    Returns
    -------
    (trial_oligos, summary_rows)
      trial_oligos : dict  {trial_num: [{"ID", "Sequence (DNA)", "Length"}, ...]}
      summary_rows : list  [{"Trial", "Tm", "Len", "Score", ...}, ...]
                 Score 오름차순 정렬, 동일 Score 시 TmRange 오름차순 2차 정렬 적용.
    파싱 실패 또는 파일 없음 시 ({}, []) 반환.
    """
    log_path = Path(log_path)
    if not log_path.exists():
        return {}, []

    trial_oligos: dict = {}
    summary_rows: list = []

    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
        current_trial = None
        is_parsing_oligo = False
        is_parsing_summary = False

        for line in lines:
            trial_match = re.search(r'PARAMETERS FOR TRIAL\s+(\d+)', line)
            if trial_match:
                current_trial = int(trial_match.group(1))
                is_parsing_oligo = False
                continue

            if "oligonucleotides need to be synthesized" in line:
                is_parsing_oligo = True
                if current_trial not in trial_oligos:
                    trial_oligos[current_trial] = []
                continue

            if "FINAL SUMMARY FOR" in line:
                is_parsing_oligo = False
                is_parsing_summary = True
                continue

            if is_parsing_oligo:
                stripped = line.strip()
                if "FINAL SUMMARY" in line or stripped == "":
                    is_parsing_oligo = False
                    continue
                if stripped.startswith("-"):
                    continue
                parts = stripped.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    trial_oligos.setdefault(current_trial, []).append({
                        "ID": int(parts[0]),
                        "Sequence (DNA)": parts[1],
                        "Length": int(parts[-1])
                    })
                continue

            if is_parsing_summary:
                if line.strip() == "" or line.startswith("-") or line.startswith("|") or "#" in line:
                    continue
                parts = line.strip().split()
                if len(parts) >= 11 and parts[0].isdigit():
                    try:
                        pipe_idx = parts.index("|")
                        summary_rows.append({
                            "Trial":     int(parts[0]),
                            "Tm":        float(parts[1]),
                            "Len":       int(parts[2]),
                            "Score":     float(parts[pipe_idx + 1]),
                            "TmRange":   float(parts[pipe_idx + 2]),
                            "Short":     int(parts[pipe_idx + 3]),
                            "Long":      int(parts[pipe_idx + 4]),
                            "#Oligos":   int(parts[pipe_idx + 5]),
                            "#Repeat":   int(parts[pipe_idx + 6]),
                            "#Misprime": int(parts[pipe_idx + 7]),
                        })
                    except (ValueError, IndexError):
                        continue

    except Exception:
        return {}, []

    # Score 오름차순 정렬, Score가 동일(특히 0.000)할 경우 TmRange 오름차순 2차 정렬
    summary_rows.sort(key=lambda r: (r["Score"], r["TmRange"]))
    return trial_oligos, summary_rows


# ============================================================================
# 5) SWARM_UTILS — 스왐 디렉토리 순회 및 PDB 포즈 수집
# ============================================================================

def iter_swarm_dirs(job_dir: Path) -> list:
    """job_dir 내 swarm_N 디렉토리를 N 오름차순으로 정렬하여 반환한다."""
    job_dir = Path(job_dir)
    return sorted(
        (d for d in job_dir.glob("swarm_*") if d.is_dir()),
        key=lambda d: int(re.search(r'swarm_(\d+)', d.name).group(1))
        if re.search(r'swarm_(\d+)', d.name) else 0
    )


def collect_all_swarm_clusters(job_dir: Path) -> dict:
    """
    job_dir 내 모든 swarm_N 디렉토리를 순회하여 cluster.repr 데이터를 수집한다.

    Returns
    -------
    dict: {"swarm_0": DataFrame, "swarm_1": DataFrame, ...}
    """
    result = {}
    for swarm_dir in iter_swarm_dirs(job_dir):
        df = parse_cluster_repr(swarm_dir)
        if not df.empty:
            result[swarm_dir.name] = df
    return result


def get_swarm_pdb_files(job_dir: Path) -> dict:
    """
    job_dir 내 swarm_N 디렉토리별로 lightdock_*.pdb 파일 목록을 반환한다.

    Returns
    -------
    dict: {"swarm_0": [Path(...), ...], "swarm_1": [...], ...}
    """
    result = {}
    for swarm_dir in iter_swarm_dirs(job_dir):
        pdbs = sorted(
            swarm_dir.glob("lightdock_*.pdb"),
            key=lambda f: int(re.search(r'lightdock_(\d+)', f.name).group(1))
            if re.search(r'lightdock_(\d+)', f.name) else 0
        )
        if pdbs:
            result[swarm_dir.name] = pdbs
    return result


# ============================================================================
# 6) RESTRAINTS — 핫스팟 원자 제한 도킹 유틸리티
# ============================================================================

def parse_atom_input(raw_input: str) -> list:
    """
    사용자가 입력한 원자 번호 문자열을 파싱하여 (체인ID, Atom Serial Number) 쌍 목록을 반환한다.

    입력 예시: "A425, A426, A427, B12"
    출력: [('A', 425), ('A', 426), ('A', 427), ('B', 12)]
    """
    pattern = re.compile(r'([A-Za-z])(\d+)')
    matches = pattern.findall(raw_input.upper())
    if not matches:
        raise ValueError(
            f"원자 번호 형식을 인식하지 못했습니다: '{raw_input}'\n"
            "올바른 형식 예시: A425, A426, A427"
        )
    return [(chain, int(serial)) for chain, serial in matches]


def parse_pdb_atoms(pdb_path: Path) -> list:
    """
    PDB 파일의 ATOM 레코드를 고정 폭(fixed-width) 슬라이싱으로 파싱한다.

    Returns
    -------
    list of dict: {"serial": int, "name": str, "resname": str, "chain": str, "resseq": int}
    """
    pdb_path = Path(pdb_path)
    atoms = []
    if not pdb_path.exists():
        return atoms

    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if not line.startswith("ATOM  "):
                    continue
                try:
                    serial  = int(line[6:11].strip())
                    name    = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain   = line[21].strip()
                    resseq  = int(line[22:26].strip())
                    atoms.append({
                        "serial":  serial,
                        "name":    name,
                        "resname": resname,
                        "chain":   chain,
                        "resseq":  resseq,
                    })
                except (ValueError, IndexError):
                    continue
    except Exception:
        pass

    return atoms


def generate_restraints_list(
    receptor_pdb: Path,
    raw_input: str,
    output_path: Path,
    role: str = "R",
) -> tuple:
    """
    사용자 입력 원자 번호를 PDB 파일과 대조하여 LightDock restraints.list 파일을 생성한다.

    Returns
    -------
    (lines_written: list[str], not_found: list[tuple])
    """
    atom_inputs = parse_atom_input(raw_input)
    all_atoms   = parse_pdb_atoms(receptor_pdb)

    atom_map = {(a["chain"], a["serial"]): (a["resname"], a["resseq"]) for a in all_atoms}

    residue_map = {}
    not_found = []

    for chain, serial in atom_inputs:
        key = (chain, serial)
        if key in atom_map:
            resname, resseq = atom_map[key]
            res_key = (chain, resseq)
            if res_key not in residue_map:
                residue_map[res_key] = resname
        else:
            not_found.append(key)

    lines_written = []
    for (chain, resseq), resname in residue_map.items():
        lines_written.append(f"{role} {chain}.{resname}.{resseq}")

    with open(output_path, "w") as f:
        f.write("\n".join(lines_written) + ("\n" if lines_written else ""))

    return lines_written, not_found


# ============================================================================
# 7) SCRIPTS — 외부 도구용 스크립트 문자열 생성기
# ============================================================================

def build_pymol_heatmap_script(receptor_pdb: str, swarm_df: "pd.DataFrame", output_png: str) -> str:
    """
    PyMOL headless .pml 스크립트 문자열을 반환한다.

    - receptor.pdb: surface 모드, 투명도 20%
    - 각 swarm 좌표에 pseudoatom (b-factor = Scoring)
    - spectrum b, blue_red 로 컬러링
    - ray 렌더 후 PNG 저장
    """
    lines = [
        f"load {receptor_pdb}, receptor",
        "show surface, receptor",
        "set transparency, 0.20, receptor",
        "color grey80, receptor",
        "",
        "# Swarm center pseudoatoms",
    ]
    for _, row in swarm_df.iterrows():
        name = f"swarm_{int(row['Swarm'])}"
        lines.append(
            f"pseudoatom {name}, pos=[{row['x']:.3f},{row['y']:.3f},{row['z']:.3f}]"
        )
        lines.append(f"alter {name}, b={row['Scoring']:.5f}")

    lines += [
        "",
        "select all_swarms, pseudoatom",
        "show spheres, all_swarms",
        "set sphere_scale, 1.5, all_swarms",
        "spectrum b, blue_red, all_swarms",
        "",
        "bg_color white",
        "orient",
        "ray 1920, 1080",
        f"png {output_png}, dpi=150",
        "quit",
    ]
    return "\n".join(lines)


def build_chimerax_script(
    complex_pdb: str,
    interface_res: dict,
    mpbind_residues: list,
    dx_file: str = None,
    output_dir: str = ".",
) -> str:
    """
    ChimeraX .cxc 스크립트 문자열을 반환한다.

    - 복합체 cartoon 로드
    - 인터페이스 잔기 stick 모드 (거리 < 4.0Å)
    - hbonds 표시
    - MPBind 잔기 gold highlight
    - dx_file 있으면 coulombic surface
    - 360도 회전 MP4 저장 + 2K PNG 저장
    """
    lines = [
        f"open {complex_pdb}",
        "preset 'overall look' 'publication 2 (depth-cued)'",
        "",
    ]

    if dx_file:
        lines += [
            f"open {dx_file}",
            "volume #2 style surface",
            "coulombic #1 surfaces #2",
            "",
        ]

    chains = interface_res.get("chains", [])
    pairs  = interface_res.get("pairs", [])
    if pairs and len(chains) >= 2:
        ca, cb = chains[0], chains[1]
        res_a_ids = ",".join(str(p[1]) for p in pairs)
        res_b_ids = ",".join(str(p[4]) for p in pairs)
        lines += [
            "# 인터페이스 잔기 (거리 < 4.0Å)",
            f"style /{ca}:{res_a_ids} stick",
            f"style /{cb}:{res_b_ids} stick",
            "hbonds",
            "",
        ]

    if mpbind_residues:
        for chain, resid in mpbind_residues:
            lines.append(f"color /{chain}:{resid} gold")
        lines.append("")

    png_path = f"{output_dir}/chimerax_interface.png"
    mp4_path = f"{output_dir}/chimerax_360.mp4"
    lines += [
        "view",
        f"save {png_path} width 2048 height 2048 supersample 3",
        "",
        "# 360도 회전 MP4",
        "movie record",
        "turn y 2 180",
        f"movie encode {mp4_path} quality high",
        "",
        "exit",
    ]
    return "\n".join(lines)


def build_dnaworks_input(params: dict) -> str:
    """
    파라미터 dict에서 DNAWORKS.inp 전체 텍스트를 생성하여 반환한다.

    Parameters
    ----------
    params : dict 키 목록
      - title, logfile, melting_temp, oligo_length, codon_org,
        sodium_conc, magnesium_conc, repeat_limit, n_solutions,
        misprime_check, protein_seq
    """
    title          = params.get("title", "DNAWorks Job")
    logfile        = params.get("logfile", "LOGFILE.txt")
    melting_temp   = params.get("melting_temp", 60)
    oligo_length   = params.get("oligo_length", 50)
    codon_org      = params.get("codon_org", "ecoli2")
    sodium_conc    = params.get("sodium_conc", 0.05)
    magnesium_conc = params.get("magnesium_conc", 0.002)
    repeat_limit   = params.get("repeat_limit", 8)
    n_solutions    = params.get("n_solutions", 5)
    misprime       = params.get("misprime_check", True)
    protein_seq    = params.get("protein_seq", "")

    lines = [
        f'title "{title}"',
        f'logfile "{logfile}"',
        f'melting low {melting_temp}',
        f'length low {oligo_length}',
        f'codon {codon_org}',
        f'concentration sodium {sodium_conc} magnesium {magnesium_conc}',
        f'repeat {repeat_limit}',
        f'solutions {n_solutions}',
    ]
    if misprime:
        lines.append('misprime 18 tip 6 max 8')
    lines += [
        'protein',
        f'  {protein_seq.strip()}',
        '//',
    ]
    return "\n".join(lines) + "\n"


# ============================================================================
# 8) ANALYSIS — 도킹 결과 과학 분석
# ============================================================================

def find_interface_residues(complex_pdb: Path, cutoff: float = 4.0) -> dict:
    """
    BioPython으로 복합체 PDB의 체인 간 거리 < cutoff Å인 잔기 쌍을 탐지한다.

    Returns
    -------
    dict: {
      "chains": [chain_id, ...],
      "pairs":  [(chain_A, resid_A, resname_A, chain_B, resid_B, resname_B), ...]
    }
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", str(complex_pdb))
        model = next(structure.get_models())
        chains = list(model.get_chains())
        if len(chains) < 2:
            return {"chains": [c.id for c in chains], "pairs": []}

        chain_a, chain_b = chains[0], chains[1]
        res_a = [r for r in chain_a if r.get_id()[0] == " "]

        pairs = []
        seen_a = set()
        atoms_b = list(chain_b.get_atoms())
        for ra in res_a:
            atoms_a = list(ra.get_atoms())
            if not atoms_a:
                continue
            min_dist = float("inf")
            closest_atom_b = None
            for atom_a in atoms_a:
                for atom_b in atoms_b:
                    try:
                        d = (atom_a.get_vector() - atom_b.get_vector()).norm()
                        if d < min_dist:
                            min_dist = d
                            closest_atom_b = atom_b
                    except Exception:
                        continue
            if min_dist < cutoff and closest_atom_b is not None:
                key_a = (chain_a.id, ra.get_id()[1])
                if key_a not in seen_a:
                    seen_a.add(key_a)
                    rb_res = closest_atom_b.get_parent()
                    pairs.append((chain_a.id, ra.get_id()[1], ra.get_resname(),
                                  chain_b.id, rb_res.get_id()[1], rb_res.get_resname()))
        return {"chains": [chain_a.id, chain_b.id], "pairs": pairs}
    except Exception as e:
        return {"chains": [], "pairs": [], "error": str(e)}


def build_energy_landscape_df(job_dir: Path) -> pd.DataFrame:
    """
    rank_by_scoring.list + 전체 swarm cluster.repr 병합.

    컬럼: Swarm, Glowworm, RMSD, Scoring, PDB, ClusterID, Population
    """
    rank_file = Path(job_dir) / "rank_by_scoring.list"
    df_rank = parse_rank_list(rank_file)
    if df_rank.empty:
        return pd.DataFrame()

    cluster_map = {}
    all_swarm_clusters = collect_all_swarm_clusters(job_dir)
    for swarm_name, df_cl in all_swarm_clusters.items():
        swarm_m = re.search(r'swarm_(\d+)', swarm_name)
        if not swarm_m:
            continue
        swarm_idx = int(swarm_m.group(1))
        for _, crow in df_cl.iterrows():
            cluster_map[(swarm_idx, int(crow["Best Glowworm"]))] = (
                int(crow["Cluster ID"]),
                int(crow["Population"]),
            )

    df_rank["ClusterID"] = df_rank.apply(
        lambda r: cluster_map.get((int(r["Swarm"]), int(r["Glowworm"])), (None, None))[0],
        axis=1,
    )
    df_rank["Population"] = df_rank.apply(
        lambda r: cluster_map.get((int(r["Swarm"]), int(r["Glowworm"])), (None, None))[1],
        axis=1,
    )

    return df_rank


def parse_swarm_coordinates(rank_file: Path) -> pd.DataFrame:
    """
    rank_by_scoring.list에서 각 Swarm의 최고 Scoring 포즈의 (x, y, z, score)를 추출한다.

    Returns
    -------
    pd.DataFrame — 컬럼: [Swarm, x, y, z, Scoring]
    """
    rank_file = Path(rank_file)
    if not rank_file.exists():
        return pd.DataFrame()

    coord_pattern = re.compile(r'\(\s*([+-]?\d+\.?\d*)\s*,\s*([+-]?\d+\.?\d*)\s*,\s*([+-]?\d+\.?\d*)')

    rows = []
    try:
        with open(rank_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("Swarm"):
                    continue
                swarm_match = re.match(r'\s*(\d+)', line)
                coord_match = coord_pattern.search(line)
                score_match = re.search(r'(\S+)\s*$', line)
                if not (swarm_match and coord_match and score_match):
                    continue
                try:
                    rows.append({
                        "Swarm":   int(swarm_match.group(1)),
                        "x":       float(coord_match.group(1)),
                        "y":       float(coord_match.group(2)),
                        "z":       float(coord_match.group(3)),
                        "Scoring": float(score_match.group(1)),
                    })
                except (ValueError, IndexError):
                    continue
    except Exception:
        return pd.DataFrame()

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.loc[df.groupby("Swarm")["Scoring"].idxmax()].reset_index(drop=True)
    return df


# ============================================================================
# 9) JOB — 잡 디렉토리 상태 분석
# ============================================================================

def get_job_summary(job_dir: Path) -> dict:
    """
    잡 디렉토리의 파이프라인 완료 상태를 분석해 요약 dict를 반환한다.

    Returns
    -------
    dict: {
      "steps": [{"name", "status", "result_path", "log_path", "timing_key", "step_num"}, ...],
      "timing": dict
    }
    """
    job_dir = Path(job_dir)

    step_defs = [
        ("Setup",    job_dir / "setup.json",
         job_dir / "step2_setup.log",     "step2", 2),
        ("시뮬레이션", job_dir / "swarm_0" / "gso_100.out",
         job_dir / "step3_simulation.log", "step3", 3),
        ("포즈 추출", job_dir / "swarm_0" / "lightdock_0.pdb",
         next(iter(sorted(job_dir.glob("step4_generate_swarm_*.log"))), None),
         "step4", 4),
        ("랭킹 생성", job_dir / "rank_by_scoring.list",
         job_dir / "step5_rank.log",      "step5", 5),
    ]

    steps = []
    for name, result_path, log_path, timing_key, step_num in step_defs:
        result_ok = result_path.exists()
        log_ok = (log_path.exists() if isinstance(log_path, Path) else False)
        if result_ok:
            status = "done"
        elif log_ok:
            status = "error"
        else:
            status = "pending"
        steps.append({
            "name": name, "status": status,
            "result_path": result_path, "log_path": log_path,
            "timing_key": timing_key, "step_num": step_num,
        })

    timing = {}
    tf = job_dir / "step_timing.json"
    if tf.exists():
        try:
            timing = json.loads(tf.read_text())
        except Exception:
            pass

    return {"steps": steps, "timing": timing}


# ============================================================================
# 10) PIPELINE — LightDock 도킹 파이프라인 오케스트레이터
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


def run_bsas_clustering(
    job_dir: Path,
    rmsd_cutoff: float = 4.0,
    progress_cb=None,
) -> dict:
    """
    모든 swarm_N 디렉토리에서 lgd_cluster_bsas.py를 실행한다.

    Returns
    -------
    dict: {"total": int, "failed": list[str], "skipped": list[str]}
    """
    job_dir = Path(job_dir)
    swarm_dirs = iter_swarm_dirs(job_dir)
    total = len(swarm_dirs)
    failed = []
    skipped = []

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

    Returns
    -------
    dict: {"success": bool, "steps": {...}, "warnings": list, "errors": list}
    """
    job_dir    = Path(job_dir)
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
    _start2 = _now_kst()

    if resume_mode and (job_dir / "setup.json").exists():
        _progress("step2", "skip")
    else:
        _progress("step2", "running")

        for f in job_dir.glob("lightdock_*.pdb"):
            f.unlink()
        for d in [job_dir / "init", *job_dir.glob("swarm_*")]:
            if d.exists():
                shutil.rmtree(d)
        old_setup = job_dir / "setup.json"
        if old_setup.exists():
            old_setup.unlink()

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

        cmd2 = ["lightdock3_setup.py", actual_rec, actual_lig, "--noxt", "--noh", "--now"]
        if anm_option:
            cmd2.append("-anm")
        cmd2.extend(restraints_arg)

        rc2 = run_subprocess(cmd2, cwd=job_dir,
                             log_path=job_dir / "step2_setup.log",
                             on_line=log_cb)
        _end2 = _now_kst()
        save_timing(job_dir, "step2", _start2, _end2)
        dur2 = format_duration_min(f"{(_end2 - _start2).total_seconds()}s")

        if not (job_dir / "setup.json").exists():
            result["errors"].append("Setup 실패. step2_setup.log를 확인하세요.")
            _progress("step2", "error", dur2)
            return result
        _progress("step2", "done", dur2)

    # ── STEP 3: GSO 시뮬레이션 ───────────────────────────────────────────────
    _start3 = _now_kst()
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
        _end3 = _now_kst()
        save_timing(job_dir, "step3", _start3, _end3)
        _progress("step3", "done",
                  format_duration_min(f"{(_end3 - _start3).total_seconds()}s"))

    # ── STEP 4: 포즈 추출 + BSAS ─────────────────────────────────────────────
    _start4 = _now_kst()
    if resume_mode and _step4_done(job_dir):
        _progress("step4", "skip")
    else:
        _progress("step4", "running")

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

        _end4 = _now_kst()
        save_timing(job_dir, "step4", _start4, _end4)
        _progress("step4", "done",
                  format_duration_min(f"{(_end4 - _start4).total_seconds()}s"))

    if not (resume_mode and _clustering_done(job_dir)):
        run_bsas_clustering(job_dir)

    # ── STEP 5: 랭킹 ────────────────────────────────────────────────────────
    _start5 = _now_kst()
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
        _end5 = _now_kst()
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
            filter_rc = run_subprocess(
                ["lgd_filter_restraints.py", "rank_by_scoring.list",
                 "restraints.list", "A", "B", "--fnat", "0.0"],
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


# ============================================================================
# 11) CLI — argparse 기반 엔트리포인트
# ============================================================================

if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="LightDock + DNAWorks 통합 도킹 파이프라인 (단일 파일 버전 v14, KST)"
    )
    parser.add_argument("--job-dir",    required=True, type=Path, help="작업 결과 디렉토리")
    parser.add_argument("--inputs-dir", required=True, type=Path, help="입력 PDB 디렉토리")
    parser.add_argument("--rec-name",   required=True, help="수용체 PDB 파일명")
    parser.add_argument("--lig-name",   required=True, help="리간드 PDB 파일명")
    parser.add_argument("--no-pdb-cleaning", action="store_true",
                        help="PDB 클리닝(HETATM 제거)을 건너뜀")
    parser.add_argument("--anm", action="store_true",
                        help="lightdock3_setup.py에 -anm 플래그 추가 (ProDy 필요)")
    parser.add_argument("--substitute-nonstandard", action="store_true",
                        help="비표준 아미노산을 표준으로 치환")
    parser.add_argument("--use-pdbfixer", action="store_true",
                        help="PDBFixer로 누락 원자 보완")
    parser.add_argument("--use-restraints", action="store_true",
                        help="--hotspot-input 기반 restraints.list 생성")
    parser.add_argument("--hotspot-input", default="",
                        help='핫스팟 원자 번호 문자열 (예: "A425, A426, B12")')
    parser.add_argument("--hotspot-role", choices=["R", "L"], default="R",
                        help='핫스팟 역할: R(수용체) 또는 L(리간드)')
    parser.add_argument("--aws", action="store_true",
                        help="AWS 환경용 코어 수(전체-1) 사용")
    args = parser.parse_args()

    def _progress(step_id, status, duration=None):
        msg = f"[{step_id}] {status}"
        if duration:
            msg += f" ({duration})"
        print(msg, flush=True)

    def _log(line):
        print(line, flush=True)

    def _warn(msg):
        print(f"[WARN] {msg}", flush=True)

    status = get_tool_status()
    print(f"[env] lightdock={status['lightdock']} "
          f"dnaworks={status['dnaworks']} "
          f"pymol={status['pymol']}", flush=True)

    result = run_full_docking_pipeline(
        job_dir=args.job_dir,
        inputs_dir=args.inputs_dir,
        rec_name=args.rec_name,
        lig_name=args.lig_name,
        pdb_cleaning=not args.no_pdb_cleaning,
        anm_option=args.anm,
        substitute_nonstandard=args.substitute_nonstandard,
        use_pdbfixer=args.use_pdbfixer,
        use_restraints=args.use_restraints,
        hotspot_input=args.hotspot_input,
        hotspot_role=args.hotspot_role,
        is_aws=args.aws,
        progress_cb=_progress,
        log_cb=_log,
        warn_cb=_warn,
    )

    if result["warnings"]:
        print(f"\n[warnings] {len(result['warnings'])} warning(s)")
        for w in result["warnings"]:
            print(f"  - {w}")
    if result["errors"]:
        print(f"\n[errors]")
        for e in result["errors"]:
            print(f"  - {e}")

    sys.exit(0 if result["success"] else 1)
