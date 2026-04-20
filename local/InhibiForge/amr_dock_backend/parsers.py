"""
amr_dock_backend/parsers.py
──────────────────────────────────────────────────────────────
LightDock / DNAWorks 결과 파일 파싱 유틸리티.

책임 범위:
  - rank_by_scoring.list / rank_by_luciferin.list 등 LightDock 랭킹 파일 → DataFrame
  - swarm_N/cluster.repr (BSAS 클러스터링 결과) → DataFrame
  - filtered/rank_filtered.list (fnat 필터링 결과) → DataFrame
  - DNAWorks LOGFILE.txt (TRIAL·올리고·FINAL SUMMARY) → (trial_oligos, summary_rows)

의존성:
  - 외부: pandas, re
  - 내부: 없음
"""

import re
from pathlib import Path

import pandas as pd


def parse_rank_list(rank_file: Path) -> pd.DataFrame:
    """
    LightDock 랭킹 파일(rank_by_scoring.list 등)을 DataFrame으로 변환한다.

    좌표 컬럼 '(x,y,z,...)' 내 공백이 있어 단순 split 시 컬럼 수가 달라지므로,
    괄호 구간을 'COORDS' 토큰으로 치환한 뒤 파싱한다.

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

    cluster.repr 실제 포맷 (콜론 구분):
        ClusterID:Population:BestScore:BestGlowworm:RepresentativePDB
        예: 0:18:22.65923:159:lightdock_159.pdb

    반환 컬럼: Cluster ID, Population, Best Score, Best Glowworm, Representative PDB
    파일이 없거나 파싱 실패 시 빈 DataFrame 반환.
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

    포맷: swarm_N_G.pdb  <Scoring>  <fnat>
    반환 컬럼: Swarm, Glowworm, Scoring, fnat
    파일이 없거나 파싱 실패 시 빈 DataFrame 반환.
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
      summary_rows : list  [{"Trial", "Tm", "Len", "Score", "TmRange",
                              "Short", "Long", "#Oligos", "#Repeat", "#Misprime"}, ...]
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
            # TRIAL 번호 감지
            trial_match = re.search(r'PARAMETERS FOR TRIAL\s+(\d+)', line)
            if trial_match:
                current_trial = int(trial_match.group(1))
                is_parsing_oligo = False
                continue

            # 올리고 리스트 시작
            if "oligonucleotides need to be synthesized" in line:
                is_parsing_oligo = True
                if current_trial not in trial_oligos:
                    trial_oligos[current_trial] = []
                continue

            # FINAL SUMMARY 시작
            if "FINAL SUMMARY FOR" in line:
                is_parsing_oligo = False
                is_parsing_summary = True
                continue

            # 올리고 파싱 중
            if is_parsing_oligo:
                stripped = line.strip()
                if "FINAL SUMMARY" in line or stripped == "":
                    is_parsing_oligo = False
                    continue
                if stripped.startswith("-"):
                    continue
                parts = stripped.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    try:
                        trial_oligos.setdefault(current_trial, []).append({
                            "ID": int(parts[0]),
                            "Sequence (DNA)": parts[1],
                            "Length": int(parts[-1])
                        })
                    except (ValueError, IndexError):
                        continue
                continue

            # FINAL SUMMARY 점수 행 파싱
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

    except Exception as e:
        # 부분적으로 파싱된 결과라도 반환 (완전히 포기하지 않음)
        if not summary_rows:
            return {}, []
        # Score 오름차순 정렬, Score가 동일(특히 0.000)할 경우 TmRange 오름차순 2차 정렬
        summary_rows.sort(key=lambda r: (r["Score"], r["TmRange"]))
        return trial_oligos, summary_rows

    # Score 오름차순 정렬, Score가 동일(특히 0.000)할 경우 TmRange 오름차순 2차 정렬
    summary_rows.sort(key=lambda r: (r["Score"], r["TmRange"]))
    return trial_oligos, summary_rows
