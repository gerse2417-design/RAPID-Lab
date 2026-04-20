"""
amr_dock_backend/analysis.py
──────────────────────────────────────────────────────────────
도킹 결과 과학 분석 유틸리티 (Streamlit 의존 없음).

책임 범위:
  - 복합체 PDB의 체인 간 인터페이스 잔기 탐지 (BioPython 기반)
  - rank_by_scoring.list + cluster.repr 병합으로 에너지 지형 DataFrame 구성
  - rank_by_scoring.list에서 스왐별 GSO 좌표 추출 (Heatmap용)

의존성:
  - 외부: BioPython (Bio.PDB), pandas, re
  - 내부: amr_dock_backend.parsers, amr_dock_backend.swarm_utils
"""

import re
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser

from .parsers import parse_rank_list, parse_cluster_repr
from .swarm_utils import collect_all_swarm_clusters


def find_interface_residues(complex_pdb: Path, cutoff: float = 4.0) -> dict:
    """
    BioPython으로 복합체 PDB의 체인 간 거리 < cutoff Å인 잔기 쌍을 탐지한다.

    Parameters
    ----------
    complex_pdb : 복합체 PDB 파일 경로 (chain A = receptor, chain B = ligand)
    cutoff      : 인터페이스 거리 기준 (Å), 기본 4.0

    Returns
    -------
    dict: {
      "chains": [chain_id, ...],
      "pairs":  [(chain_A, resid_A, resname_A, chain_B, resid_B, resname_B), ...]
    }
    오류 발생 시 {"chains": [], "pairs": [], "error": str} 반환.
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

    각 포즈에 해당 클러스터 ID와 Population을 조인하여
    에너지 지형 산점도(RMSD vs Scoring) 데이터를 준비한다.

    컬럼: Swarm, Glowworm, RMSD, Scoring, PDB, ClusterID, Population
    rank 파일이 없거나 비어 있으면 빈 DataFrame 반환.
    """
    rank_file = Path(job_dir) / "rank_by_scoring.list"
    df_rank = parse_rank_list(rank_file)
    if df_rank.empty:
        return pd.DataFrame()

    # cluster.repr 전체 수집 → {(swarm_idx, glowworm_id): (cluster_id, population)}
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

    Coordinates 컬럼 형식: (x, y, z, r1, r2, ...) — 앞 3개가 translation 좌표.

    Returns
    -------
    pd.DataFrame — 컬럼: [Swarm, x, y, z, Scoring]
    Swarm 당 최고 Scoring인 행 1개만 반환.
    파싱 실패 시 빈 DataFrame 반환.
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
