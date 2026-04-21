"""
amr_dock_backend/swarm_utils.py
──────────────────────────────────────────────────────────────
LightDock 스왐(swarm) 디렉토리 순회 및 PDB 포즈 목록 수집 유틸리티.

책임 범위:
  - swarm_N 디렉토리를 N 오름차순으로 정렬하여 반환하는 공통 헬퍼
  - 전체 스왐의 cluster.repr 데이터 수집
  - 스왐별 lightdock_*.pdb 포즈 파일 목록 수집

의존성:
  - 외부: pandas, re
  - 내부: amr_dock_backend.parsers (parse_cluster_repr)
"""

import re
from pathlib import Path

import pandas as pd

from .parsers import parse_cluster_repr


def iter_swarm_dirs(job_dir: Path) -> list:
    """
    job_dir 내 swarm_N 디렉토리를 N 오름차순으로 정렬하여 반환한다.

    코드베이스 전반에 반복되던 정렬 패턴을 단일 함수로 통합.
    swarm_N 형식이 아닌 항목은 0으로 처리하여 오류 없이 건너뛴다.
    """
    job_dir = Path(job_dir)
    return sorted(
        (d for d in job_dir.glob("swarm_*") if d.is_dir()),
        key=lambda d: int(re.search(r'swarm_(\d+)', d.name).group(1))
        if re.search(r'swarm_(\d+)', d.name) else 0
    )


def collect_all_swarm_clusters(job_dir: Path) -> dict:
    """
    job_dir 내 모든 swarm_N 디렉토리를 순회하여
    각 스왐의 cluster.repr 데이터를 수집한다.

    Returns
    -------
    dict: {"swarm_0": DataFrame, "swarm_1": DataFrame, ...}
    cluster.repr 가 없는 스왐은 결과에서 제외된다.
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
    PDB 파일이 없는 스왐은 결과에서 제외된다.
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
