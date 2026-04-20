"""
amr_dock_backend/env.py
──────────────────────────────────────────────────────────────
외부 도구 설치 감지 및 실행 환경 상수.

책임 범위:
  - LightDock / DNAWorks / PyMOL 설치 여부 확인
  - 도구 실행 경로 상수 (환경변수 override 지원)
  - CPU 코어 수 계산 (로컬/AWS 모드)

의존성:
  - 외부: subprocess, os, shutil, math
  - 내부: 없음

참고:
  - 경로 상수는 환경변수 DNAWORKS_BIN, PYMOL_BIN으로 override 가능.
  - get_tool_status()는 호출 시점에 1회 감지하므로 모듈 로드 시 부팅 지연 없음.
"""

import math
import os
import shutil
import subprocess
from pathlib import Path

# ── 경로 상수 (환경변수 override 지원) ──────────────────────────────────────

DNAWORKS_EXECUTABLE = Path(
    os.environ.get("DNAWORKS_BIN", "/home/sooyeon/amr/DNAWorks/dnaworks")
)

CHIMERAX_BIN = os.environ.get("CHIMERAX_BIN", "/usr/bin/chimerax-daily")

_PYMOL_FIXED = Path(os.environ.get("PYMOL_BIN", "/usr/bin/pymol"))


# ── CPU 유틸 ─────────────────────────────────────────────────────────────────

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


# ── 도구 설치 감지 ───────────────────────────────────────────────────────────

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

    호출 시점에 subprocess를 실행하므로 모듈 로드 시 부팅 지연이 없다.

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
