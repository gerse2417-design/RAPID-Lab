"""공통 inputs/outputs 루트 해석 헬퍼.

RFdiffusion 계열 파이프라인과 LightDock/DNAWorks 파이프라인이 동일한
inputs/outputs 디렉터리를 공유하기 위한 최소한의 유틸리티.

기본값은 ``/home/sooyeon/amr/data`` 를 루트로 사용하고,
``AMR_WORKSPACE`` / ``AMR_INPUT_DIR`` / ``AMR_OUTPUT_DIR`` 환경변수로 오버라이드할 수 있다.
"""

from __future__ import annotations

import os
from pathlib import Path

DEFAULT_WORKSPACE = Path("/home/sooyeon/amr/data")


def workspace_root() -> Path:
    return Path(os.environ.get("AMR_WORKSPACE", str(DEFAULT_WORKSPACE)))


def input_root() -> Path:
    return Path(os.environ.get("AMR_INPUT_DIR", str(workspace_root() / "inputs")))


def output_root() -> Path:
    return Path(os.environ.get("AMR_OUTPUT_DIR", str(workspace_root() / "outputs")))


def job_dir(job_name: str) -> Path:
    """신규 job 의 출력 디렉터리를 반환. 없으면 생성한다."""
    d = output_root() / job_name
    d.mkdir(parents=True, exist_ok=True)
    return d


def ensure_roots() -> None:
    """inputs/outputs 루트가 존재하지 않으면 생성."""
    input_root().mkdir(parents=True, exist_ok=True)
    output_root().mkdir(parents=True, exist_ok=True)
