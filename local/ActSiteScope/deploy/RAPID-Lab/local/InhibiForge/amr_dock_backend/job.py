"""
amr_dock_backend/job.py
──────────────────────────────────────────────────────────────
LightDock 잡(job) 디렉토리 상태 분석 유틸리티.

책임 범위:
  - 잡 디렉토리를 스캔하여 각 파이프라인 단계의 완료 여부 판단
  - 타이밍 정보 로드

의존성:
  - 외부: json
  - 내부: 없음
"""

import json
from pathlib import Path


def get_job_summary(job_dir: Path) -> dict:
    """
    잡 디렉토리의 파이프라인 완료 상태를 분석해 요약 dict를 반환한다.

    각 단계는 산출물 파일(result_path) 존재 여부로 완료(done)/에러(error)/대기(pending)를 판정.
    타이밍 정보는 step_timing.json에서 로드한다.

    Returns
    -------
    dict: {
      "steps": [
        {"name": str, "status": "done"|"error"|"pending",
         "result_path": Path, "log_path": Path or None,
         "timing_key": str, "step_num": int},
        ...
      ],
      "timing": dict  (step_timing.json 내용, 없으면 {})
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
