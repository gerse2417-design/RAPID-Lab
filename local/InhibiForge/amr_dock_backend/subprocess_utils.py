"""
amr_dock_backend/subprocess_utils.py
──────────────────────────────────────────────────────────────
순수 subprocess 실행 래퍼 및 타이밍 유틸리티 (Streamlit 의존 없음).

책임 범위:
  - 외부 명령어 실행 + 실시간 라인 콜백 + 로그 파일 저장
  - 단계별 실행 시간을 step_timing.json에 누적 저장
  - 실행 시간 문자열 포맷 변환 ("1442.02s" → "24.03min")

의존성:
  - 외부: subprocess, datetime, json
  - 내부: 없음
"""

import datetime
import json
import subprocess
from pathlib import Path

# KST(UTC+9) 타임존 — step_timing.json 등 모든 타임스탬프를 KST 로 명시.
KST = datetime.timezone(datetime.timedelta(hours=9))


def now_kst() -> datetime.datetime:
    """KST 기준 timezone-aware 현재 시각."""
    return datetime.datetime.now(KST)


def _to_kst(dt: datetime.datetime) -> datetime.datetime:
    """naive 는 KST 로 간주해 tz 를 부여하고, aware 는 KST 로 변환한다."""
    if dt.tzinfo is None:
        return dt.replace(tzinfo=KST)
    return dt.astimezone(KST)


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
