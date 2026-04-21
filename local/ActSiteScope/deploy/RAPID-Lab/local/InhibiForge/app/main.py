#!/usr/bin/env python3
"""
AMR Pipeline CLI Orchestrator
─────────────────────────────────────────────────────────────
Docker / AWS CLI 배포용 통합 파이프라인 실행기.

서브커맨드:
  rfdiff     RFdiffusion + ProteinMPNN + AlphaFold (GPU)
  lightdock  LightDock + DNAWorks 도킹 (CPU)
  all        rfdiff → lightdock 순차 실행

사용 예시:
  python main.py rfdiff --name my_job --contigs "100" --num_designs 5 --step all
  python main.py lightdock --job-dir /data/output/dock --inputs-dir /data/input --rec-name rec.pdb --lig-name lig.pdb
  python main.py all --rfdiff-args "--name test --contigs 100" --lightdock-args "--job-dir /data/output --inputs-dir /data/input --rec-name r.pdb --lig-name l.pdb"
"""

import os
import sys
import shlex
import argparse
from pathlib import Path

from dotenv import load_dotenv

# ── .env 파일 자동 로드 ────────────────────────────────────────────────────
_env_file = os.environ.get("ENV_FILE", "")
if _env_file and Path(_env_file).is_file():
    load_dotenv(_env_file)
elif Path("/config/.env.aws").is_file():
    load_dotenv("/config/.env.aws")
elif Path(Path(__file__).parent.parent / "config" / ".env.aws").is_file():
    load_dotenv(Path(__file__).parent.parent / "config" / ".env.aws")


def run_rfdiff(extra_argv: list[str]):
    """RFdiffusion 파이프라인 실행 (backends/run_make_ari_v4.py)."""
    print("\n" + "=" * 60)
    print("[Pipeline] RFdiffusion 파이프라인 시작")
    print("=" * 60)

    original_argv = sys.argv
    sys.argv = ["run_make_ari_v4.py"] + extra_argv
    try:
        from backends.run_make_ari_v4 import main as rfdiff_main
        rfdiff_main()
    finally:
        sys.argv = original_argv


def run_lightdock(extra_argv: list[str]):
    """LightDock 파이프라인 실행 (backends/run_lightdock_dnaworks_v14.py)."""
    print("\n" + "=" * 60)
    print("[Pipeline] LightDock 도킹 파이프라인 시작")
    print("=" * 60)

    # argparse 직접 파싱하여 run_full_docking_pipeline() 호출
    dock_parser = argparse.ArgumentParser()
    dock_parser.add_argument("--job-dir", required=True, type=Path)
    dock_parser.add_argument("--inputs-dir", required=True, type=Path)
    dock_parser.add_argument("--rec-name", required=True)
    dock_parser.add_argument("--lig-name", required=True)
    dock_parser.add_argument("--no-pdb-cleaning", action="store_true")
    dock_parser.add_argument("--anm", action="store_true")
    dock_parser.add_argument("--substitute-nonstandard", action="store_true")
    dock_parser.add_argument("--use-pdbfixer", action="store_true")
    dock_parser.add_argument("--use-restraints", action="store_true")
    dock_parser.add_argument("--hotspot-input", default="")
    dock_parser.add_argument("--hotspot-role", choices=["R", "L"], default="R")
    dock_parser.add_argument("--aws", action="store_true")
    args = dock_parser.parse_args(extra_argv)

    # AWS 환경변수가 설정되어 있으면 자동으로 --aws 활성화
    if os.environ.get("LIGHTDOCK_AWS", "").lower() in ("true", "1", "yes"):
        args.aws = True

    from backends.run_lightdock_dnaworks_v14 import (
        run_full_docking_pipeline,
        get_tool_status,
    )

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

    if not result["success"]:
        sys.exit(1)


def main():
    _HELP = """\
AMR Pipeline CLI — Docker/AWS 배포용 통합 실행기

사용법:
  main.py {rfdiff,lightdock,all} [인자 ...]

서브커맨드별 사용 예시:
  main.py rfdiff --name my_job --contigs "100" --num_designs 5
  main.py lightdock --job-dir /data/output/dock --inputs-dir /data/input --rec-name rec.pdb --lig-name lig.pdb
  main.py all --rfdiff-args "--name test --contigs 100" --lightdock-args "--job-dir /data/out --inputs-dir /data/in --rec-name r.pdb --lig-name l.pdb"
"""

    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(_HELP)
        sys.exit(0 if len(sys.argv) >= 2 else 1)

    pipeline = sys.argv[1]
    extra_argv = sys.argv[2:]

    if pipeline == "rfdiff":
        run_rfdiff(extra_argv)

    elif pipeline == "lightdock":
        run_lightdock(extra_argv)

    elif pipeline == "all":
        all_parser = argparse.ArgumentParser(description="rfdiff → lightdock 순차 실행")
        all_parser.add_argument("--rfdiff-args", type=str, required=True,
                                help='RFdiffusion 인자 문자열 (예: "--name test --contigs 100")')
        all_parser.add_argument("--lightdock-args", type=str, required=True,
                                help='LightDock 인자 문자열 (예: "--job-dir /data/out --inputs-dir /data/in --rec-name r.pdb --lig-name l.pdb")')
        args = all_parser.parse_args(extra_argv)

        rfdiff_argv = shlex.split(args.rfdiff_args)
        lightdock_argv = shlex.split(args.lightdock_args)

        print("[Pipeline] 전체 파이프라인 실행: rfdiff → lightdock\n")
        run_rfdiff(rfdiff_argv)
        print("\n[Pipeline] RFdiffusion 완료. LightDock으로 전환합니다.\n")
        run_lightdock(lightdock_argv)
        print("\n[Pipeline] 전체 파이프라인 완료.")

    else:
        print(f"[Error] 알 수 없는 서브커맨드: {pipeline}")
        print("사용 가능: rfdiff, lightdock, all")
        sys.exit(1)


if __name__ == "__main__":
    main()
