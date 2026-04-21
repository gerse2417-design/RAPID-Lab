#!/usr/bin/env bash
# ============================================================================
# InhibiForge 통합 이미지 엔트리포인트
#   RUN_MODE=ui   → Streamlit UI 기동 (기본)
#   RUN_MODE=cli  → app/main.py 호출 (인자는 그대로 전달)
# ============================================================================
set -euo pipefail

# ── ChimeraX 심볼릭 링크 생성 ────────────────────────────────────────────────
#   compose 에서 /usr/lib/ucsf-chimerax-daily 를 마운트. chimerax-daily 명령을
#   /usr/bin 에서 쓸 수 있도록 runtime 에 symlink 생성 ($ORIGIN=/usr/lib/.../bin
#   으로 해석되어야 RPATH 로 libpython3.11.so.1.0 을 찾을 수 있음 — /usr/bin 에
#   실제 파일을 bind-mount 하면 $ORIGIN 이 /usr/bin 으로 잘못 풀려 실패).
if [ -x /usr/lib/ucsf-chimerax-daily/bin/ChimeraX ] && [ ! -e /usr/bin/chimerax-daily ]; then
    ln -s /usr/lib/ucsf-chimerax-daily/bin/ChimeraX /usr/bin/chimerax-daily
fi

MODE="${RUN_MODE:-ui}"

case "${MODE}" in
  ui)
    exec python3 -m streamlit run /app/ui/main.py \
      --server.port "${STREAMLIT_SERVER_PORT:-8501}" \
      --server.address "${STREAMLIT_SERVER_ADDRESS:-0.0.0.0}" \
      --server.headless true
    ;;
  cli)
    exec python3 /app/main.py "$@"
    ;;
  *)
    echo "ERROR: RUN_MODE must be 'ui' or 'cli' (got '${MODE}')" >&2
    exit 1
    ;;
esac
