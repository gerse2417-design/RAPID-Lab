#!/bin/bash
# ActSiteScope Dual-Mode Entrypoint

# Default RUN_MODE to 'ui' if not set
MODE=${RUN_MODE:-ui}

if [ "$MODE" = "cli" ]; then
    echo "[*] ActSiteScope CLI Mode Active"
    # Pass all remaining arguments to main.py
    exec python3 /app/app/main.py "$@"
else
    echo "[*] ActSiteScope UI Mode Active"
    # Start Streamlit
    exec streamlit run /app/app/streamlit_app.py --server.port 8501 --server.address 0.0.0.0
fi
