#!/bin/bash

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE=$1
OUTPUT_DIR="af2_results"
VENV_PATH="$SCRIPT_DIR/af2_env"

if [ -z "$INPUT_FILE" ]; then
    echo "Usage: $0 <input.faa>"
    exit 1
fi

echo "🚀 Starting AlphaFold2 Prediction Pipeline..."
echo "Input File: $INPUT_FILE"
echo "Output Directory: $OUTPUT_DIR (relative to script dir)"

# 1. Run Prediction (Docker)
python3 "$SCRIPT_DIR/af2_runner.py" "$INPUT_FILE" --output "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "❌ Prediction failed."
    exit 1
fi

echo "✅ Prediction finished successfully!"

# 2. Instructions for Visualization
echo ""
echo "🎨 To visualize the results, run the following command:"
echo "----------------------------------------------------"
echo "source $VENV_PATH/bin/activate"
echo "streamlit run $SCRIPT_DIR/visualizer.py --server.port 8501"
echo "----------------------------------------------------"
echo "Then open your browser and go to http://localhost:8501"
