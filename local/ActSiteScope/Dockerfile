# --- ActSiteScope CLI/GPU Dockerfile ---
# NVIDIA CUDA 12.1 runtime base image
FROM nvidia/cuda:12.1.0-runtime-ubuntu22.04

# 1. System Dependencies
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    python3.10 \
    python3-pip \
    openjdk-17-jre-headless \
    libxrender1 \
    libxext6 \
    libsm6 \
    libgl1-mesa-glx \
    libglib2.0-0 \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# 2. Workspace Setup
WORKDIR /app

# 3. Python Dependency Installation
COPY app/requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# 4. Copy Source Code & Binaries
# We copy everything, but .dockerignore will filter out local data/logs/UI
COPY . .

# 5. Internal Tool Environment Variables
ENV INPUT_DIR=/app/inputs \
    RESULT_DIR=/app/results \
    SCRATCH_DIR=/tmp/hotspot_scratch \
    VINA_PATH=/app/scripts/bin/vina \
    P2RANK_PATH=/app/scripts/bin/p2rank/prank \
    JRE_PATH=/usr/bin/java \
    FPOCKET_PATH=/app/scripts/bin/fpocket/fpocket \
    PDB2PQR_PATH=pdb2pqr \
    APBS_PATH=/app/bin/apbs_bin/bin/apbs \
    MPBIND_DIR=/app/scripts/MPBind \
    PYTHON_EXECUTABLE=python3 \
    PYTHONPATH=/app \
    RUN_MODE=ui

# 6. Set Execution Permissions
RUN chmod +x /app/scripts/bin/vina \
    /app/scripts/bin/p2rank/prank \
    /app/scripts/bin/fpocket/fpocket \
    /app/bin/apbs_bin/bin/apbs \
    /app/app/entrypoint.sh || true

# 7. Pre-create Directories
RUN mkdir -p /app/inputs /app/results /tmp/hotspot_scratch

# 8. Entrypoint & Expose
EXPOSE 8501
ENTRYPOINT ["/bin/bash", "/app/app/entrypoint.sh"]
