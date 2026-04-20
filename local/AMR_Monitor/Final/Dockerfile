# Use NVIDIA CUDA base image for GPU support
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

# Define environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV NUMBA_CACHE_DIR=/tmp

# Update system & install dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    wget \
    git \
    build-essential \
    kalign \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Alias python3 to python
RUN ln -s /usr/bin/python3 /usr/bin/python || true

# Create working directory
WORKDIR /app

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install dependencies from app/requirements.txt
COPY app/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# Install ColabFold directly for AlphaFold2 support
RUN pip install --no-cache-dir "colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold"
RUN pip install --no-cache-dir "jax[cuda12]>=0.4.23" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html || true

# Copy project files
COPY app/ /app/app/
COPY config/ /app/config/
COPY af2_module/ /app/af2_module/
COPY database.py pipeline.py diamond prodigal /app/
# Copy the database optionally or assume mounted/pre-filled in radar folder
COPY radar/ /app/radar/

# Ensure binaries are executable and pipeline directories are writable
RUN chmod +x /app/diamond /app/prodigal || true && \
    mkdir -p /app/radar/pipeline && \
    chmod -R 777 /app/radar/pipeline /app/radar/database || true

# Entrypoint for CLI Execution
ENTRYPOINT ["python", "/app/app/main.py"]
