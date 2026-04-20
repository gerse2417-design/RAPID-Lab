# InhibiForge — Docker CLI 배포 가이드 (v1.1 / Stage 1)

CLI 파이프라인(`app/main.py`)을 컨테이너로 배포하기 위한 가이드.
Streamlit UI까지 포함한 통합 이미지(`inhibiforge-full:v3.0`)는 [README_UI.md](README_UI.md)를 참조.

**Stage 1 변경 요약 (v1.0 → v1.1, 통합은 v2.1 → v3.0)**
- `PIPELINE_SCRIPT` 기본값: `run_make_ari_v4.py` → **`run_make_ari_v5.py`**
- 신규 ENV `LIGHTDOCK_SCRIPT=/app/backends/run_lightdock_dnaworks_v14.py`
- 신규 ENV `INHIBIFORGE_PROJECTS_DIR=/projects` + 동명 볼륨
- `common/amr_paths.py` 모듈 이미지에 포함 (`input_root`/`output_root`/`job_dir`)
- ChimeraX/DNAWorks/pymol 이미지 미포함 — 호스트 설치 + 바이너리 마운트

---

## 디렉토리 구조

```
deploy/inhibiforge-cli/
├── Dockerfile.cli              # ⭐ CLI 통합 이미지 (GPU+CPU, UI 없음) → inhibiforge-cli:v1.1
├── Dockerfile.InhibiForge      # 통합 이미지 (UI+CLI, RUN_MODE 토글)  → inhibiforge-full:v3.0
├── Dockerfile.gpu              # (Legacy) RFdiffusion 단독 이미지
├── Dockerfile.cpu              # (Legacy) LightDock 단독 이미지
├── docker-compose.yml          # UI/CLI 서비스 정의 (full 이미지용)
├── entrypoint.sh               # RUN_MODE=ui|cli 분기
├── .dockerignore
├── README_CLI.md               # ← 본 문서
├── README_UI.md                # 통합 UI 이미지 가이드
│
├── app/                        # CLI 코드 (이미지 안: /app)
│   ├── main.py                 # argparse 진입점 (rfdiff / lightdock / all)
│   ├── requirements-app.txt    # 통합 의존성 (CLI + UI 공용)
│   ├── requirements-gpu.txt    # (Legacy) Dockerfile.gpu 전용
│   ├── requirements-cpu.txt    # (Legacy) Dockerfile.cpu 전용
│   └── backends/
│       ├── _chain_helpers.py                 # 공유 헬퍼 — contigs 파싱 + protocol 감지
│       ├── run_make_ari_v4.py                # (Legacy) RFdiffusion + ProteinMPNN + AF2 (v2.1 유산)
│       ├── run_make_ari_v5.py                # ⭐ Stage 1 RFdiffusion 러너 (Projects 레이아웃, KST)
│       ├── run_mpnn.py                       # ProteinMPNN 단계
│       ├── run_af2.py                        # AlphaFold-Multimer 단계
│       └── run_lightdock_dnaworks_v14.py     # ⭐ Stage 1 LightDock+DNAWorks 러너
│
├── amr_dock_backend/           # LightDock 백엔드 패키지 (이미지 안: /app/amr_dock_backend)
│   └── pipeline.py, job.py, ...
│
├── common/                     # ⭐ Stage 1 공통 모듈 (이미지 안: /app/common)
│   ├── __init__.py
│   └── amr_paths.py            # input_root() / output_root() / ensure_roots() / job_dir()
│
├── ui/                         # Streamlit UI (full 이미지에만 사용)
│   └── lib/
│       ├── pipeline_runner.py       # UI → 백엔드 subprocess + per-job env 주입
│       └── project_cleanup.py       # 성공 시 중간 산출물 정리 (pLDDT Top5 보존)
│
├── config/
│   └── .env.aws                # AWS 자격증명 등 (★ 이미지에 포함되지 않음, --env-file로만 주입)
│
└── test_data/
    ├── input/                  # 샘플 입력 PDB
    └── outputs/                # 결과물 (fallback)
```

---

## 도커 이미지 종류

| 이미지 | Dockerfile | 구성 | 진입점 | 용도 |
|--------|-----------|------|--------|------|
| `inhibiforge-cli:v1.1` | `Dockerfile.cli` | GPU+CPU CLI (UI 없음) | `python3 /app/main.py` | 가벼운 CLI 전용 배포 (서버 자동화, AWS 배치) |
| `inhibiforge-full:v3.0` | `Dockerfile.InhibiForge` | UI + GPU + CPU | `/entrypoint.sh` (`RUN_MODE` 분기) | UI 포함 단일 이미지. CLI 모드도 가능 |
| `amr-rfdiff:latest` | `Dockerfile.gpu` | RFdiffusion 단독 (Legacy) | RFdiffusion 직접 호출 | 분리형 GPU 서버 (구버전) |
| `amr-lightdock:latest` | `Dockerfile.cpu` | LightDock 단독 (Legacy) | LightDock 직접 호출 | 분리형 CPU 서버 (구버전) |

새 배포에서는 `inhibiforge-cli:v1.1` 또는 `inhibiforge-full:v3.0` 사용을 권장.

---

## 빌드

### CLI 통합 이미지 (권장)

```bash
cd /home/sooyeon/amr/deploy/inhibiforge-cli
docker build -f Dockerfile.cli -t inhibiforge-cli:v1.1 .
```

### 통합본(UI 포함)

```bash
docker build -f Dockerfile.InhibiForge -t inhibiforge-full:v3.0 .
```

### Legacy (단독 이미지)

```bash
docker build -f Dockerfile.gpu -t amr-rfdiff:latest .
docker build -f Dockerfile.cpu -t amr-lightdock:latest .
```

---

## 실행 (`inhibiforge-cli:v1.1`)

ENTRYPOINT가 `python3 /app/main.py`로 고정되어 있어, 컨테이너 인자가 그대로 서브커맨드/옵션이 됩니다.

### RFdiffusion → ProteinMPNN → AlphaFold (GPU 단계)

```bash
docker run --rm --gpus all \
  --env-file config/.env.aws \
  -v /home/sooyeon/amr/RFdiffusion:/mnt/ebs/rfdiffusion \
  -v $(pwd)/test_data/input:/data/input \
  -v $(pwd)/test_data/outputs:/data/output \
  -v /home/sooyeon/amr/LightdockDockQ/InhibiForge_app/Projects:/projects \
  inhibiforge-cli:v1.1 \
    rfdiff \
      --name my_job \
      --contigs "100" \
      --num_designs 5 \
      --iterations 50 \
      --step all
```

### LightDock + DNAWorks (CPU 단계)

```bash
docker run --rm \
  --env-file config/.env.aws \
  -v /home/sooyeon/amr/RFdiffusion:/mnt/ebs/rfdiffusion \
  -v $(pwd)/test_data/input:/data/input \
  -v $(pwd)/test_data/outputs:/data/output \
  -v /home/sooyeon/amr/LightdockDockQ/InhibiForge_app/Projects:/projects \
  inhibiforge-cli:v1.1 \
    lightdock \
      --job-dir /data/output/2UUY \
      --inputs-dir /data/input \
      --rec-name 2UUY_rec.pdb \
      --lig-name 2UUY_lig.pdb
```

### 전체 파이프라인 (rfdiff → lightdock 순차)

```bash
docker run --rm --gpus all \
  --env-file config/.env.aws \
  -v /home/sooyeon/amr/RFdiffusion:/mnt/ebs/rfdiffusion \
  -v $(pwd)/test_data/input:/data/input \
  -v $(pwd)/test_data/outputs:/data/output \
  -v /home/sooyeon/amr/LightdockDockQ/InhibiForge_app/Projects:/projects \
  inhibiforge-cli:v1.1 \
    all \
      --rfdiff-args "--name test --contigs 100 --num_designs 5" \
      --lightdock-args "--job-dir /data/output/test --inputs-dir /data/input --rec-name rec.pdb --lig-name lig.pdb"
```

### 도움말

```bash
docker run --rm inhibiforge-cli:v1.1 --help
docker run --rm inhibiforge-cli:v1.1 rfdiff --help
docker run --rm inhibiforge-cli:v1.1 lightdock --help
docker run --rm inhibiforge-cli:v1.1 all --help
```

### 통합 이미지로 동일 동작

`inhibiforge-full:v3.0`도 `RUN_MODE=cli`만 주면 위와 동일하게 동작합니다.

```bash
docker run --rm --gpus all \
  -e RUN_MODE=cli \
  --env-file config/.env.aws \
  -v ... \
  inhibiforge-full:v3.0 \
    rfdiff --name my_job --contigs "100" --num_designs 5
```

---

## tar.gz 배포 (오프라인/AWS 전송)

### 아카이브 위치

```
/home/sooyeon/amr/deploy/
├── inhibiforge-cli:v1.1.tar.gz     (~11 GB)
└── inhibiforge-full:v3.0.tar.gz    (빌드 직후 생성)
```

### 새로 만들기

```bash
docker save inhibiforge-cli:v1.1  | gzip > /home/sooyeon/amr/deploy/inhibiforge-cli:v1.1.tar.gz
docker save inhibiforge-full:v3.0 | gzip > /home/sooyeon/amr/deploy/inhibiforge-full:v3.0.tar.gz
```

### 다른 머신에서 로드

```bash
gunzip -c inhibiforge-cli:v1.1.tar.gz | docker load
docker images | grep inhibiforge      # 등록 확인
```

### S3 경유 배포 예 (tar.gz 방식)

```bash
aws s3 cp /home/sooyeon/amr/deploy/inhibiforge-cli:v1.1.tar.gz s3://my-bucket/images/

# 타깃 EC2에서
aws s3 cp s3://my-bucket/images/inhibiforge-cli:v1.1.tar.gz - | gunzip | docker load
```

> AWS 정식 배포는 **ECR push** 경로를 권장합니다. 자세한 절차는
> [aws/README_AWS.md](aws/README_AWS.md) (준비 예정) 참조.

---

## 환경변수

`.env.aws`(또는 `--env-file`로 주입한 파일)와 이미지 내장 기본값을 함께 사용합니다.
**이미지 안에는 어떤 시크릿도 굽지 않았습니다** — `.env.aws`는 런타임에만 주입됩니다.

### 공통

| 변수 | 설명 | 이미지 기본값 |
|------|------|---------------|
| `ENV_FILE` | `.env` 파일 경로 (자동 로드) | `/config/.env.aws` (있을 때만) |
| `VENV_CMD` | venv 활성화 명령 (Docker에선 빈값) | `""` |
| `LOG_LEVEL` | 로그 레벨 | `INFO` |

### RFdiffusion (GPU)

| 변수 | 설명 | 이미지 기본값 |
|------|------|---------------|
| `PIPELINE_SCRIPT` | UI/CLI가 호출할 RFdiffusion 러너 (Stage 1) | `/app/backends/run_make_ari_v5.py` |
| `RFDIFF_BASE_DIR` | RFdiffusion 베이스 (볼륨) | `/mnt/ebs/rfdiffusion` |
| `RFDIFF_CODE_DIR` | RFdiffusion 코드 디렉토리 | `${RFDIFF_BASE_DIR}/RFdiffusion` |
| `RFDIFF_INPUT_DIR` | 입력 PDB (fallback) | `/data/input` |
| `RFDIFF_OUTPUT_DIR` | 출력 디렉토리 (fallback) | `/data/output` |
| `COLABDESIGN_DIR` | ColabDesign 리포 | `${RFDIFF_BASE_DIR}/ColabDesign_repo` |
| `AF2_CACHE_DIR` | AlphaFold 캐시 | `${RFDIFF_BASE_DIR}/af2_cache` |
| `AF2_PARAM_DIR` | AlphaFold 파라미터 | `${AF2_CACHE_DIR}/colabfold` |
| `CUSTOM_LIB_DIR` | 커스텀 공유 라이브러리 | `/mnt/ebs/lib` |
| `GPU_MEMORY_FRACTION` | GPU VRAM 사용 비율 | `0.7` |
| `E3NN_JIT_COMPILE` | e3nn JIT 비활성화 | `0` |

### Stage 1 공통 (Projects 레이아웃)

| 변수 | 설명 | 이미지 기본값 |
|------|------|---------------|
| `INHIBIFORGE_PROJECTS_DIR` | 프로젝트 루트 (호스트 `Projects/` → 컨테이너 `/projects`) | `/projects` |
| `AMR_WORKSPACE` | per-job 워크스페이스 (런처가 주입) | (런타임) |
| `AMR_INPUT_DIR` | per-job 입력 (런처가 주입) | (런타임) |
| `AMR_OUTPUT_DIR` | per-job 출력 (런처가 주입) | (런타임) |
| `INHIBIFORGE_PROJECT_DIR` | per-job 프로젝트 디렉토리 (런처가 주입) | (런타임) |
| `PYTHONPATH` | `common.amr_paths` import 경로 포함 | `...:/app` |

`AMR_*` / `INHIBIFORGE_PROJECT_DIR`는 `ui/lib/pipeline_runner.py`가 파이프라인
실행 시마다 주입합니다. 사용자가 직접 설정할 필요는 없지만 CLI에서 Stage 1
러너를 단독 실행할 때는 `common.amr_paths.ensure_roots()`와 호환되도록
`AMR_WORKSPACE`를 명시적으로 주면 안전합니다.

### LightDock (CPU)

| 변수 | 설명 | 이미지 기본값 |
|------|------|---------------|
| `LIGHTDOCK_SCRIPT` | UI/CLI가 호출할 LightDock+DNAWorks 러너 (Stage 1) | `/app/backends/run_lightdock_dnaworks_v14.py` |
| `LD_BASE_DIR` | LightDock 백엔드 루트 (`amr_dock_backend` import 경로) | `/app` |
| `LIGHTDOCK_INPUT_DIR` | 입력 PDB (fallback) | `/data/input` |
| `LIGHTDOCK_OUTPUT_DIR` | 출력 디렉토리 (fallback) | `/data/output` |
| `DNAWORKS_BIN` | DNAWorks 바이너리 (호스트 설치 + 마운트) | `/mnt/ebs/dnaworks/dnaworks` |
| `CHIMERAX_BIN` | ChimeraX 실행 파일 (호스트 설치 + 마운트) | `/usr/bin/chimerax` |
| `PYMOL_BIN` | PyMOL 실행 파일 (호스트 설치 + 마운트) | `/usr/bin/pymol` |
| `LIGHTDOCK_AWS` | AWS 코어 계산 모드 | `true` |

> **중요**: `DNAWORKS_BIN` / `CHIMERAX_BIN` / `PYMOL_BIN`은 이미지에 **포함되지 않은**
> 호스트 바이너리입니다. 호스트에 설치하고 `docker-compose.yml`의 볼륨 마운트
> 주석을 해제하거나 `-v` 옵션으로 주입하세요. 자세한 설치법은
> [README_UI.md](README_UI.md)의 "ChimeraX 설치 가이드" 섹션 참조.

> **참고 — multi-chain 자동 처리**
> Step 2/3 (`run_mpnn.py`, `run_af2.py`)은 `--contigs` 인자(오케스트레이터가 자동 전달)를
> [`app/backends/_chain_helpers.py`](app/backends/_chain_helpers.py)로 파싱해
> **binder / partial / fixbb 프로토콜을 자동 감지**합니다. RFdiffusion이 출력한
> 다중 체인 PDB(target+binder 등)도 정상 처리되고, 이전에 발생하던
> `ValueError: Only single chain PDBs are supported` 문제가 해결됩니다.

### 통합 이미지(full) 전용

| 변수 | 설명 | 이미지 기본값 |
|------|------|---------------|
| `RUN_MODE` | `ui` 또는 `cli` | `ui` |
| `STREAMLIT_SERVER_PORT` | Streamlit 포트 | `8501` |
| `STREAMLIT_SERVER_ADDRESS` | Streamlit 바인드 | `0.0.0.0` |

---

## S3 데이터 연동

### 입력 다운로드

```bash
aws s3 sync s3://my-bucket/input ./test_data/input
```

### 출력 업로드

```bash
aws s3 sync ./test_data/outputs s3://my-bucket/output/${JOB_NAME}
```

---

## EBS 볼륨 마운트 가이드

### GPU 서버 EBS 구조

```
/mnt/ebs/rfdiffusion/
├── RFdiffusion/          # RFdiffusion 코드 (git clone)
├── ColabDesign_repo/     # ColabDesign (git clone)
├── af2_cache/
│   └── colabfold/        # AlphaFold 파라미터 (~5 GB)
└── ...

/mnt/ebs/lib/             # 커스텀 공유 라이브러리
```

### CPU 서버 EBS 구조

```
/mnt/ebs/dnaworks/
└── dnaworks              # DNAWorks 바이너리 (Fortran 컴파일)
```

### EBS 마운트 (EC2)

```bash
sudo mkdir -p /mnt/ebs
sudo mount /dev/xvdf /mnt/ebs

# 재부팅 시 자동 마운트
echo '/dev/xvdf /mnt/ebs ext4 defaults,nofail 0 2' | sudo tee -a /etc/fstab
```

---

## 트러블슈팅

- **`ValueError: Only single chain PDBs are supported when chain_id not specified`**
  Step 2 (ProteinMPNN) 또는 Step 3 (AF2)에서 발생. `_chain_helpers.py` 도입 이전 빌드 잔재.
  최신 `inhibiforge-cli:v1.1` (또는 `inhibiforge-full:v3.0`)로 재빌드 후 다시 실행하세요.
  ```bash
  cd /home/sooyeon/amr/deploy/inhibiforge-cli
  docker build -f Dockerfile.cli -t inhibiforge-cli:v1.1 .
  ```
  이미지가 `_chain_helpers.py`를 포함했는지 빠른 확인:
  ```bash
  docker run --rm --entrypoint ls inhibiforge-cli:v1.1 /app/backends/_chain_helpers.py
  ```

- **`ModuleNotFoundError: No module named 'common'`** (Stage 1 전용)
  `common/` 모듈이 이미지에 포함되지 않았거나 `PYTHONPATH`에 `/app`이 빠졌을 때 발생.
  ```bash
  docker run --rm --entrypoint ls inhibiforge-cli:v1.1 /app/common/amr_paths.py
  docker run --rm --entrypoint printenv inhibiforge-cli:v1.1 PYTHONPATH
  ```

- **ChimeraX / DNAWorks / pymol `Command not found`**
  호스트 설치본이 컨테이너에 마운트되지 않은 경우. 호스트에 각 도구를 설치한 뒤
  `docker-compose.yml`의 바이너리 마운트 주석을 해제하거나 `-v`로 직접 주입.
  [README_UI.md](README_UI.md)의 "ChimeraX 설치 가이드" 섹션 참고.

- **AF2 파라미터 디렉토리 (`--param_dir`) 누락**
  `mk_af_model`이 항상 호출되므로 binder가 아닌 fixbb 모드에서도 `param_dir`가 유효해야 합니다.
  볼륨 마운트(`-v /mnt/ebs/rfdiffusion:/mnt/ebs/rfdiffusion`)에 `af2_cache/colabfold/`가 들어있는지 확인.

---

## 관련 문서

- [README_UI.md](README_UI.md) — `inhibiforge-full:v3.0` (UI 포함) 운영 가이드, `RUN_MODE` 토글, docker-compose 사용법
- `Dockerfile.cli`, `Dockerfile.InhibiForge` — 빌드 스펙
- `entrypoint.sh` — UI/CLI 모드 분기 로직
- [`app/backends/_chain_helpers.py`](app/backends/_chain_helpers.py) — contigs 파싱 + 프로토콜 감지 헬퍼
- [`app/backends/run_make_ari_v5.py`](app/backends/run_make_ari_v5.py) — Stage 1 RFdiffusion 러너
- [`app/backends/run_lightdock_dnaworks_v14.py`](app/backends/run_lightdock_dnaworks_v14.py) — Stage 1 LightDock+DNAWorks 러너
- [`common/amr_paths.py`](common/amr_paths.py) — Stage 1 경로 해석
- **AWS 배포**: [aws/README_AWS.md](aws/README_AWS.md) (준비 예정) — ECR/EFS/EC2/SSM
