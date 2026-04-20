# InhibiForge — Stage 3 (억제 단백질 설계 파이프라인)

RFdiffusion → ProteinMPNN → AlphaFold2-multimer → LightDock → DNAWorks 까지
단일 Streamlit UI / CLI 에서 실행하는 통합 이미지입니다.

## Docker Hub 이미지

이 디렉토리 소스로 빌드·배포된 공식 이미지:

**`syparkbioinfo/rapid-lab:inhibiforge-v3.1`** → https://hub.docker.com/r/syparkbioinfo/rapid-lab/tags

```bash
docker pull syparkbioinfo/rapid-lab:inhibiforge-v3.1
```

v3.1 변경점(v3.0 대비):
- RAPID Lab 허브 연동 "🏠 서비스 홈" 버튼 (환경변수 `HUB_URL`)
- RFdiffusion 후보군 갯수 상한(`max_value=50`) 제거 — 무제한 입력
- LightDock 결과 페이지 "🔍 모든 포즈 보기 (fnat 필터 무시)" 토글 추가

## 빠른 실행

### A) docker run 단독
```bash
docker run --rm --gpus all -p 8501:8501 \
  --shm-size=16gb \
  -e RUN_MODE=ui \
  -e HUB_URL=http://localhost:8502 \
  -v /path/to/RFdiffusion:/mnt/ebs/rfdiffusion \
  -v /path/to/DNAWorks:/mnt/ebs/dnaworks:ro \
  -v /path/to/Projects:/projects \
  -v /usr/lib/ucsf-chimerax-daily:/usr/lib/ucsf-chimerax-daily:ro \
  syparkbioinfo/rapid-lab:inhibiforge-v3.1
```

### B) docker compose
본 디렉토리의 [`docker-compose.yml`](./docker-compose.yml) 은 InhibiForge 단독 실행용.
RAPID Lab 통합 스택(hub + amr-monitor + actsite + inhibiforge)으로 띄우려면
리포 루트의 통합 compose(별도 설치 파이프라인)를 참고할 것.

```bash
cd local/InhibiForge
docker compose up -d
# http://localhost:8501 접속
```

`RUN_MODE=cli` 로 전환하면 동일 이미지에서 [`app/main.py`](./app/main.py) argparse 기반
파이프라인이 실행됩니다. 자세한 사용법은 [`README_CLI.md`](./README_CLI.md), UI는
[`README_UI.md`](./README_UI.md) 참고.

## 소스에서 재빌드

```bash
cd local/InhibiForge
docker build -f Dockerfile.InhibiForge -t syparkbioinfo/rapid-lab:inhibiforge-v3.2 .
```

이미지 내부 포함:
- Streamlit UI (`ui/`)
- CLI 진입점 (`app/main.py`)
- 도킹 백엔드 (`amr_dock_backend/`, `app/backends/`)
- 파이프라인 스크립트 (`run_make_ari_v5.py`, `run_lightdock_dnaworks_v14.py`, `run_af2.py`, `run_mpnn.py`)
- Python 3.10 / torch 2.1 / DGL 1.1 / JAX 0.4.30 / OpenMM / lightdock

이미지 내 **미포함** (호스트에서 마운트):
- RFdiffusion 모델 가중치 (`/mnt/ebs/rfdiffusion`)
- ColabFold AF2 params + JIT cache (첫 실행 시 자동 다운로드)
- DNAWorks Fortran 바이너리 (`/mnt/ebs/dnaworks`)
- ChimeraX Daily (`/usr/lib/ucsf-chimerax-daily`)

## 환경변수 요약

| Key | 기본값 | 용도 |
|---|---|---|
| `RUN_MODE` | `ui` | `ui` 또는 `cli` |
| `HUB_URL` | `http://localhost:8502` | 사이드바 "서비스 홈" 버튼 링크 |
| `AF2_PARAM_DIR` | `/mnt/ebs/rfdiffusion/af2_cache/colabfold` | AF2 파라미터 디렉토리 |
| `AF2_JAX_CACHE_DIR` | `/tmp/jax_cache` | JAX 컴파일 캐시 경로 (shm_size 16GB 요구) |
| `DNAWORKS_BIN` | `/mnt/ebs/dnaworks/dnaworks` | DNAWorks 실행 파일 |
| `RFDIFF_BASE_DIR` | `/mnt/ebs/rfdiffusion` | RFdiffusion 루트 |

## 관련 링크

- Docker Hub: https://hub.docker.com/r/syparkbioinfo/rapid-lab
- Upstream RFdiffusion: https://github.com/RosettaCommons/RFdiffusion
- Upstream ColabFold: https://github.com/sokrypton/ColabFold
- Upstream LightDock: https://github.com/lightdock/lightdock
- DNAWorks: https://hpcwebapps.cit.nih.gov/dnaworks/
