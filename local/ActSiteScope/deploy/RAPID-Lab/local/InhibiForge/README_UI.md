# InhibiForge — UI 배포 가이드 (v3.0)

통합 이미지 `inhibiforge-full:v3.0`은 동일 컨테이너에서 **Streamlit UI**와
**CLI 파이프라인**을 모두 제공합니다. 환경변수 `RUN_MODE`로 모드를 전환합니다.

**Stage 1 변경 요약 (v2.1 → v3.0)**
- RFdiffusion 러너: `run_make_ari_v4.py` → **`run_make_ari_v5.py`** (Projects 레이아웃, KST)
- LightDock 러너: v13 → **`run_lightdock_dnaworks_v14.py`**
- 모든 프로젝트 I/O 단일 폴더 통일: `Projects/<프로젝트명>/`, 컨테이너 경로 `/projects`
- `common/amr_paths.py` 모듈 이미지에 포함 (`input_root()`, `output_root()`, `job_dir()`)
- ChimeraX/DNAWorks/pymol은 **이미지에 포함되지 않음** — 호스트 설치 후 바이너리 마운트

| `RUN_MODE` | 동작 | 진입점 |
|------------|------|--------|
| `ui` (기본) | Streamlit 서버 기동, 8501 노출 | `streamlit run /app/ui/main.py` |
| `cli` | `app/main.py` argparse 호출 | `python3 /app/main.py "$@"` |

CLI 전용 이미지(`inhibiforge-cli:v1.1`)는 그대로 별도 유지되며, 통합 이미지의
`RUN_MODE=cli`는 동일 코드를 한 컨테이너에서 같이 쓰고 싶을 때 사용합니다.

---

## 1. 빠른 시작 (UI)

```bash
cd /home/sooyeon/amr/deploy/inhibiforge-cli
docker compose up -d
```

브라우저에서 **http://localhost:8501** 접속.

종료:
```bash
docker compose down
```

---

## 2. 빠른 시작 (CLI 모드)

`compose run`으로 일회성 실행:
```bash
cd /home/sooyeon/amr/deploy/inhibiforge-cli
docker compose run --rm inhibiforge-cli \
    rfdiff --name my_job --contigs "100" --num_designs 5 --step all
```

또는 `docker run` 직접:
```bash
docker run --rm --gpus all \
  -e RUN_MODE=cli \
  --env-file config/.env.aws \
  -v /home/sooyeon/amr/RFdiffusion:/mnt/ebs/rfdiffusion \
  -v $(pwd)/test_data/input:/data/input \
  -v $(pwd)/test_data/outputs:/data/output \
  -v /home/sooyeon/amr/LightdockDockQ/InhibiForge_app/Projects:/projects \
  inhibiforge-full:v3.0 \
    lightdock --job-dir /data/output/my_job \
              --inputs-dir /data/input \
              --rec-name rec.pdb --lig-name lig.pdb
```

---

## 3. 환경변수

`.env.aws`(또는 `--env-file`로 주입한 파일)와 이미지 내장 기본값을 같이 사용합니다.
이미지 안에는 시크릿이 굽혀있지 않습니다.

| 변수 | 기본값 (이미지) | 의미 |
|------|----------------|------|
| `RUN_MODE` | `ui` | `ui` 또는 `cli` |
| `STREAMLIT_SERVER_PORT` | `8501` | Streamlit 포트 |
| `STREAMLIT_SERVER_ADDRESS` | `0.0.0.0` | Streamlit 바인드 주소 |
| `RFDIFF_BASE_DIR` | `/mnt/ebs/rfdiffusion` | 호스트 RFdiffusion 마운트 위치 |
| `RFDIFF_INPUT_DIR` | `/data/input` | RFdiffusion 입력 (v4 호환 / fallback) |
| `RFDIFF_OUTPUT_DIR` | `/data/output` | RFdiffusion 출력 (v4 호환 / fallback) |
| `PIPELINE_SCRIPT` | `/app/backends/run_make_ari_v5.py` | UI가 호출하는 RFdiffusion 러너 (Stage 1 v5) |
| `LIGHTDOCK_SCRIPT` | `/app/backends/run_lightdock_dnaworks_v14.py` | UI가 호출하는 LightDock+DNAWorks 러너 (Stage 1 v14) |
| `INHIBIFORGE_PROJECTS_DIR` | `/projects` | Stage 1 프로젝트 루트 (`Projects/<name>/` 영속 저장) |
| `LD_BASE_DIR` | `/app` | LightDock 백엔드 루트 (`amr_dock_backend` import 경로) |
| `LIGHTDOCK_INPUT_DIR` | `/data/input` | LightDock 입력 (fallback) |
| `LIGHTDOCK_OUTPUT_DIR` | `/data/output` | LightDock 출력 (fallback) |
| `VENV_CMD` | `""` | Docker에선 빈값 (시스템 Python 사용) |
| `PYTHONPATH` | `...:/app` | `common.amr_paths` import가 동작하도록 `/app` 포함 |

`.env.aws`에는 AWS 자격증명 등만 넣고, 위 경로는 이미지 기본값을 그대로 두면 됩니다.
Stage 1 러너(`run_make_ari_v5.py`, `run_lightdock_dnaworks_v14.py`)는
`pipeline_runner.py`가 per-job으로 `AMR_WORKSPACE`/`AMR_INPUT_DIR`/`AMR_OUTPUT_DIR`/
`INHIBIFORGE_PROJECT_DIR`를 주입하므로 사용자가 직접 설정할 필요는 없습니다.

---

## 4. 볼륨 매핑 정리

| 컨테이너 경로 | 호스트 경로 (compose 기본) | 용도 |
|---------------|---------------------------|------|
| `/mnt/ebs/rfdiffusion` | `/home/sooyeon/amr/RFdiffusion` | RFdiffusion 코드 + 가중치 (대용량, 이미지 외부) |
| `/data/input` | `./test_data/input` | 사용자 입력 PDB / FASTA (fallback) |
| `/data/output` | `./test_data/outputs` | 결과물 (fallback, v4 호환) |
| `/projects` | `/home/sooyeon/amr/LightdockDockQ/InhibiForge_app/Projects` | **Stage 1 프로젝트 영속 저장소**. 컨테이너를 재시작/삭제해도 프로젝트가 유지됨 |
| `/usr/bin/chimerax-daily` (주석) | 호스트 ChimeraX 바이너리 | `split_complex` 단계에서 필요 시 주석 해제 (아래 ChimeraX 섹션 참조) |

### 🛠 ChimeraX 설치 가이드

본 파이프라인의 `split_complex` 단계에서는 ChimeraX가 필수입니다. 아래 공식
링크에서 운영체제에 맞는 최신 버전을 설치해 주세요.

공식 설치 가이드: https://github.com/RBVI/ChimeraX

- 참고 사항: 설치 후 시스템 경로(PATH)에 `chimerax` 명령어가 등록되어 있는지 확인해 주세요.
- 오류 방지: 설치 도중 문제가 발생하거나 파이프라인에서 `Command not found` 에러가
  발생하면, 관리자에게 문의하기 전 환경 변수 설정을 먼저 확인 바랍니다.

**컨테이너에서 호스트 ChimeraX 사용법**: `docker-compose.yml`의

```yaml
# - /usr/bin/chimerax-daily:/usr/bin/chimerax-daily:ro
# - /opt/UCSF:/opt/UCSF:ro
```

두 줄을 주석 해제하면 호스트 설치본이 그대로 컨테이너 안에서 호출됩니다.
배포판에 따라 설치 디렉토리 경로가 다를 수 있으므로 경로를 실제 호스트에
맞춰 조정하세요. DNAWorks / pymol 도 동일한 방식(호스트 설치 + 바이너리 마운트)으로
처리합니다.

### Projects 볼륨 설명

Stage 1 리팩터링으로 **프로젝트 단위 폴더**가 모든 산출물의 루트가 됩니다.
호스트의 `InhibiForge_app/Projects/<project_name>/` 아래에 RFdiffusion 출력,
ProteinMPNN seqs, AF2 PDB, LightDock 결과, pLDDT Top5, DNAWorks 등이 누적됩니다.
컨테이너는 이 디렉토리를 `/projects`로 마운트하므로:

- UI의 "파일 보관함" 페이지가 호스트 프로젝트를 그대로 노출
- 컨테이너 재시작/삭제 후에도 결과 유지
- 여러 호스트/인스턴스가 동일 NFS·EFS·SMB를 마운트하면 프로젝트 공유 가능

---

## 5. CLI ↔ UI 모드 전환 요약

```bash
# UI (기본)
docker compose up -d                          # → :8501
docker compose down

# CLI 일회성
docker compose run --rm inhibiforge-cli rfdiff --name test --contigs 100

# 같은 이미지를 직접 docker run 으로
docker run --rm -e RUN_MODE=ui  ... inhibiforge-full:v3.0
docker run --rm -e RUN_MODE=cli ... inhibiforge-full:v3.0 rfdiff --name test ...
```

---

## 6. 트러블슈팅

- **`ValueError: Only single chain PDBs are supported when chain_id not specified`** (Step 2/3)
  `_chain_helpers.py` 도입 이전 빌드의 잔재. 최신 이미지로 재빌드 후 재실행.
  ```bash
  cd /home/sooyeon/amr/deploy/inhibiforge-cli
  docker build -f Dockerfile.InhibiForge -t inhibiforge-full:v3.0 .
  ```
  헬퍼 존재 여부 빠른 확인: `docker run --rm --entrypoint ls inhibiforge-full:v3.0 /app/backends/_chain_helpers.py`
- **UI는 뜨는데 RFdiffusion v3/v4를 찾는다고 에러**:
  컨테이너 환경변수 `PIPELINE_SCRIPT`가 `/app/backends/run_make_ari_v5.py`로 설정되어 있는지 확인.
  `docker exec <id> printenv PIPELINE_SCRIPT`.
- **`ModuleNotFoundError: No module named 'common'`**:
  `common/` 모듈이 이미지에 포함되지 않았거나 `PYTHONPATH`에 `/app`이 없는 경우.
  ```bash
  docker run --rm --entrypoint ls inhibiforge-full:v3.0 /app/common/amr_paths.py
  docker run --rm --entrypoint printenv inhibiforge-full:v3.0 PYTHONPATH
  ```
- **ChimeraX `Command not found`**:
  호스트에 ChimeraX가 설치되어 있지 않거나 compose의 바이너리 마운트 주석이 해제되지 않음.
  위 **ChimeraX 설치 가이드** 섹션 참조.
- **AF2 (Step 3) `Failed to add kernel node to a CUDA graph: CUDA_ERROR_INVALID_VALUE`**:
  실제 원인은 JAX/CUDA 가 아니라 **Docker 기본 `/dev/shm=64MB` 가 AF2-multimer(binder)
  shared memory 에 부족**해서 발생. 드라이버가 올리는 에러 메시지가
  `cuGraphAddKernelNode` 로 거슬러 올라오면서 오해하기 쉽지만 근본은 shm.
  `docker-compose.yml` 의 `shm_size: "16gb"` 로 해결됨 (v3.0 기본값).
  직접 `docker run` 시엔 `--shm-size=16g` 옵션 필요.
- **DNAWorks / pymol 호출 실패**:
  ChimeraX와 동일 — 호스트 설치 + compose 볼륨 마운트 필요.
- **`/projects` 비어 있음**:
  compose의 `/home/sooyeon/amr/LightdockDockQ/InhibiForge_app/Projects:/projects`
  마운트가 정상인지 `docker inspect` 로 확인. 호스트 경로가 다르면 compose에서 조정.
- **GPU 인식 실패**:
  `docker run --gpus all` 누락 또는 호스트의 nvidia-container-toolkit 미설치.
- **포트 8501 충돌**:
  compose의 `ports: "8501:8501"` 좌측 값을 다른 포트(예: `"18501:8501"`)로 변경.
- **`config/.env.aws` 파일 없음**:
  AWS 인증이 필요 없으면 빈 파일을 만들거나 compose의 `env_file:` 항목을 제거.

## 7. 관련 파일

- [`app/backends/run_make_ari_v5.py`](app/backends/run_make_ari_v5.py) — Stage 1 RFdiffusion 러너
- [`app/backends/run_lightdock_dnaworks_v14.py`](app/backends/run_lightdock_dnaworks_v14.py) — Stage 1 LightDock+DNAWorks 러너
- [`app/backends/_chain_helpers.py`](app/backends/_chain_helpers.py) — `run_mpnn.py`/`run_af2.py`가 공유하는 contigs 파싱 + 프로토콜(binder/partial/fixbb) 감지
- [`common/amr_paths.py`](common/amr_paths.py) — Stage 1 공통 경로 해석 (`input_root`, `output_root`, `job_dir`)
- [`ui/lib/pipeline_runner.py`](ui/lib/pipeline_runner.py) — UI → 백엔드 subprocess 호출, per-job env 주입
- [`ui/lib/project_cleanup.py`](ui/lib/project_cleanup.py) — 성공 시 중간 산출물 정리 (pLDDT Top5 보존)
- [`Dockerfile.InhibiForge`](Dockerfile.InhibiForge) — 통합 이미지 빌드 스펙
- [`entrypoint.sh`](entrypoint.sh) — `RUN_MODE=ui|cli` 분기
- [README_CLI.md](README_CLI.md) — CLI 전용 이미지(`inhibiforge-cli:v1.1`) 가이드, 환경변수 전체 표
- **AWS 배포**: [`aws/README_AWS.md`](aws/README_AWS.md) — ECR/EFS/EC2/SSM 단계별 런북
