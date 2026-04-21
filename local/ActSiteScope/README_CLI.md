# ActSiteScope CLI & Docker User Guide

본 프로젝트는 대시보드 시스템의 핵심 분석 파이프라인을 CLI(Command Line Interface)로 제공하여, 대량의 단백질 구조 데이터를 AWS 클라우드 등 배포 환경에서 효율적으로 처리할 수 있도록 지원합니다.

## 1. 주요 특징
- **통합 분석**: M-CSA, AI-Based Binding(MPBind), P2Rank, fPocket, APBS 데이터의 컨센서스 통합.
- **GPU 가속**: NVIDIA CUDA 12.1 기반 AI 추론(MPBind) 지원.
- **배포 최적화**: AWS S3 입력/출력 및 EFS 연동 지원.
- **UI 제외**: 분석 엔진 자원 최적화를 위해 Streamlit 등 UI 패키지를 제외한 경량화 이미지.

## 2. Docker 빌드 및 실행

### 이미지 빌드
```bash
docker build -t actsitescope-cli .
```

### 실행 예시 (로컬 파일 처리)
```bash
docker run --gpus all \
  -v $(pwd)/inputs:/app/inputs \
  -v $(pwd)/results:/app/results \
  actsitescope-cli \
  --pdb /app/inputs/target_protein.pdb \
  --ligands /app/inputs/drug1.sdf /app/inputs/drug2.mol2
```

### 실행 예시 (S3 파일 처리)
AWS 인증 정보가 환경 변수로 설정되어 있어야 합니다.
```bash
docker run --gpus all \
  -e AWS_ACCESS_KEY_ID=YOUR_KEY \
  -e AWS_SECRET_ACCESS_KEY=YOUR_SECRET \
  -v $(pwd)/results:/app/results \
  actsitescope-cli \
  --pdb s3://my-bucket/inputs/protein.pdb \
  --ligands s3://my-bucket/inputs/cmpd1.sdf
```

## 3. CLI 상세 옵션
| 옵션 | 설명 | 비고 |
| :--- | :--- | :--- |
| `--pdb` | 분석할 단백질 PDB 경로 (로컬 또는 S3) | **필수** |
| `--uniprot` | M-CSA 조회를 위한 UniProt ID | 선택 사항 |
| `--ligands` | 도킹 검증용 리간드 파일 목록 (Space 구분) | 선택 사항 |
| `--env` | 환경 변수 설정 파일 경로 | 기본: `config/.env.aws` |
| `--output` | 결과 저장 경로 오버라이드 | 기본: `.env` 설정값 |
| `--th_mp` | AI 결합 확률 임계값 (0.0~1.0) | 기본: 0.5 |
| `--th_p2r` | P2Rank 포켓 상위 % 기준 | 기본: 30 |

## 4. 환경 변수 설정 (`.env.aws`)
`config/.env.aws` 파일을 통해 분석 파라미터 및 도구 경로를 제어할 수 있습니다.
- `INPUT_DIR`, `RESULT_DIR`: 데이터 입출력 경로.
- `TH_APBS`: 정전기적 특성 반영 상위 % 임계값.

## 5. 결과물 확인
종료 후 지정된 결과 디렉토리에 다음 파일들이 생성됩니다:
- `Final_Hotspot_Report.txt`: 통합 분석 결과 및 도킹 에너지 요약 리포트.
- `site_0_out.pdb`: 리간드 도킹 최적 포즈 (Consensus 부위 기반).
- `*.dx`: 3D 시각화용 정전기장 맵 파일.
