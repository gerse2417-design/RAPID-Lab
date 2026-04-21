# ActSiteScope CLI & Docker User Guide (v2.2)

본 프로젝트는 대시보드 시스템의 핵심 분석 파이프라인을 CLI(Command Line Interface)로 제공하여, 대량의 단백질 구조 데이터를 AWS 클라우드 등 배포 환경에서 효율적으로 처리할 수 있도록 지원합니다.

## 1. 주요 특징
- **통합 분석**: M-CSA, AI-Based Binding(MPBind), P2Rank, fPocket, APBS 데이터의 컨센서스 통합.
- **GPU 가속**: NVIDIA CUDA 12.1 기반 AI 추론(MPBind) 지원.
- **지능형 검색**: M-CSA 조회 시 UniProt ID 외에 **PDB ID** 직접 입력 지원 (PDBe API 연동).

## 2. Docker 실행 및 로드

### 이미지 로드 (최초 1회)
```bash
docker load -i actsitescope-full-v2.2.tar.gz
```

### 실행 예시 (로컬 파일 처리)
```bash
docker run --gpus all \
  -v $(pwd)/inputs:/app/inputs \
  -v $(pwd)/results:/app/results \
  actsitescope-full:v2.2 \
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
  actsitescope-full:v2.2 \
  --pdb s3://my-bucket/inputs/protein.pdb \
  --ligands s3://my-bucket/inputs/cmpd1.sdf
```

## 3. CLI 상세 옵션
| 옵션 | 설명 | 비고 |
| :--- | :--- | :--- |
| `--pdb` | 분석할 단백질 PDB 경로 (로컬 또는 S3) | **필수** |
| `--uniprot` | M-CSA 조회를 위한 UniProt ID 또는 **PDB ID** | 지능형 자동 매핑 |
| `--ligands` | 도킹 검증용 리간드 파일 목록 (Space 구분) | 선택 사항 |
| `--env` | 환경 변수 설정 파일 경로 | 기본: `config/.env.aws` |

## 4. 환경 변수 설정 (`.env.aws`)
`config/.env.aws` 파일을 통해 분석 파라미터 및 도구 경로를 제어할 수 있습니다.

## 5. 결과물 확인
종료 후 지정된 결과 디렉토리에 다음 파일들이 생성됩니다:
- `Final_Hotspot_Report.txt`: 통합 분석 결과 및 도킹 에너지 요약 리포트.
- `site_0_out.pdb`: 리간드 도킹 최적 포즈 (Consensus 부위 기반).
