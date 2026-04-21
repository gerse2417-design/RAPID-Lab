# Bioinformatics Pipeline CLI v1.0 실행 가이드 (AWS/On-premise용)

본 가이드는 컨테이너화된 Bioinformatics Pipeline CLI를 실행하고 결과를 확인하는 방법을 설명합니다.

## 1. 이미지 로드 및 환경 설정
이미지 파일(`.tar.gz`)을 서버에 로드하고 환경 변수 파일을 준비합니다.

```bash
# 이미지 로드
gunzip -c my-pipeline-cli-v1.0.tar.gz | docker load

# 환경 변수 확인 (전달받은 config/.env.aws 파일 사용)
ls config/.env.aws
```

## 2. 필수 구성 요소 (Volume Mount)
분석을 위해 다음 두 개의 경로가 마운트되어야 합니다.
1.  **데이터베이스 (`/app/radar/database`)**: 대용량 BOARDS DB가 위치한 경로 (기본 약 7GB)
2.  **데이터 입출력 (`/data`)**: 입력 FASTA 파일과 결과물이 저장될 경로

## 3. 실행 명령어 (Best Practice)
다음 명령어는 **GPU(RTX 40 시리즈 포함) 호환성 패치**와 **일반 사용자 권한 설정**이 모두 적용된 최적의 실행 방식입니다.

```bash
docker run --rm --gpus all \
  --user $(id -u):$(id -g) \
  --env-file config/.env.aws \
  -v /path/to/your/input_data:/data \
  -v /path/to/your/database:/app/radar/database \
  my-pipeline-cli:v1.0 \
  --input /data/target_genome.fasta \
  --output /data/results \
  --run-alphafold
```

### 주요 옵션 설명
*   `--user $(id -u):$(id -g)`: 컨테이너 내부 생성 파일의 소유권을 현재 호스트 사용자로 설정하여 권한 문제를 방지합니다.
*   `--run-alphafold`: Stage 1(RADAR) 이후 Stage 2(AlphaFold2) 단백질 구조 예측을 실행합니다.

## 4. 결과물 위치
호스트 서버의 마운트된 경로(결과 폴더)에서 다음을 확인할 수 있습니다.
*   `*.tsv`: 항생제 내성 분석 리포트
*   `target_protein.faa`: 핵심 히트 단백질 서열
*   `af2_results/`: AlphaFold2 PDB 결과 및 신뢰도 분석 데이터
