# ActSiteScope Studio (UI & CLI Dual Mode)

본 컨테이너는 대량 데이터 처리를 위한 **CLI(배치) 모드**와 편리한 분석을 위한 **UI(웹 대시보드) 모드**를 모두 지원합니다.

## 1. UI 모드 실행 (추천)

Docker Compose를 사용하면 가장 간편하게 웹 스튜디오를 실행할 수 있습니다.

### 실행 방법
```bash
docker-compose up -d
```
*   **접속 주소**: `http://localhost:8501`
*   **특징**: 마우스 클릭만으로 PDB 업로드, AI 분석, 도킹 검증을 단계별로 수행할 수 있습니다.

## 2. CLI 모드 실행 (배치 처리)

자동화 파이프라인이나 서버에서 명령어로 실행할 때 사용합니다.

### 실행 방법
```bash
docker run --rm --gpus all \
  -e RUN_MODE=cli \
  -v $(pwd)/inputs:/app/inputs \
  -v $(pwd)/results:/app/results \
  actsitescope-cli:v1.1 \
  --pdb /app/inputs/protein.pdb --uniprot P62593
```
*   `-e RUN_MODE=cli`: CLI 모드로 전환합니다.
*   나머지 옵션은 기존 `README_CLI.md`와 동일합니다.

## 3. 모드 전환 원리
컨테이너 시작 시 `RUN_MODE` 환경 변수를 확인합니다.
- `RUN_MODE=ui` (기본값): Streamlit 웹 서버 실행 (8501 포트)
- `RUN_MODE=cli`: Python 분석 스크립트 직접 실행 (인자값 필요)

## 4. 결과물 확인
UI 모드에서 분석 완료 후 상단의 **"통합 분석 리포트 다운로드"** 버튼을 누르거나, 볼륨 마운트된 로컬의 `results/` 폴더에서 `Final_Hotspot_Report.txt`를 확인할 수 있습니다.
