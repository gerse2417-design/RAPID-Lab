# ActSiteScope Studio (UI & CLI Dual Mode)

본 컨테이너는 대량 데이터 처리를 위한 **CLI(배치) 모드**와 편리한 분석을 위한 **UI(웹 대시보드) 모드**를 모두 지원합니다. **v2.2** 버전에는 모든 분석 로직 및 UI 버그를 최종적으로 수정한 안정화 버전입니다.

## 0. 이미지 로드 (최초 1회)
전달받은 외장하드의 압축 파일을 Docker 시스템으로 불러옵니다.
```bash
docker load -i actsitescope-full-v2.2.tar.gz
```

## 1. UI 모드 실행 (추천)
Docker Compose를 사용하면 가장 간편하게 웹 스튜디오를 실행할 수 있습니다.
```bash
docker-compose up -d
```
*   **접속 주소**: `http://localhost:8501`
*   **특징**: M-CSA 표준 모티프 기반 검색, EP 수치 정상화, 도킹 검증 탭 최적화 등 최신 UI 로직 반영.

## 2. CLI 모드 실행 (배치 처리)
```bash
docker run --rm --gpus all \
  -e RUN_MODE=cli \
  -v $(pwd)/inputs:/app/inputs \
  -v $(pwd)/results:/app/results \
  actsitescope-full:v2.2 \
  --pdb /app/inputs/protein.pdb --uniprot P62593
```

## 3. 주요 업데이트 내역

### v2.2 (최신 - 안정화 버전)
- **EP 수치 정상화**: M-CSA 표준 모티프 잔기들의 정전기적 전위(EP) 값이 0으로 출력되던 버그 2종 수정.
- **도킹 검증 뷰어 최적화**: 3D 뷰어에서 M-CSA 잔기 레이블이 중복으로 표시되던 현상 제거.
- **통합 엔진 데이터 정합성**: M-CSA 앵커 잔기의 EP 수치가 최종 통합 분석 결과에도 정확히 반영되도록 수정.

### v2.1
- **M-CSA 지능형 검색**: UniProt ID뿐만 아니라 **PDB ID(예: 1m40)** 직접 입력 시에도 자동 매핑 조회 지원.
- **Anchor 라벨**: 통합 분석 결과 테이블에서 핵심 M-CSA 잔기에 빨간색 **(Anchor)** 강조 표시.
- **리간드 정렬**: 다중 리간드 입력 시 사용자가 업로드한 순서대로 결과 유지.

## 4. 결과물 확인
UI 상단의 **"통합 분석 리포트 다운로드"** 버튼 또는 로컬 `results/` 폴더의 리포트 파일 확인.

