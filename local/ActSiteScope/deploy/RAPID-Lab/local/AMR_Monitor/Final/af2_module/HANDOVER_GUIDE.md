# 🧬 AlphaFold2 Structure 모듈 인수인계 가이드

본 문서는 `af2_structure` 모듈을 다른 자동화 파이프라인에 연동하고자 하는 작업자를 위한 연동 가이드 및 인수인계 명세입니다.

## 1. 모듈 개요 (Module Overview)
`af2_structure`는 단백질 아미노산 서열(`FASTA` 형식)을 입력받아 **AlphaFold2 (ColabFold)** 알고리즘을 통해 고정밀 3D 구조를 예측하는 독립 모듈입니다. 외부 파이프라인에서 래퍼(Wrapper)나 쉘 스크립트를 통해 손쉽게 호출할 수 있도록 설계되었습니다.

- **위치**: `/home/kcpak/projects/af2_structure/`
- **핵심 엔진**: `ColabFold (Dockerized)`
- **연동 인터페이스**: CLI (Shell Script / Python)

---

## 2. 환경 구축 (Environment Setup)

### 🐍 파이썬 의존성 (Host Side)
호스트 머신의 파이썬 환경 구축을 위해 제공된 `requirements.txt`를 활용하십시오.
```bash
# 가상환경 생성 및 활성화 권장
python3 -m venv af2_env
source af2_env/bin/activate
pip install -r requirements.txt
```

### 🐳 인프라 요구사항 (Infrastructure)
본 모듈은 **GPU 가속**을 위해 NVIDIA Docker 환경이 필수적입니다.
- **Docker Image**: `ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2`
- **GPU 권장 사양**: NVIDIA RTX 3090 / A100 이상 (VRAM 12GB+ 권장)
- **가속**: CUDA 지원 필수 (`--gpus all` 옵션 사용 가능 상태)

---

## 3. 연동 인터페이스 (Interface Definition)

### 📥 입력 (Input)
- **파일 형식**: `.faa` 또는 `.fasta` (표준 아미노산 서열)
- **경로 전달**: 실행 스크립트의 첫 번째 인자로 입력 파일의 절대 경로 또는 상대 경로를 전달합니다.

### 🚀 실행 명령 (Execution)
쉘 스크립트를 통한 표준 호출 방식입니다:
```bash
sh run_pipeline.sh <input_file_path.faa>
```
파이썬 내부에서 직접 호출하려는 경우 `af2_runner.py`의 `run_alphafold_docker` 함수를 임포트하여 활용할 수 있습니다.

### 📤 출력 (Output)
모든 결과는 `af2_results/` 디렉터리에 생성됩니다 (출력 경로는 `--output` 인자로 변경 가능).
- **최종 구조**: `*_unrelaxed_rank_001_*.pdb` (가장 신뢰도가 높은 모델)
- **품질 지표**: `*_predicted_aligned_error_v1.json` (PAE 정보)
- **로그**: `af2_history.log` 에 실행 결과(성공여부/소요시간)가 기록됩니다.

---

## 4. 연동 시 주의사항 (Integration Checklist)

> [!IMPORTANT]
> **캐시 관리**: `af2_cache` 및 `af2_jax_cache` 디렉터리는 모델 가중치와 컴파일 데이터를 저장합니다. 이 디렉터리를 유지하면 두 번째 실행부터 속도가 대폭 향상되므로 자동화 파이프라인 구축 시 이 경로를 고정적으로 마운트하십시오.

> [!WARNING]
> **동시성 제어**: AF2는 대량의 GPU VRAM을 점유합니다. 동시에 여러 개의 인스턴스를 실행할 경우 VRAM Out of Memory(OOM)가 발생할 수 있으므로 큐(Queue) 시스템을 통한 순차 처리를 권장합니다.

---
**작성자**: Antigravity (AI Coding Assistant)  
**마지막 업데이트**: 2026-04-05
