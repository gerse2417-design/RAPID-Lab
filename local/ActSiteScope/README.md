# ActSiteScope v2.2 사용 가이드

본 저장소에 업로드된 파일들은 ActSiteScope 대시보드 시스템 및 분석 파이프라인의 **핵심 소스 코드**들입니다.
빠른 소스코드 확인 및 버전 관리를 위해, 11GB 이상의 대용량 AI 모델 가중치(`pytorch_model.bin` 등) 및 의존성 라이브러리가 제외된 **경량화된 소스 코드** 위주로 업로드되어 있습니다.

원활한 시스템 사용 및 구동을 위해 다음 두 가지 설치 방법 중 하나를 선택하여 진행해 주시기 바랍니다.

---

## 🚀 방법 1: Full Version (Docker 패키지) 다운로드 및 구동 (적극 권장)

초기 환경 설정이나 대용량 모델 파일 개별 다운로드 없이, 즉시 실행 가능한 가장 안정적이고 추천하는 방식입니다.

### 1. 풀 버전 패키지 다운로드
아래의 Google Drive 링크에 접속하여 `actsitescope_v2.2_gdrive.tar` (약 9.6GB) 파일을 다운로드합니다.
🔗 **[ActSiteScope v2.2 풀 버전 다운로드 (Google Drive)](https://drive.google.com/drive/u/0/folders/1X7TCGVYN_XEULKxXBy7FDzzOCzbRZ3GG)**

*(해당 압축 파일에는 Ubuntu/Docker 환경에서 즉시 실행 가능한 11GB 모델 내장 Docker 이미지 및 관련 스크립트가 모두 포함되어 있습니다.)*

### 2. 패키지 압축 해제 및 이미지 로드
Linux 환경(Ubuntu 권장)에서 다운로드 받은 파일의 압축을 풀고 Docker 이미지를 로드합니다.
```bash
# 1) 압축 해제
tar -xvf actsitescope_v2.2_gdrive.tar
cd gdrive_release/

# 2) Docker 이미지 로드 (최초 1회, 5~10분 소요될 수 있음)
docker load -i actsitescope-full-v2.2-image.tar.gz
```

### 3. 시스템 실행 (UI 대시보드 모드)
미리 준비된 `docker-compose.yml`을 이용하여 즉시 서버를 구동합니다.
```bash
docker-compose up -d
```
이후 웹 브라우저에서 `http://localhost:8501` 에 접속하시면 대시보드를 바로 사용하실 수 있습니다.

---

## 🛠️ 방법 2: 소스 코드를 로컬 환경에서 직접 구동 (개발자용)

본 GitHub 저장소의 코드를 직접 수정하거나 로컬 Python 환경에서 개발 목적으로 구동하실 때 사용하는 방법입니다.

### 1. 저장소 Clone 및 가상환경 설정
```bash
# 본 저장소를 로컬로 복제합니다.
git clone https://github.com/chlee19990109-cloud/RAPID-Lab.git
cd RAPID-Lab/local/ActSiteScope

# Python 3.10 이상의 가상 환경 생성 및 진입
python3 -m venv venv
source venv/bin/activate

# 필수 라이브러리 설치
pip install -r app/requirements.txt
```

### 2. 누락된 AI 모델 가중치(Weights) 추가 다운로드
본 저장소에는 용량 문제로 MPBind의 ProtT5 모델 가중치가 포함되어 있지 않습니다.
위의 [Google Drive 링크](https://drive.google.com/drive/u/0/folders/1X7TCGVYN_XEULKxXBy7FDzzOCzbRZ3GG) 내 모델 가중치 폴더(또는 관련 배포처)에서 `pytorch_model.bin` (약 11GB) 등 필수 모델 파일들을 다운로드 받은 후, 아래의 지정된 경로에 직접 배치하셔야 합니다.

- **모델 배치 경로**: `scripts/MPBind/src/ProtT5/prot_t5_xl_uniref50/pytorch_model.bin`

*(그 외 P2Rank 모델 파일이나 외부 CLI 바이너리(Vina, FPocket, APBS) 등도 OS 환경에 맞게 추가 설정이 필요할 수 있습니다.)*

### 3. 시스템 실행
필수 모델 및 라이브러리 세팅이 모두 끝난 후 아래 명령어로 Streamlit 대시보드를 실행합니다.
```bash
python -m streamlit run app.py
```
