# ActSiteScope v2.2 GitHub 업로드 및 병합 가이드

본 가이드는 10GB 규모의 ActSiteScope v2.2 전체 배포판을 GitHub Release의 용량 제한(파일당 2GB) 내에서 배포하고 사용하는 방법을 안내합니다.

## 1. 업로드 방법 (개발자용)

1. GitHub Repository의 **Releases** 페이지로 이동합니다.
2. `Create a new release`를 클릭하고 태그를 `v2.2`로 설정합니다.
3. 다음의 분할된 파일들을 Assets 영역에 모두 드래그하여 업로드합니다:
   - `actsitescope-full-v2.2.tar.gz.partaa`
   - `actsitescope-full-v2.2.tar.gz.partab`
   - `actsitescope-full-v2.2.tar.gz.partac`
   - `actsitescope-full-v2.2.tar.gz.partad`
   - `actsitescope-full-v2.2.tar.gz.partae`
4. 모든 파일이 업로드된 것을 확인한 후 `Publish release`를 클릭합니다.

## 2. 다운로드 및 파일 병합 방법 (사용자용)

본 배포판은 대용량 파일이므로 분할되어 업로드되었습니다. 사용을 위해서는 다운로드 후 다시 하나로 합쳐야 합니다.

### (리눅스/WSL 환경)
1. 위 5개의 파일을 모두 동일한 폴더에 다운로드합니다.
2. 터미널에서 다음 명령어를 실행하여 파일들을 하나로 합칩니다:
   ```bash
   cat actsitescope-full-v2.2.tar.gz.part* > actsitescope-full-v2.2.tar.gz
   ```
3. 합쳐진 파일의 압축을 해제합니다:
   ```bash
   tar -xzvf actsitescope-full-v2.2.tar.gz
   ```

### (윈도우 환경)
1. 7-Zip 또는 WinRAR과 같은 압축 관리 프로그램을 설치합니다.
2. `partaa` 파일을 마우스 오른쪽 버튼으로 클릭하고 "여기서 압축 풀기"를 선택하면 자동으로 나머지 파트들을 인식하여 하나의 폴더로 압축을 해제합니다.

---
**주의**: 모든 파트 파일이 완벽하게 다운로드되어야 병합 및 압축 해제가 가능합니다.
