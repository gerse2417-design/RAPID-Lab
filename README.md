# RAPID-Lab
Resistance-Averting Protein Interface Design Lab: 내성 억제 단백질 인터페이스 디자인
**By ResistBreakers**: TA 이충환(팀장), TA 박관철, DA 박수연, AA 설민석, AA 김용기

[![Docker Hub](https://img.shields.io/badge/Docker_Hub-syparkbioinfo%2Frapid--lab-2496ED?logo=docker&logoColor=white)](https://hub.docker.com/r/syparkbioinfo/rapid-lab)
## Docker 이미지 접근

Docker Hub 리포: https://hub.docker.com/r/syparkbioinfo/rapid-lab

```bash
docker pull syparkbioinfo/rapid-lab:<tag>
```

사용 가능한 태그:

| 서비스 | 태그 | 크기 |
|---|---|---|
| Hub (랜딩 페이지) | `hub-v1.0` | 770 MB |
| AMR Monitor | `amr-monitor-v2.0` | 17.4 GB |
| ActSiteScope | `actsite-v2.2` | 31.4 GB |
| InhibiForge | `inhibiforge-v3.1` | 30.6 GB |

예:
```bash
docker pull syparkbioinfo/rapid-lab:inhibiforge-v3.1
```

또는 `docker-compose.yml` 안의 `image:` 필드가 이미 올바른 태그를 가리키므로 
`docker compose pull` 한 번이면 4개 다 받을 수 있음.
