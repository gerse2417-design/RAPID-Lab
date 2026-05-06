# RAPID-Lab
**Resistance-Averting Protein Interface Design Lab**: 내성 억제 단백질 인터페이스 디자인

[![Docker Hub](https://img.shields.io/badge/Docker_Hub-syparkbioinfo%2Frapid--lab-2496ED?logo=docker&logoColor=white)](https://hub.docker.com/r/syparkbioinfo/rapid-lab)

---

## 1. 프로젝트 개요

**RAPID-Lab** (Resistance-Averting Protein Interface Design + Lab) 은 단순한 데이터 분석의 차원을 넘어, 단백질의 3차원 구조적 완성도를 깊이 있게 이해하고, 이를 바탕으로 내성을 무력화할 수 있는 **내성 억제 단백질을 설계(Design)** 하는 전 과정을 지원하는 통합 환경입니다.

> **Our Vision**  
> *"억제 단백질 설계를 통해 항생제의 유효 수명을 회복하고, 인류의 의료 자산을 지켜냅니다."*

### 프로젝트 기간
**2026.02.24 ~ 2026.04.20** (종료)

### 솔루션 제안 배경

<p align="center">
  <img src="[RAPID-Lab/docs/solution-background.PNG](https://github.com/chlee19990109-cloud/RAPID-Lab/blob/main/docs/solution-background.PNG)" alt="RAPID-Lab 솔루션 제안 배경" width="900"/>
</p>

### 솔루션 컨셉

<p align="center">
  <img src="[RAPID-Lab/docs/solution-.png](https://github.com/chlee19990109-cloud/RAPID-Lab/blob/main/docs/solution-concept.PNG)" alt="RAPID-Lab 솔루션 컨셉" width="900"/>
</p>

---

## 2. 팀원 구성

**ResistBreakers**

| 역할 | 이름 |
|---|---|
| TA (팀장) | 이충환 |
| TA | 박관철 |
| DA | 박수연 |
| AA | 설민석 |
| AA | 김용기 |

---

## 3. 시스템 아키텍처

<p align="center">
  <img src="docs/architecture.png" alt="RAPID-Lab 시스템 아키텍처" width="900"/>
</p>
<p align="center"><em>그림. RAPID-Lab 4-서비스 통합 아키텍처</em></p>

---

## 4. 기술 스택
| 분류 | 사용한 기술 |
|---|---|
| **구조 뷰어/검증** | PyMOL · ChimeraX(complex split 자동화) |
| **WGS 기반 AMR 유전자 탐지·구조 예측** | Prodigal · BLAST · AlphaFold2 |
| **활성·결합부위 예측** | MPBind · P2Rank · fPocket · APBS · PDB2PQR · M-CSA REST API |
| **단백질 설계 SOTA 모델** | RFdiffusion · ProteinMPNN | 
| **설계 검증 및 도킹 · DNA 합성** | AF2-Multimer · LightDock · DNAWorks |
| **인프라** | Python · PyTorch · JAX · Docker · Streamlit · FastAPI · Linux |

---

## 5. 주요 서비스 구성 및 기능

| 서비스 | 핵심 기능 |
|---|---|
| **Hub** | 각 서비스로 이동하는 메인 허브 페이지 |
| **AMR Monitor** | WGS 데이터 기반 항생제 내성 유전자 탐지 및 단백질 구조 예측 |
| **ActSiteScope** | 항생제 내성 단백질 구조 기반 활성부위(결합부위) 탐지 및 예측 |
| **InhibiForge** | RFdiffusion 및 ProteinMPNN 기반 맞춤형 내성 유발 단백질 억제제 설계 및 검증 |

---

## 6. 시연 영상

- 🎬 **[RAPID-Lab-AMR Monitor-시연 영상](https://drive.google.com/file/d/1wmS5g5lzZUKFVNpOrpABbef2RDHGTXcn/view?usp=drive_link)**
- 🎬 **[RAPID-Lab-ActSiteScope-시연 영상](https://drive.google.com/file/d/1FLI0JQnl_QaqB92QdPVLufUXSYuQjtjD/view?usp=drive_link)**
- 🎬 **[RAPID-Lab-InhibiForge-시연 영상](https://drive.google.com/file/d/1v3yChhV8af4Hzs_9Toa_Fnvt8fiU-GZ4/view?usp=drive_link)**


---

## 7. 시스템 요구사항

| 항목 | 요구 사양 |
|---|---|
| **OS** | Ubuntu 24.04 LTS |
| **GPU** | NVIDIA RTX 4070 16GB (개발 환경) |
| **CUDA** | 11.8 이상 |
| **Docker** | 24.x 이상 |
| **Docker Compose** | v2.20 이상 |
| **RAM** | 32GB 이상 권장 |
| **디스크** | 100GB 이상 (4개 서비스 이미지 합계 약 80GB + 데이터 영역) |

> RAPID-Lab 솔루션의 모든 서비스는 GPU 사용을 권장하며, NVIDIA Container Toolkit이 설치되어 있어야 합니다.

---

## 8. 설치 방법 (Docker Hub)

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

또는 `docker-compose.yml` 안의 `image:` 필드가 이미 올바른 태그를 가리키므로 `docker compose pull` 한 번이면 4개를 모두 받을 수 있습니다.

---

## 9. 리포지토리 폴더 구조

```
RAPID-Lab/
├── README.md                   # 본 문서
├── docker-compose.yml          # 4개 서비스 통합 오케스트레이션
├── docs/                       # 프로젝트 문서 및 다이어그램
│   ├── architecture.png
│   ├── solution-background.png
│   └── solution-concept.png
└── local/                      # 서비스별 소스 코드
    ├── hub/                    # 메인 허브 랜딩 페이지(추가 예정)
    ├── Final/amr-monitor/            # AMR Monitor 서비스
    ├── actsite-scope/          # ActSiteScope 서비스
    └── inhibiforge/            # InhibiForge 서비스
```

`local/` 하위 각 서비스 폴더에는 해당 서비스의 소스 코드, Dockerfile, 의존성 정의 파일이 포함되어 있으며, `docs/` 폴더에는 시스템 아키텍처 다이어그램을 비롯한 프로젝트 관련 문서가 포함됩니다.

---

## 10. 트러블슈팅 (FAQ)

**Q1. `docker compose up` 실행 시 GPU를 인식하지 못합니다.**  
A. NVIDIA Container Toolkit 설치 여부와 호스트에서 `nvidia-smi` 정상 동작을 먼저 확인하세요. 이후 아래 명령으로 컨테이너에서 GPU가 보이는지 점검합니다.
```bash
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi
```

**Q2. 이미지 pull 도중 디스크 공간 부족 에러가 발생합니다.**  
A. 4개 이미지 합계가 약 80GB이므로 `/var/lib/docker` 가 위치한 파티션에 충분한 여유 공간이 필요합니다. `docker system df` 로 사용량을 확인하고, 필요 시 `docker system prune -a` 로 미사용 리소스를 정리하세요.

**Q3. 포트 충돌(`port is already allocated`) 이 납니다.**  
A. `docker-compose.yml` 의 `ports:` 필드에서 호스트 포트를 다른 값으로 변경하거나, `lsof -i :<포트>` 로 점유 프로세스를 확인 후 종료하세요.

**Q4. RTX 4070 16GB 외의 GPU에서도 동작하나요?**  
A. CUDA 11.8 호환 GPU라면 동작 가능하나, AlphaFold2 / RFdiffusion 사용 서비스(ActSiteScope, InhibiForge)는 VRAM 12GB 이상을 권장합니다.

---

## 11. 참고 문헌 (References)

Barhoon, M., & Mahdiuni, H. (2025). Exploring protein–protein docking tools: Comprehensive insights into traditional and deep-learning approaches. Journal of Chemical Information and Modeling, 65, 6446–6469. https://doi.org/10.1021/acs.jcim.2025

Campitelli, P., Modi, T., & Ozkan, S. B. (2025). Dissecting allosteric mutations for antibiotic resistance by time-dependent linear response theory. Journal of Chemical Theory and Computation. Advance online publication.

Cortina, G. A., & Kasson, P. M. (2018). Predicting allostery and microbial drug resistance with molecular simulations. Current Opinion in Structural Biology, 52, 80–86. https://doi.org/10.1016/j.sbi.2018.09.001

Dauparas, J., Anishchenko, I., Bennett, N., Bai, H., Ragotte, R. J., Milles, L. F., Wicky, B. I. M., Courbet, A., de Haas, R. J., Bethel, N., Leung, P. J. Y., Huddy, T. F., Pellock, S., Tischer, D., Chan, F., Koepnick, B., Nguyen, H., Qin, A., Laber, A., … Baker, D. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. Science, 378(6615), 49–56. https://doi.org/10.1126/science.abl7484

Hoover, D. M., & Lubkowski, J. (2002). DNAWorks: An automated method for designing oligonucleotides for PCR-based gene synthesis. Nucleic Acids Research, 30(10), e43. https://doi.org/10.1093/nar/30.10.e43

Hoover, D. (2012). Using DNAWorks in designing oligonucleotides for PCR-based gene synthesis. In J. Peccoud (Ed.), Gene synthesis: Methods and protocols (Methods in Molecular Biology, Vol. 852, pp. 209–223). Springer. https://doi.org/10.1007/978-1-61779-564-0_16

Jiménez-García, B., Roel-Touris, J., Romero-Durana, M., Vidal, M., Jiménez-González, D., & Fernández-Recio, J. (2018). LightDock: A new multi-scale approach to protein–protein docking. Bioinformatics, 34(1), 49–55. https://doi.org/10.1093/bioinformatics/btx555

Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., Tunyasuvunakool, K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A., Meyer, C., Kohl, S. A. A., Ballard, A. J., Cowie, A., Romera-Paredes, B., Nikolov, S., Jain, R., Adler, J., … Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596, 583–589. https://doi.org/10.1038/s41586-021-03819-2

Jurrus, E., Engel, D., Star, K., Monson, K., Brandi, J., Felberg, L. E., Brookes, D. H., Wilson, L., Chen, J., Liles, K., Chun, M., Li, P., Gohara, D. W., Dolinsky, T., Konecny, R., Koes, D. R., Nielsen, J. E., Head-Gordon, T., Geng, W., … Baker, N. A. (2018). Improvements to the APBS biomolecular solvation software suite. Protein Science, 27, 112–128. https://doi.org/10.1002/pro.3280

Ko, S., Kim, J., Lim, J., Lee, S.-M., Park, J. Y., Woo, J., Scott-Nevros, Z. K., Kim, J. R., Yoon, H., & Kim, D. (2024). Blanket antimicrobial resistance gene database with structural information, BOARDS, provides insights on historical landscape of resistance prevalence and effects of mutations in enzyme structure. mSystems, 9(1), e00943-23. https://doi.org/10.1128/msystems.00943-23

Krivák, R., & Hoksza, D. (2018). P2Rank: Machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. Journal of Cheminformatics, 10, Article 39. https://doi.org/10.1186/s13321-018-0285-8

Lu, S., Hu, L., Lin, H., Judge, A., Rivera, P., Palaniappan, M., Sankaran, B., Wang, J., Prasad, B. V. V., & Palzkill, T. (2022). An active site loop toggles between conformations to control antibiotic hydrolysis and inhibition potency for CTX-M β-lactamase drug-resistance enzymes. Nature Communications, 13, Article 6726. https://doi.org/10.1038/s41467-022-34564-3

Martins, F. G., Santos, H. A., & Sousa, S. F. (2026). A review of current computational tools for peptide–protein docking. Journal of Computational Chemistry. Advance online publication.

Meng, E. C., Goddard, T. D., Pettersen, E. F., Couch, G. S., Pearson, Z. J., Morris, J. H., & Ferrin, T. E. (2023). UCSF ChimeraX: Tools for structure building and analysis. Protein Science, 32(11), e4792. https://doi.org/10.1002/pro.4792

Mirdita, M., Schütze, K., Moriwaki, Y., Heo, L., Ovchinnikov, S., & Steinegger, M. (2022). ColabFold: Making protein folding accessible to all. Nature Methods, 19, 679–682. https://doi.org/10.1038/s41592-022-01488-1

Roel-Touris, J., Bonvin, A. M. J. J., & Jiménez-García, B. (2020). LightDock goes information-driven. Bioinformatics, 36(3), 950–952. https://doi.org/10.1093/bioinformatics/btz642

Rosignoli, S., & Paiardini, A. (2022). Boosting the full potential of PyMOL with structural biology plugins. Biomolecules, 12(12), Article 1764. https://doi.org/10.3390/biom12121764

Wang, Y., Boadu, F., & Cheng, J. (2025). MPBind: A multitask protein binding site predictor using protein language models and equivariant GNNs. Bioinformatics, 41(11), btaf589. https://doi.org/10.1093/bioinformatics/btaf589

Watson, J. L., Juergens, D., Bennett, N. R., Trippe, B. L., Yim, J., Eisenach, H. E., Ahern, W., Borst, A. J., Ragotte, R. J., Milles, L. F., Wicky, B. I. M., Hanikel, N., Pellock, S. J., Courbet, A., Sheffler, W., Wang, J., Venkatesh, P., Sappington, I., Vázquez Torres, S., … Baker, D. (2023). De novo design of protein structure and function with RFdiffusion. Nature, 620, 1089–1100. https://doi.org/10.1038/s41586-023-06415-8

Wu, M.-H., Xie, Z., & Zhi, D. (2025). A Folding-Docking-Affinity framework for protein-ligand binding affinity prediction. Communications Chemistry, 8, Article 61. https://doi.org/10.1038/s42004-025-01506-1

이병진. (2023). 인공지능 기반 단백질 설계 도구. KOSEN Report. 스탠퍼드대학교 유전학과.

---
