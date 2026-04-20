"""
amr_dock_backend/scripts.py
──────────────────────────────────────────────────────────────
외부 도구용 스크립트 문자열 생성기 (Streamlit 의존 없음).

책임 범위:
  - PyMOL 헤드리스 렌더링용 .pml 스크립트 생성 (GSO Heatmap)
  - ChimeraX 인터페이스 분석용 .cxc 스크립트 생성
  - DNAWorks 입력 파일(DNAWORKS.inp) 텍스트 생성

의존성:
  - 외부: pandas (타입 힌트용)
  - 내부: 없음
"""

import pandas as pd


def build_pymol_heatmap_script(receptor_pdb: str, swarm_df: "pd.DataFrame", output_png: str) -> str:
    """
    PyMOL headless .pml 스크립트 문자열을 반환한다.

    - receptor.pdb: surface 모드, 투명도 20%
    - 각 swarm 좌표에 pseudoatom (b-factor = Scoring)
    - spectrum b, blue_red 로 컬러링
    - ray 렌더 후 PNG 저장

    Parameters
    ----------
    receptor_pdb : 수용체 PDB 파일 절대 경로 문자열
    swarm_df     : parse_swarm_coordinates() 반환값 (Swarm, x, y, z, Scoring 컬럼)
    output_png   : 출력 PNG 절대 경로 문자열
    """
    lines = [
        f"load {receptor_pdb}, receptor",
        "show surface, receptor",
        "set transparency, 0.20, receptor",
        "color grey80, receptor",
        "",
        "# Swarm center pseudoatoms",
    ]
    for _, row in swarm_df.iterrows():
        name = f"swarm_{int(row['Swarm'])}"
        lines.append(
            f"pseudoatom {name}, pos=[{row['x']:.3f},{row['y']:.3f},{row['z']:.3f}]"
        )
        lines.append(f"alter {name}, b={row['Scoring']:.5f}")

    lines += [
        "",
        "select all_swarms, pseudoatom",
        "show spheres, all_swarms",
        "set sphere_scale, 1.5, all_swarms",
        "spectrum b, blue_red, all_swarms",
        "",
        "bg_color white",
        "orient",
        "ray 1920, 1080",
        f"png {output_png}, dpi=150",
        "quit",
    ]
    return "\n".join(lines)


def build_chimerax_script(
    complex_pdb: str,
    interface_res: dict,
    mpbind_residues: list,
    dx_file: str = None,
    output_dir: str = ".",
) -> str:
    """
    ChimeraX .cxc 스크립트 문자열을 반환한다.

    - 복합체 cartoon 로드
    - 인터페이스 잔기 stick 모드 (거리 < 4.0Å)
    - hbonds 표시
    - MPBind 잔기 gold highlight
    - dx_file 있으면 coulombic surface
    - 360도 회전 MP4 저장
    - 2K PNG 저장

    Parameters
    ----------
    complex_pdb      : 복합체 PDB 파일 경로
    interface_res    : find_interface_residues() 반환값
    mpbind_residues  : [(chain, resid), ...] 형태의 MPBind 잔기 목록
    dx_file          : APBS .dx 파일 경로 (선택)
    output_dir       : 출력 파일 저장 디렉토리
    """
    lines = [
        f"open {complex_pdb}",
        "preset 'overall look' 'publication 2 (depth-cued)'",
        "",
    ]

    if dx_file:
        lines += [
            f"open {dx_file}",
            "volume #2 style surface",
            "coulombic #1 surfaces #2",
            "",
        ]

    # 인터페이스 잔기 stick
    chains = interface_res.get("chains", [])
    pairs = interface_res.get("pairs", [])
    if pairs and len(chains) >= 2:
        ca, cb = chains[0], chains[1]
        res_a_ids = ",".join(str(p[1]) for p in pairs)
        res_b_ids = ",".join(str(p[4]) for p in pairs)
        lines += [
            "# 인터페이스 잔기 (거리 < 4.0Å)",
            f"style /{ca}:{res_a_ids} stick",
            f"style /{cb}:{res_b_ids} stick",
            "hbonds",
            "",
        ]

    # MPBind 잔기 gold highlight
    if mpbind_residues:
        for chain, resid in mpbind_residues:
            lines.append(f"color /{chain}:{resid} gold")
        lines.append("")

    # 2K 이미지 및 360도 회전 MP4 저장
    png_path = f"{output_dir}/chimerax_interface.png"
    mp4_path = f"{output_dir}/chimerax_360.mp4"
    lines += [
        "view",
        f"save {png_path} width 2048 height 2048 supersample 3",
        "",
        "# 360도 회전 MP4",
        "movie record",
        "turn y 2 180",
        f"movie encode {mp4_path} quality high",
        "",
        "exit",
    ]
    return "\n".join(lines)


def build_dnaworks_input(params: dict) -> str:
    """
    파라미터 dict에서 DNAWORKS.inp 전체 텍스트를 생성하여 반환한다.

    Parameters
    ----------
    params : dict 키 목록
      - title        : str   작업 제목
      - logfile      : str   로그 파일 이름 (기본: "LOGFILE.txt")
      - melting_temp : int   목표 융해 온도 (°C)
      - oligo_length : int   올리고 길이 (nt)
      - codon_org    : str   코돈 역번역 생물종 (예: "ecoli2")
      - sodium_conc  : float 나트륨 이온 농도 (M)
      - magnesium_conc: float 마그네슘 이온 농도 (M)
      - repeat_limit : int   반복 패턴 서열 제한 (nt)
      - n_solutions  : int   설계 후보 수
      - misprime_check: bool 오결합 방어 적용 여부
      - hairpin_check : bool 헤어핀 발생 검사 적용 여부
      - protein_seq  : str   아미노산 서열 (1-letter code)
    """
    title        = params.get("title", "DNAWorks Job")
    logfile      = params.get("logfile", "LOGFILE.txt")
    melting_temp = params.get("melting_temp", 60)
    oligo_length = params.get("oligo_length", 50)
    codon_org    = params.get("codon_org", "ecoli2")
    sodium_conc  = params.get("sodium_conc", 0.05)
    magnesium_conc = params.get("magnesium_conc", 0.002)
    repeat_limit = params.get("repeat_limit", 8)
    n_solutions  = params.get("n_solutions", 5)
    misprime     = params.get("misprime_check", True)
    hairpin      = params.get("hairpin_check", True)
    protein_seq  = params.get("protein_seq", "")

    lines = [
        f'title "{title}"',
        f'logfile "{logfile}"',
        f'melting low {melting_temp}',
        f'length low {oligo_length}',
        f'codon {codon_org}',
        f'concentration sodium {sodium_conc} magnesium {magnesium_conc}',
        f'repeat {repeat_limit}',
        f'solutions {n_solutions}',
    ]
    if hairpin:
        lines.append('hairpin')
    if misprime:
        lines.append('misprime 18 tip 6 max 8')
    lines += [
        'protein',
        f'  {protein_seq.strip()}',
        '//',
    ]
    return "\n".join(lines) + "\n"
