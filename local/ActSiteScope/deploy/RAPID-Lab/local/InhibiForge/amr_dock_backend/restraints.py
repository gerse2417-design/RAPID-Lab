"""
amr_dock_backend/restraints.py
──────────────────────────────────────────────────────────────
LightDock 핫스팟 원자 제한 도킹(restraints) 관련 유틸리티.

책임 범위:
  - 사용자 입력 원자 번호 문자열 파싱 (예: "A425, A426") → (체인, 번호) 쌍 목록
  - PDB 파일 ATOM 레코드 고정폭 파싱 → 원자 정보 목록
  - 원자 번호 매칭을 통한 restraints.list 파일 생성

의존성:
  - 외부: re
  - 내부: 없음
"""

import re
from pathlib import Path


def parse_atom_input(raw_input: str) -> list:
    """
    사용자가 입력한 원자 번호 문자열을 파싱하여 (체인ID, Atom Serial Number) 쌍 목록을 반환한다.

    입력 예시: "A425, A426, A427, B12"
    출력: [('A', 425), ('A', 426), ('A', 427), ('B', 12)]

    Raises
    ------
    ValueError : 인식 가능한 원자 번호가 하나도 없을 때
    """
    pattern = re.compile(r'([A-Za-z])(\d+)')
    matches = pattern.findall(raw_input.upper())
    if not matches:
        raise ValueError(
            f"원자 번호 형식을 인식하지 못했습니다: '{raw_input}'\n"
            "올바른 형식 예시: A425, A426, A427"
        )
    return [(chain, int(serial)) for chain, serial in matches]


def parse_pdb_atoms(pdb_path: Path) -> list:
    """
    PDB 파일의 ATOM 레코드를 고정 폭(fixed-width) 슬라이싱으로 파싱한다.
    BioPython 불필요 — 표준 라이브러리만 사용.

    PDB 고정 폭 포맷 (1-indexed):
        1-6   : 레코드 타입 (ATOM  / HETATM)
        7-11  : Atom Serial Number
        13-16 : Atom Name
        18-20 : Residue Name
        22    : Chain ID
        23-26 : Residue Sequence Number

    Returns
    -------
    list of dict: {"serial": int, "name": str, "resname": str, "chain": str, "resseq": int}
    파싱 실패 줄은 조용히 건너뜀.
    """
    pdb_path = Path(pdb_path)
    atoms = []
    if not pdb_path.exists():
        return atoms

    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if not line.startswith("ATOM  "):
                    continue
                try:
                    serial  = int(line[6:11].strip())
                    name    = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain   = line[21].strip()
                    resseq  = int(line[22:26].strip())
                    atoms.append({
                        "serial":  serial,
                        "name":    name,
                        "resname": resname,
                        "chain":   chain,
                        "resseq":  resseq,
                    })
                except (ValueError, IndexError):
                    continue
    except Exception:
        pass

    return atoms


def generate_restraints_list(
    receptor_pdb: Path,
    raw_input: str,
    output_path: Path,
    role: str = "R",
) -> tuple:
    """
    사용자 입력 원자 번호를 PDB 파일과 대조하여 LightDock restraints.list 파일을 생성한다.

    처리 흐름:
      1. parse_atom_input()  → (체인ID, Atom Serial Number) 쌍 목록
      2. parse_pdb_atoms()   → PDB ATOM 레코드 파싱
      3. (체인, 원자번호) 매칭 → 잔기명·잔기번호 획득
      4. 동일 잔기 중복 제거 후 "R A.LEU.25" 형식으로 저장

    Parameters
    ----------
    receptor_pdb : 수용체(또는 리간드) PDB 파일 경로
    raw_input    : 사용자 입력 원자 번호 문자열 (예: "A425, A426")
    output_path  : 생성할 restraints.list 파일 경로
    role         : "R" (수용체) 또는 "L" (리간드)

    Returns
    -------
    (lines_written: list[str], not_found: list[tuple])
    """
    atom_inputs = parse_atom_input(raw_input)
    all_atoms   = parse_pdb_atoms(receptor_pdb)

    atom_map = {(a["chain"], a["serial"]): (a["resname"], a["resseq"]) for a in all_atoms}

    residue_map = {}
    not_found = []

    for chain, serial in atom_inputs:
        key = (chain, serial)
        if key in atom_map:
            resname, resseq = atom_map[key]
            res_key = (chain, resseq)
            if res_key not in residue_map:
                residue_map[res_key] = resname
        else:
            not_found.append(key)

    lines_written = []
    for (chain, resseq), resname in residue_map.items():
        lines_written.append(f"{role} {chain}.{resname}.{resseq}")

    with open(output_path, "w") as f:
        f.write("\n".join(lines_written) + ("\n" if lines_written else ""))

    return lines_written, not_found
