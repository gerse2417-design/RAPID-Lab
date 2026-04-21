"""
amr_dock_backend/pdb_utils.py
──────────────────────────────────────────────────────────────
PDB 파일 전처리 유틸리티.

책임 범위:
  - 비표준 아미노산 → 표준 아미노산 치환 (HETATM → ATOM)
  - PDBFixer를 이용한 누락 원자 보완
  - LightDock ANM 호환을 위한 HETATM 제거 클리닝
  - PDB 파일에서 아미노산 서열(1-letter code) 추출
  - 위 단계를 묶은 전처리 파이프라인 오케스트레이터

의존성:
  - 외부: BioPython (Bio.PDB, Bio.SeqUtils), openmm/pdbfixer (선택)
  - 내부: 없음
"""

import shutil
import tempfile
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# PDBFixer 선택적 의존성 가드 (미설치 시에도 모듈 정상 로드)
try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    PDBFIXER_AVAILABLE = True
except ImportError:
    PDBFIXER_AVAILABLE = False

# 비표준 아미노산 → 표준 아미노산 매핑 테이블
NONSTANDARD_RESIDUE_MAP = {
    "MSE": ("MET", "ATOM"),   # Selenomethionine
    "HYP": ("PRO", "ATOM"),   # Hydroxyproline
    "CSE": ("CYS", "ATOM"),   # Selenocysteine
    "TPO": ("THR", "ATOM"),   # Phosphothreonine
    "SEP": ("SER", "ATOM"),   # Phosphoserine
    "PTR": ("TYR", "ATOM"),   # Phosphotyrosine
    "MLY": ("LYS", "ATOM"),   # N-methyl lysine
    "CME": ("CYS", "ATOM"),   # S,S-(2-hydroxyethyl)thiocysteine
    "OCS": ("CYS", "ATOM"),   # Cysteinesulfonic acid
    "KCX": ("LYS", "ATOM"),   # Lysine NZ-carboxylic acid
}


def get_sequence_from_pdb(pdb_path) -> str:
    """PDB 파일에서 아미노산 서열(1-letter code)을 추출."""
    if not pdb_path or not Path(pdb_path).exists():
        return ""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        full_seq = []
        for model in structure:
            for chain in model:
                chain_seq = ""
                for residue in chain:
                    if residue.get_id()[0] == " ":  # standard amino acids only
                        try:
                            aa = seq1(residue.get_resname())
                            if aa and aa != "X" and aa != "*":
                                chain_seq += aa
                        except Exception:
                            continue
                if chain_seq:
                    full_seq.append(chain_seq)
        return "".join(full_seq)
    except Exception as e:
        return f"Error: {str(e)}"


def substitute_nonstandard_residues(input_path, output_path, residue_map=None) -> dict:
    """
    PDB 파일의 HETATM 레코드 중 알려진 비표준 아미노산을
    표준 아미노산 ATOM 레코드로 치환한다.

    Returns
    -------
    dict: {"substitutions": {원래잔기명: 치환횟수}}
    """
    if residue_map is None:
        residue_map = NONSTANDARD_RESIDUE_MAP

    substitutions = {}
    input_path = Path(input_path)
    if not input_path.exists():
        return {"substitutions": substitutions}

    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            if line[:6] == "HETATM":
                resname = line[17:20].strip()
                if resname in residue_map:
                    canonical, new_rec = residue_map[resname]
                    line = f"{new_rec:<6}" + line[6:17] + f"{canonical:<3}" + line[20:]
                    substitutions[resname] = substitutions.get(resname, 0) + 1
            fout.write(line)

    return {"substitutions": substitutions}


def clean_pdb_for_lightdock(input_path, output_path) -> bool:
    """
    LightDock ANM 계산 시 원자 수 불일치 오류 방지:
    표준 아미노산(ATOM)만 남기고 HETATM(물, 이온 등)을 제거.
    오류 발생 시 RuntimeError를 raise하며, 호출부에서 처리한다.
    """
    input_path = Path(input_path)
    if not input_path.exists():
        return False
    try:
        with open(input_path, "r") as fin, open(output_path, "w") as fout:
            for line in fin:
                if line.startswith("ATOM  "):
                    fout.write(line)
                elif line.startswith("TER") or line.startswith("END"):
                    fout.write(line)
        return True
    except Exception as e:
        raise RuntimeError(f"PDB 클리닝 중 오류: {e}") from e


def run_pdbfixer_pipeline(input_path, output_path,
                          add_missing_atoms=True, remove_heterogens=True) -> dict:
    """
    OpenMM PDBFixer로 PDB 파일의 누락 원자(백본/측쇄)를 보완하고
    잔여 HETATM을 제거한다.

    pdbfixer 미설치 또는 오류 발생 시 success=False를 반환하며
    파이프라인을 중단하지 않는다 (graceful degradation).

    Returns
    -------
    dict: {"success": bool, "missing_atoms": int, "skipped_reason": str or None}
    """
    if not PDBFIXER_AVAILABLE:
        return {
            "success": False,
            "missing_atoms": 0,
            "skipped_reason": "pdbfixer 미설치 — 이 단계를 건너뜁니다.",
        }
    try:
        fixer = PDBFixer(filename=str(input_path))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        n_atoms = sum(len(v) for v in fixer.missingAtoms.values())
        if add_missing_atoms:
            fixer.addMissingAtoms()
        if remove_heterogens:
            fixer.removeHeterogens(keepWater=False)
        with open(output_path, "w") as fout:
            PDBFile.writeFile(fixer.topology, fixer.positions, fout)
        return {"success": True, "missing_atoms": n_atoms, "skipped_reason": None}
    except Exception as exc:
        return {
            "success": False,
            "missing_atoms": 0,
            "skipped_reason": f"pdbfixer 오류: {exc}",
        }


def preprocess_pdb_for_docking(input_path, output_path,
                                substitute_nonstandard=False,
                                run_pdbfixer=False,
                                strict_hetatm_removal=True,
                                logger=None) -> dict:
    """
    LightDock 도킹 전 PDB 전처리 파이프라인 오케스트레이터.

    실행 순서:
      1. 비표준 아미노산 치환  (substitute_nonstandard=True 시)
      2. PDBFixer 누락 원자 보완 (run_pdbfixer=True 시)
      3. 엄격 HETATM 제거       (strict_hetatm_removal=True 시, 기본값)

    Parameters
    ----------
    logger : callable or None
        로그 메시지를 받을 콜백 (예: st.write, print). None이면 무시.

    Returns
    -------
    dict: {"substitutions": dict, "pdbfixer_result": dict or None, "warnings": list}
    """
    summary = {"substitutions": {}, "pdbfixer_result": None, "warnings": []}

    with tempfile.TemporaryDirectory() as tmpdir:
        current = str(input_path)

        # Step 1: 비표준 아미노산 치환
        if substitute_nonstandard:
            step1_out = f"{tmpdir}/step1_subst.pdb"
            r = substitute_nonstandard_residues(current, step1_out)
            summary["substitutions"] = r["substitutions"]
            if logger and r["substitutions"]:
                logger("치환: " + ", ".join(
                    f"{k}→{NONSTANDARD_RESIDUE_MAP[k][0]}({v}건)"
                    for k, v in r["substitutions"].items()
                ))
            current = step1_out

        # Step 2: PDBFixer 누락 원자 보완
        if run_pdbfixer:
            step2_out = f"{tmpdir}/step2_fixed.pdb"
            pf = run_pdbfixer_pipeline(current, step2_out)
            summary["pdbfixer_result"] = pf
            if pf["success"]:
                current = step2_out
                if logger:
                    logger(f"PDBFixer: {pf['missing_atoms']}개 원자 보완 완료")
            else:
                summary["warnings"].append(pf["skipped_reason"])
                if logger:
                    logger(f"WARNING: {pf['skipped_reason']}")

        # Step 3: 엄격 HETATM 제거 (최종 안전망)
        if strict_hetatm_removal:
            step3_out = f"{tmpdir}/step3_clean.pdb"
            clean_pdb_for_lightdock(current, step3_out)
            current = step3_out

        shutil.copy2(current, output_path)

    return summary
