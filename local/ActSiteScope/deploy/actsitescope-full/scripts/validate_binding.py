import os
import sys
import subprocess
import json
import argparse
import time

class VinaValidator:
    """
    AutoDock Vina를 사용하여 예측된 부위의 결합력을 검증합니다.
    (Meeko 라이브러리를 사용하여 PDBQT 변환을 수행합니다.)
    """
    def __init__(self, project_root=None):
        if project_root is None:
            # 기본적으로 현재 파일(scripts/)의 상단 폴더를 프로젝트 루트로 간주
            project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        self.project_root = project_root
        self.bin_dir = os.path.join(project_root, "scripts", "bin")
        os.makedirs(self.bin_dir, exist_ok=True)
        
        # 환경 변수 지원 (AWS 배포 환경 대응)
        self.vina_path = os.environ.get('VINA_PATH') or os.path.join(self.bin_dir, "vina")
        self.venv_python = os.environ.get('PYTHON_EXECUTABLE') or os.path.join(project_root, "venv", "bin", "python3")
        
        # Docker 환경 등 venv가 없는 경우 대비
        if not os.path.exists(self.venv_python):
            self.venv_python = "python3"
        
    def check_vina(self):
        """Vina 바이너리 존재 확인"""
        return os.path.exists(self.vina_path)

    def prepare_pdbqt(self, receptor_pdb, ligand_sdf, output_dir):
        """
        Meeko를 사용하여 PDBQT 파일을 생성합니다.
        (subprocess를 통해 venv의 meeko 명령어를 실행)
        """
        os.makedirs(output_dir, exist_ok=True)
        receptor_pdbqt = os.path.join(output_dir, "receptor.pdbqt")
        ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")

        print(f"[*] Meeko를 이용한 포맷 변환 시작...")
        
        try:
            # 1. 수용체 변환 (간소화된 방식, 실제로는 전하 보정이 필요할 수 있음)
            # Meeko의 mk_receptor 등을 활용하거나 smina 등의 내장 변환 활용 가능
            # 여기서는 Vina 1.2.x가 지원하는 자동 변환 기능을 염두에 두거나 
            # Meeko를 통한 리간드 변환을 우선적으로 처리합니다.
            
            # 리간드 변환 (Meeko 활용)
            # python -m meeko.prepare_ligand -i ligand.sdf -o ligand.pdbqt
            cmd_ligand = [self.venv_python, "-m", "meeko.prepare_ligand", "-i", ligand_sdf, "-o", ligand_pdbqt]
            subprocess.run(cmd_ligand, check=True)
            
            # 수용체 변환 (Meeko 혹은 단순 PDB->PDBQT)
            # 수용체는 정적인 경우가 많으므로 단순 변환 시도
            # (Vina 1.2.x는 PDB를 직접 받을 수도 있으나 PDBQT가 안정적임)
            return receptor_pdbqt, ligand_pdbqt
        except Exception as e:
            print(f"[!] 변환 실패: {e}")
            return None, None

    def run_docking(self, receptor_pdb, ligand_sdf, center, size=[20, 20, 20], output_dir="results", site_id=0):
        """Vina 1.2.5 실행"""
        os.makedirs(output_dir, exist_ok=True)
        
        # 리간드 고유 이름 추출 (파일명 중복 방지용)
        lig_name = os.path.splitext(os.path.basename(ligand_sdf))[0]

        # PDBQT 변환 (수용체) - Robust converter
        receptor_pdbqt = os.path.join(output_dir, f"site_{site_id}_{lig_name}_receptor.pdbqt")
        cmd_receptor = [self.venv_python, os.path.join(os.path.dirname(__file__), "prepare_vina.py"), "receptor", receptor_pdb, receptor_pdbqt]
        subprocess.run(cmd_receptor, capture_output=True)
            
        # PDBQT 변환 (리간드) - Robust converter
        ligand_pdbqt = os.path.join(output_dir, f"site_{site_id}_{lig_name}_ligand.pdbqt")
        cmd_ligand = [self.venv_python, os.path.join(os.path.dirname(__file__), "prepare_vina.py"), "ligand", ligand_sdf, ligand_pdbqt]
        subprocess.run(cmd_ligand, capture_output=True)

        output_pdbqt = os.path.join(output_dir, f"site_{site_id}_{lig_name}_out.pdbqt")
        log_path = os.path.join(output_dir, f"site_{site_id}_{lig_name}_docking.log")
        
        cmd = [
            self.vina_path,
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
            "--out", output_pdbqt,
            "--exhaustiveness", "8"
        ]
        
        try:
            process = subprocess.run(cmd, capture_output=True, text=True)
            affinity = 0.0
            
            with open(log_path, "w") as f:
                f.write(process.stdout)
                if process.stderr:
                    f.write("\n--- STDERR ---\n")
                    f.write(process.stderr)
                    
            for line in process.stdout.splitlines():
                if line.strip().startswith("1"):
                    parts = line.split()
                    if len(parts) > 1:
                        affinity = float(parts[1])
                    break
            
            # PDBQT -> PDB 변환 (시각화용)
            output_pdb = output_pdbqt.replace(".pdbqt", ".pdb")
            # 단순 변환 (실제로는 정교한 파싱 필요하지만 여기선 파일 존재 확인 위주)
            if os.path.exists(output_pdbqt):
                with open(output_pdbqt, "r") as f_in, open(output_pdb, "w") as f_out:
                    for line in f_in:
                        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("CONECT"):
                            f_out.write(line)

            return {"site_id": site_id, "affinity": affinity, "output_pdb": output_pdb}
        except Exception as e:
            return {"site_id": site_id, "error": str(e)}

    def run_batch_docking_parallel(self, receptor_pdb, ligand_sdf, hotspots, max_workers=4):
        """여러 hotspot에 대해 병렬로 도킹을 수행합니다."""
        from concurrent.futures import ThreadPoolExecutor
        
        results = {}
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for i, hs in enumerate(hotspots):
                futures.append(executor.submit(
                    self.run_docking, receptor_pdb, ligand_sdf, hs["center"], site_id=i
                ))
            
            for future in futures:
                res = future.result()
                if "site_id" in res:
                    results[res["site_id"]] = res
        return results

if __name__ == "__main__":
    # 스크립트 직접 실행 로직
    pass
