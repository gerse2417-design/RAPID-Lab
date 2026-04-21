import os
import sys
import subprocess
import json
import argparse
import glob
import shutil
import uuid

class MPBindWrapper:
    """
    MPBind (Multi-Scale Protein Binding Site Prediction) 실제 실행을 위한 래퍼 클래스입니다.
    """
    def __init__(self, project_root=None):
        if project_root is None:
            # 기본적으로 현재 파일(scripts/)의 부모 폴더를 프로젝트 루트로 간주
            project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        self.project_root = os.path.abspath(project_root)
        
        # 환경 변수 지원 (Docker/AWS 배포 시 중요)
        self.mpbind_dir = os.environ.get('MPBIND_DIR') or os.path.join(self.project_root, "scripts", "MPBind")
        self.inference_script = os.path.join(self.mpbind_dir, "experiment", "inference.py")
        
        # Docker 환경에서는 시스템 python3를 사용하거나 환경 변수로 지정된 경로 사용
        self.venv_python = os.environ.get('PYTHON_EXECUTABLE') or os.path.join(self.project_root, "venv", "bin", "python3")
        if not os.path.exists(self.venv_python):
            # 런타임에 venv가 없으면 시스템 python3 시도
            self.venv_python = "python3"

    def check_environment(self):
        """
        실행에 필요한 모델 가중치 및 도구가 있는지 확인합니다.
        """
        # 1. 2차 구조 분석 도구(pydssp) 확인
        try:
            # venv_python을 사용하여 pydssp 로드 가능 여부 확인
            test_cmd = [self.venv_python, "-c", "import pydssp; import Bio.PDB; print('OK')"]
            res = subprocess.run(test_cmd, capture_output=True, text=True, timeout=5)
            if res.returncode != 0 or "OK" not in res.stdout:
                return False, "필수 파이썬 패키지(pydssp, biopython)가 설치되지 않았습니다."
        except Exception as e:
            return False, f"환경 검사 중 오류: {str(e)}"
            
        # 2. ProtT5 가중치 확인 (5.3GB)
        prot_t5_path = os.path.join(self.mpbind_dir, "src", "ProtT5", "prot_t5_xl_uniref50")
        if not os.path.exists(prot_t5_path):
            return False, f"ProtT5 모델이 {prot_t5_path}에 없습니다."
            
        return True, "준비 완료"

    def predict(self, pdb_path, output_dir, version="2", bind_types=[0, 3]):
        """
        AI 모델을 사용하여 결합 부위를 예측합니다.
        """
        # 0. 절대 경로 변환
        output_dir = os.path.abspath(output_dir)
        pdb_path = os.path.abspath(pdb_path)
        os.makedirs(output_dir, exist_ok=True)
        
        # 1. 고유 임시 입력 폴더 생성 (동시 실행 방지)
        unique_id = str(uuid.uuid4())[:8]
        rel_input_name = f"tmp_in_{unique_id}"
        tmp_input_dir = os.path.join(self.mpbind_dir, rel_input_name)
        os.makedirs(tmp_input_dir, exist_ok=True)
        
        try:
            # PDB 복사 (이름에서 공백 제거하여 안전하게 처리)
            dest_pdb_name = os.path.basename(pdb_path).replace(" ", "_")
            dest_pdb = os.path.join(tmp_input_dir, dest_pdb_name)
            shutil.copy2(pdb_path, dest_pdb)
            
            # 2. 결과 출력 폴더 설정
            rel_output_name = "results"
            real_output_dir = os.path.join(tmp_input_dir, rel_output_name)
            os.makedirs(real_output_dir, exist_ok=True)

            # 명령어 구성
            cmd = [
                self.venv_python, "inference.py",
                "--input", rel_input_name,
                "--output", rel_output_name,
                "--version", version,
                "--binding_type"
            ] + list(map(str, bind_types))

            print(f"[*] MPBind 실행 커맨드: {' '.join(cmd)}")
            
            env = os.environ.copy()
            # bin 폴더도 PATH에 추가 (mkdssp 등을 위해)
            env["PATH"] = f"{os.path.join(self.mpbind_dir, 'src', 'feature_extraction')}:" + env.get("PATH", "")
            env["PYTHONPATH"] = f"{self.mpbind_dir}:{os.path.join(self.mpbind_dir, 'experiment')}:" + env.get("PYTHONPATH", "")
            
            process = subprocess.run(
                cmd, 
                cwd=os.path.join(self.mpbind_dir, "experiment"),
                capture_output=True,
                text=True,
                env=env
            )
            
            if process.returncode != 0:
                error_log = os.path.join(output_dir, "mpbind_error.log")
                with open(error_log, "w") as ef:
                    ef.write(f"COMMAND: {' '.join(cmd)}\n")
                    ef.write(f"STDOUT:\n{process.stdout}\n\nSTDERR:\n{process.stderr}")
                print(f"❌ MPBind 실행 실패 (Exit Code {process.returncode}). 상세 내용은 'mpbind_error.log'를 확인하세요.")
                return None
                
            # 결과물 파일 찾기
            result_pdb_pattern = os.path.join(real_output_dir, "*.pdb")
            predicted_files = glob.glob(result_pdb_pattern)
            
            if not predicted_files:
                return None
                
            # 결과를 결과 디렉토리로 복사 (Streamlit에서 접근 가능하게)
            final_files = []
            for pf in predicted_files:
                target = os.path.join(output_dir, os.path.basename(pf))
                shutil.copy2(pf, target)
                final_files.append(target)

            return final_files
            
        except Exception as e:
            error_log = os.path.join(output_dir, "mpbind_error.log")
            with open(error_log, "a") as ef:
                ef.write(f"\n[EXCEPTION] {str(e)}\n")
            return None
        finally:
            # 임시 폴더 삭제 (성공/실패 무관하게 정리)
            try:
                shutil.rmtree(tmp_input_dir)
            except:
                pass

    def extract_hotspots(self, predicted_pdb_path, distance_threshold=8.0):
        """
        MPBind 결과 PDB에서 고득점 잔기들을 군집화하여 개별 'Hotspot' 중심점을 찾습니다.
        """
        import numpy as np
        
        atoms = []
        with open(predicted_pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        prob = float(line[60:66].strip())
                        if prob >= 0.5:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            res_name = line[17:20].strip()
                            res_num = line[22:26].strip()
                            atoms.append({
                                "coord": np.array([x, y, z]),
                                "prob": prob,
                                "res": f"{res_name} {res_num}"
                            })
                    except:
                        continue
        
        if not atoms:
            return []

        hotspots = []
        atoms = sorted(atoms, key=lambda x: x['prob'], reverse=True)
        
        for atom in atoms:
            assigned = False
            for hs in hotspots:
                dist = np.linalg.norm(atom['coord'] - hs['center'])
                if dist < distance_threshold:
                    hs['atoms'].append(atom)
                    total_prob = sum(a['prob'] for a in hs['atoms'])
                    hs['center'] = sum(a['coord'] * a['prob'] for a in hs['atoms']) / total_prob
                    assigned = True
                    break
            
            if not assigned:
                hotspots.append({
                    "center": atom['coord'],
                    "atoms": [atom],
                    "prob": atom['prob']
                })
        
        formatted_hotspots = []
        for i, hs in enumerate(hotspots[:10]):
            unique_res = sorted(list(set(a['res'] for a in hs['atoms'])))
            avg_prob = sum(a['prob'] for a in hs['atoms']) / len(hs['atoms'])
            formatted_hotspots.append({
                "site": f"Site {i+1}",
                "center": hs['center'].tolist(),
                "prob": avg_prob,
                "key_residues": ", ".join(unique_res[:3]) + (" ..." if len(unique_res) > 3 else ""),
                "all_residues": ", ".join(unique_res),
                "residues": ", ".join(unique_res)
            })
            
        return formatted_hotspots

def get_real_mpbind_predictions(pdb_path, project_root=None):
    """
    Standalone function to run MPBind and extract hotspots.
    """
    wrapper = MPBindWrapper(project_root)
    # Output results to a subdirectory in project root or current dir
    res_dir = os.path.join(wrapper.project_root, "results")
    os.makedirs(res_dir, exist_ok=True)
    
    results = wrapper.predict(pdb_path, res_dir)
    if not results: return []
    hotspots = wrapper.extract_hotspots(results[0])
    return hotspots

def get_residue_center(pdb_str, residue_query):
    """
    Calculate the center of a residue based on PDB text.
    residue_query format: "Chain resNo" or "resNo"
    """
    import numpy as np
    coords = []
    parts = residue_query.split()
    target_chain = parts[0] if len(parts) > 1 else None
    target_res = parts[1] if len(parts) > 1 else parts[0]
    
    for line in pdb_str.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                r_num = line[22:26].strip()
                chain = line[21].strip()
                if r_num == target_res:
                    if target_chain is None or target_chain == chain:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except:
                continue
    if not coords: return None
    return np.mean(coords, axis=0).tolist()
