import os
import subprocess
import logging
import numpy as np
import sys
from Bio.PDB import PDBParser

class DXParser:
    """
    Parser for OpenDX grid files.
    Calculates electrostatic potential at arbitrary points using trilinear interpolation.
    """
    def __init__(self, dx_path):
        self.dx_path = dx_path
        self.grid = None
        self.origin = None
        self.delta = None
        self.counts = None
        self._parse()

    def _parse(self):
        if not os.path.exists(self.dx_path):
            return
        
        with open(self.dx_path, 'r') as f:
            lines = f.readlines()
            
        data = []
        parsing_data = False
        
        for line in lines:
            if line.startswith('object 1'):
                self.counts = [int(x) for x in line.split()[-3:]]
            elif line.startswith('origin'):
                self.origin = [float(x) for x in line.split()[-3:]]
            elif line.startswith('delta'):
                if not self.delta: self.delta = []
                d = [float(x) for x in line.split()[-3:]]
                self.delta.append(max(d, key=abs)) 
            elif line.startswith('object 3'):
                parsing_data = True
                continue
            elif parsing_data:
                if line.startswith('attribute'):
                    break
                data.extend([float(x) for x in line.split()])
        
        if data:
            self.grid = np.array(data).reshape(self.counts)

    def get_potential(self, x, y, z):
        if self.grid is None:
            return 0.0
        
        ix = (x - self.origin[0]) / self.delta[0]
        iy = (y - self.origin[1]) / self.delta[1]
        iz = (z - self.origin[2]) / self.delta[2]
        
        if ix < 0 or ix >= self.counts[0]-1 or \
           iy < 0 or iy >= self.counts[1]-1 or \
           iz < 0 or iz >= self.counts[2]-1:
            return 0.0
        
        x0, y0, z0 = int(ix), int(iy), int(iz)
        x1, y1, z1 = x0 + 1, y0 + 1, z0 + 1
        dx, dy, dz = ix - x0, iy - y0, iz - z0
        
        v000 = self.grid[x0, y0, z0]
        v100 = self.grid[x1, y0, z0]
        v010 = self.grid[x0, y1, z0]
        v001 = self.grid[x0, y0, z1]
        v110 = self.grid[x1, y1, z0]
        v101 = self.grid[x1, y0, z1]
        v011 = self.grid[x0, y1, z1]
        v111 = self.grid[x1, y1, z1]
        
        c00 = v000 * (1 - dx) + v100 * dx
        c01 = v001 * (1 - dx) + v101 * dx
        c10 = v010 * (1 - dx) + v110 * dx
        c11 = v011 * (1 - dx) + v111 * dx
        
        c0 = c00 * (1 - dy) + c10 * dy
        c1 = c01 * (1 - dy) + c11 * dy
        
        return c0 * (1 - dz) + c1 * dz

    def get_spatial_average_potential(self, x, y, z, radius=2.0, sigma=1.0):
        """
        Calculates the weighted average potential within a sphere around (x,y,z).
        Uses Gaussian weighting based on distance from the center.
        """
        if self.grid is None:
            return 0.0
            
        # 1. Determine grid index range for the sphere (bounding box)
        ix_c = (x - self.origin[0]) / self.delta[0]
        iy_c = (y - self.origin[1]) / self.delta[1]
        iz_c = (z - self.origin[2]) / self.delta[2]
        
        idx_radius_x = int(np.ceil(radius / abs(self.delta[0])))
        idx_radius_y = int(np.ceil(radius / abs(self.delta[1])))
        idx_radius_z = int(np.ceil(radius / abs(self.delta[2])))
        
        nx, ny, nz = self.counts
        
        x_min = max(0, int(np.floor(ix_c - idx_radius_x)))
        x_max = min(nx - 1, int(np.ceil(ix_c + idx_radius_x)))
        y_min = max(0, int(np.floor(iy_c - idx_radius_y)))
        y_max = min(ny - 1, int(np.ceil(iy_c + idx_radius_y)))
        z_min = max(0, int(np.floor(iz_c - idx_radius_z)))
        z_max = min(nz - 1, int(np.ceil(iz_c + idx_radius_z)))
        
        total_weight = 0.0
        weighted_sum = 0.0
        count = 0
        
        # 2. Iterate through index range and accumulate points within the radius
        two_sigma_sq = 2.0 * (sigma ** 2)
        
        for i in range(x_min, x_max + 1):
            for j in range(y_min, y_max + 1):
                for k in range(z_min, z_max + 1):
                    # Actual coordinate of this grid point
                    gx = self.origin[0] + i * self.delta[0]
                    gy = self.origin[1] + j * self.delta[1]
                    gz = self.origin[2] + k * self.delta[2]
                    
                    dist_sq = (gx - x)**2 + (gy - y)**2 + (gz - z)**2
                    if dist_sq <= radius**2:
                        dist = np.sqrt(dist_sq)
                        # Gaussian weight
                        weight = np.exp(-dist_sq / two_sigma_sq)
                        
                        weighted_sum += self.grid[i, j, k] * weight
                        total_weight += weight
                        count += 1
                        
        if total_weight > 0:
            return weighted_sum / total_weight
        
        # Fallback to single point if no grid points found (should not happen with R=2.0)
        return self.get_potential(x, y, z)

def get_pdb_dimensions(pdb_path):
    """
    Calculate the bounding box of a PDB file.
    Returns: (center_x, center_y, center_z), (dim_x, dim_y, dim_z)
    """
    coords = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except:
                    continue
    
    if not coords:
        return (0, 0, 0), (40, 40, 40)
    
    coords = np.array(coords)
    min_c = np.min(coords, axis=0)
    max_c = np.max(coords, axis=0)
    center = (min_c + max_c) / 2.0
    dims = max_c - min_c
    return center, dims

def run_apbs_pipeline(pdb_path, output_dir, project_root=None):
    """
    Run PDB2PQR -> APBS with high stability.
    Uses absolute paths for binaries and relative paths for files to avoid path issues.
    """
    if not os.path.exists(pdb_path):
        return None
    
    if project_root is None:
        # 현재 파일(scripts/)의 상단 폴더를 프로젝트 루트로 간주
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        
    base_name = os.path.splitext(os.path.basename(pdb_path))[0]
    # APBS는 작업 디렉토리 내부에서 상대 경로로 파일을 인식할 때 더 안정적입니다.
    pqr_rel_name = f"{base_name}.pqr"
    in_rel_name = f"{base_name}.in"
    dx_rel_name = f"{base_name}.dx"
    
    # Binary paths - 환경 변수 지원 (AWS 배포 환경 대응)
    pdb2pqr_bin = os.environ.get('PDB2PQR_PATH') or os.path.join(project_root, "venv", "bin", "pdb2pqr")
    apbs_bin = os.environ.get('APBS_PATH') or os.path.join(project_root, "bin", "apbs_bin", "bin", "apbs")
    
    # Calculate Dynamic Grid
    center, dims = get_pdb_dimensions(pdb_path)
    max_dim = max(dims)
    cglen = max_dim * 2.0 + 10.0
    fglen = max_dim * 1.5 + 5.0
    
    # 1. Run PDB2PQR (Absolute paths for CLI, output to results folder)
    pqr_full_path = os.path.join(output_dir, pqr_rel_name)
    try:
        logging.info(f"Running PDB2PQR: {pdb2pqr_bin}")
        # PDB2PQR는 절대 경로를 받아도 무리가 없음
        cmd_pqr = [pdb2pqr_bin, "--ff=PARSE", "--whitespace", "--keep-chain", pdb_path, pqr_full_path]
        res = subprocess.run(cmd_pqr, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        with open(os.path.join(output_dir, "apbs_error.log"), "w") as f:
            f.write(f"PDB2PQR FAILED\nCMD: {' '.join(cmd_pqr)}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
        return None

    # 2. Generate APBS input file (Using relative filenames for internal logic)
    apbs_in_content = f"""
read
    mol pqr {pqr_rel_name}
end
elec mg-auto
    dime 97 97 97
    cglen {cglen:.2f} {cglen:.2f} {cglen:.2f}
    fglen {fglen:.2f} {fglen:.2f} {fglen:.2f}
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    chgm spl2
    srfm smol
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    calcforce no
    calcenergy no
    write pot dx {base_name}
end
quit
"""
    in_full_path = os.path.join(output_dir, in_rel_name)
    with open(in_full_path, 'w') as f:
        f.write(apbs_in_content)

    # 3. Run APBS (CWD set to output_dir, using relative paths for input file)
    try:
        logging.info(f"Running APBS: {apbs_bin}")
        env = os.environ.copy()
        # Set LD_LIBRARY_PATH to find included libs
        ld_path = os.path.join(os.path.dirname(apbs_bin), "..", "lib")
        if os.path.exists(ld_path):
            env["LD_LIBRARY_PATH"] = ld_path + ":" + env.get("LD_LIBRARY_PATH", "")

        cmd_apbs = [apbs_bin, in_rel_name]
        res = subprocess.run(cmd_apbs, check=True, capture_output=True, text=True, env=env, cwd=output_dir)
        
        dx_full_path = os.path.join(output_dir, dx_rel_name)
        # Handle APBS adding .0 suffix
        if not os.path.exists(dx_full_path) and os.path.exists(dx_full_path + ".0"):
            os.rename(dx_full_path + ".0", dx_full_path)
            
        return dx_full_path if os.path.exists(dx_full_path) else None
            
    except subprocess.CalledProcessError as e:
        with open(os.path.join(output_dir, "apbs_error.log"), "w") as f:
            f.write(f"APBS FAILED\nCMD: {' '.join(cmd_apbs)}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
        return None

def score_pocket_electrostatics(pocket_residues_str, dx_path, pdb_path):
    """
    Calculates the average electrostatic potential with improved normalization and matching.
    """
    if not dx_path or not os.path.exists(dx_path):
        return 0.5, 0.0
    
    if not pocket_residues_str:
        return 0.5, 0.0

    try:
        parser_dx = DXParser(dx_path)
        residue_list = [r.strip().upper() for r in pocket_residues_str.split(',')]
        
        parser_pdb = PDBParser(QUIET=True)
        structure = parser_pdb.get_structure("protein", pdb_path)
        
        potentials = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_name = residue.get_resname().strip().upper()
                    res_id = residue.get_id()[1]
                    res_key = f"{res_name} {res_id}"
                    
                    if res_key in residue_list or f"{res_id}" in residue_list:
                        for atom in residue:
                            coord = atom.get_coord()
                            pot = parser_dx.get_spatial_average_potential(coord[0], coord[1], coord[2])
                            potentials.append(pot)
        
        if not potentials:
            return 0.5, 0.0
        
        avg_pot = np.mean(np.abs(potentials))
        score = min(avg_pot / 30.0, 1.0)
        
        if score == 0.5: score = 0.51
        return float(score), float(avg_pot)

    except Exception:
        return 0.5, 0.0

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_pdb = os.path.join(os.path.dirname(__file__), "..", "inputs", "clean_1m40.pdb1")
