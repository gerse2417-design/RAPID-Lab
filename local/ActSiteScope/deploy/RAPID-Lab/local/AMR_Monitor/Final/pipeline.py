import os
import sys
import json
import numpy as np
import glob
import subprocess
import re
import dateutil.tz
import shutil
from datetime import datetime
import requests
import traceback
from database import init_db, PipelineRunDB, TargetStructureDB, HotspotCandidateDB, BinderDesignDB, DockingEvaluationDB, DNAWorksDB, RadarBoardsDB

# --- 몽키 패치 로직 (DIAMOND) ---
def diamond_make_db(strain, db):
    db_name = db.replace(".faa", "")
    faa = f"radar/database/{db_name}.faa"
    dmnd = f"radar/database/{db_name}"
    print(f"🚀 [Engine] DIAMOND DB 생성: {db_name}")
    try:
        subprocess.run(["diamond", "makedb", "--in", faa, "-d", dmnd], check=True)
    except FileNotFoundError:
        print("❗ 에러: diamond 프로그램이 시스템에 설치되어 있지 않습니다!")
        raise

def diamond_blastp_run(strain, db):
    db_name = db.replace(".faa", "")
    dmnd = f"radar/database/{db_name}.dmnd"
    query = f"radar/pipeline/genome/annotation/2.anno/{strain}/test_genome.faa"
    output = f"radar/pipeline/genome/blastp/3.align/{strain}/test_genome.b6"
    os.makedirs(os.path.dirname(output), exist_ok=True)
    print(f"🚀 [Engine] DIAMOND 정렬 실행: {strain}")
    subprocess.run([
        "diamond", "blastp", "-q", query, "-d", dmnd, "-o", output,
        "--outfmt", "6", "--evalue", "1e-9", "--max-target-seqs", "1"
    ], check=True)

KNOWN_SEQUENCE = "LFTEKFCADGICFIMRAKNEIDHIFSELYSVPNCLQKPYFKLKVQELLLFLCMPLVICTPILIGFAILIPYLCFKNLEKRSIVNRLRAEQKENQQKQVVLALLIHSELFDSGFR"
BOARDS_BASE_URL = "https://sbml.unist.ac.kr/psp"

DEFAULT_CONFIG = {
    "alphafold": {
        "model_type": "auto",
        "num_models": 5,
        "num_recycle": 3,
        "recycle_early_stop_tolerance": 0.5,
        "num_ensemble": 1,
        "num_seeds": 1,
        "random_seed": 0,
        "use_dropout": False,
        "num_relax": 0
    },
    "radar": {"db_name": "BOARDS", "cutoff": 0.95, "evalue": "1e-9", "sample_name": "Enterococcus faecalis V583", "select_by": "gene", "target_list": "1980"},
}

class FinalWorkflowPipeline:
    def __init__(self, db_session, tool_paths=None, custom_config=None):
        self.session = db_session
        self.execution_log = [] # 작업 실행 로그 추적
        self.tool_paths = tool_paths or {
            'colabfold': 'colabfold_batch',
        }
        
        self.config = DEFAULT_CONFIG.copy()
        if custom_config:
            for key, value in custom_config.items():
                if key in self.config and isinstance(value, dict):
                    self.config[key].update(value)

    def _log_task(self, name, start_time, end_time, status="Success"):
        """실행 로그에 작업 정보를 추가합니다."""
        elapsed = end_time - start_time
        self.execution_log.append({
            "task": name,
            "start": start_time.strftime('%Y-%m-%d %H:%M:%S'),
            "end": end_time.strftime('%Y-%m-%d %H:%M:%S'),
            "elapsed": str(elapsed).split('.')[0], # 마이크로초 제외
            "status": status
        })

    def save_run(self, run_id, mode, fasta_path):
        run_db = PipelineRunDB(run_id=run_id, mode=mode, target_fasta_path=fasta_path)
        # Use merge instead of add to handle existing run_id and update it if necessary
        self.session.merge(run_db)
        self.session.commit()
        return run_db

    def run_command(self, cmd_list, log_prefix="", env=None, skip_on_dev=False):
        """명령어를 실행하고 표준 출력을 반환합니다. 환경 변수와 에러 캐칭을 지원합니다."""
        cmd_str = " ".join([str(c) for c in cmd_list])
        print(f"\n======================================")
        print(f"[{log_prefix}] 실제 커맨드 실행 요청:\n{cmd_str}")
        print(f"======================================")
        
        if skip_on_dev:
            print(f"[{log_prefix}] (MOCK) 로컬 개발 모드로 커맨드 건너뜀 (실제 실행시 해제)")
            return True, ""
            
        start_time = datetime.now()
        print(f"[{log_prefix}] 시작 시간: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

        try:
            # 실제 Linux 환경 구동용 (오래 걸릴 수 있음)
            res = subprocess.run(cmd_list, env=env, check=True, text=True, capture_output=True)
            
            end_time = datetime.now()
            self._log_task(log_prefix, start_time, end_time, "Success")
            elapsed = end_time - start_time
            print(f"[{log_prefix}] 종료 시간: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"[{log_prefix}] 경과 시간: {elapsed}")
            print(f"[{log_prefix}] 정상 실행 완료.")
            return True, res.stdout
        except subprocess.CalledProcessError as e:
            end_time = datetime.now()
            self._log_task(log_prefix, start_time, end_time, "Error")
            elapsed = end_time - start_time
            print(f"[{log_prefix}] 종료 시간 (오류): {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"[{log_prefix}] 경과 시간: {elapsed}")
            print(f"[{log_prefix}] 실행 오류: {e.stderr}")
            print(f"[{log_prefix}] 테스트 구동 목적으로 에러를 패스합니다. (파일 미생성 시 파이프라인 정지 가능)")
            return False, e.stdout

    def get_plddt_stats_fixed(self, pdb_path):
        plddts = []
        if not pdb_path or not os.path.exists(pdb_path): return []
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try:
                        plddts.append(float(line[60:66].strip()))
                    except: continue
        return plddts

    def load_radar_results(self, sample_name, strain_tag, output_dir="results"):
        import pandas as pd
        import numpy as np
        
        species = sample_name.split()[0].lower() if sample_name else "unknown"
        tag = strain_tag.strip().lower()

        # Path to pre-calculated or authentic generated reports
        base_dir = "radar/pipeline/result/wgs/table"
        species_dir = os.path.join(base_dir, species)
        
        # Try to locate files
        prefix = f"{species}_{tag}" if tag else species
        stat_path = os.path.join(species_dir, f"{prefix}_statistical_report.tsv")
        clin_path = os.path.join(species_dir, f"{prefix}_clinical_report.tsv")

        if not os.path.exists(clin_path):
            return False, None
            
        df_stat = pd.read_csv(stat_path, sep='\t') if os.path.exists(stat_path) else None
        df_clin = pd.read_csv(clin_path, sep='\t')
        
        if df_clin.empty:
            return False, None

        # Pick the top target
        top_row = df_clin.iloc[0]
        target_name = str(top_row["temp_id"])
        target_id_num = target_name.split("_")[0]

        # Extract Models for the top target
        struct_dir = "radar/database/structures"
        models_data = []
        for i in range(5):
            pdb_path = os.path.join(struct_dir, f"{target_id_num}_ref_rank{i}.pdb")
            if os.path.exists(pdb_path):
                plddt_list = self.get_plddt_stats_fixed(pdb_path)
                avg_plddt = np.mean(plddt_list) if plddt_list else 0.0
                
                with open(pdb_path, "r") as f:
                    pdb_content = f.read()

                models_data.append({
                    "rank": i,
                    "file_name": os.path.basename(pdb_path),
                    "pdb_path": os.path.abspath(pdb_path),
                    "plddt": avg_plddt,
                    "model_name": f"ranked_{i}",
                    "plddt_list": plddt_list,
                    "pdb_data": pdb_content
                })
        
        pae_path = os.path.join(struct_dir, f"{target_id_num}_PAE.png")
        if not os.path.exists(pae_path):
            pae_path = None

        radar_data = {
            "is_hit": True,
            "target_id": target_id_num,
            "target_name": target_name,
            "models": models_data,
            "pae_image": pae_path,
            "df_table1": df_stat,
            "df_table2": df_clin
        }
        
        return True, radar_data

    def run_actual_radar_engine(self, fasta_path):
        import pandas as pd
        engine_start = datetime.now()
        print(f"[BOARDS DB 연동] 시작 시간: {engine_start.strftime('%Y-%m-%d %H:%M:%S')}")
        
        # --- 경로 설정 ---
        RADAR_PATH = os.path.abspath("radar")
        BIN_DIR    = os.path.abspath(".")   # diamond, prodigal 바이너리가 여기 있음
        DIAMOND    = os.path.join(BIN_DIR, "diamond")
        PRODIGAL   = os.path.join(BIN_DIR, "prodigal")
        BOARDS_DB  = os.path.join(RADAR_PATH, "database", "BOARDS")   # .dmnd 이미 존재
        
        # PATH에 바이너리 디렉토리 추가
        os.environ["PATH"] += os.pathsep + BIN_DIR
        
        # --- 파라미터 파싱 ---
        rc = self.config.get('radar', {})
        manual_sample = rc.get("sample_name", "").strip()
        manual_tag    = rc.get("strain_tag",   "").strip()
        user_provided_strain = manual_tag  # 사용자가 실제로 입력한 값 보존 (자동채우기 전)
        cutoff        = float(rc.get("cutoff", 0.95))
        
        # 파일 읽기
        try:
            with open(fasta_path, "r") as f:
                raw_content = f.read()
        except Exception:
            raw_content = ""
        
        # NCBI 헤더에서 자동 파싱 (빈칸일 때만)
        lines = raw_content.splitlines()
        
        # DNA vs Protein 자동 판별 (헤더 파싱보다 먼저 수행)
        clean_seq = "".join(l for l in lines if not l.startswith(">")).upper().replace("\n","").replace(" ","")
        is_protein = any(c in clean_seq for c in "DEFHIKLMPQRSVWY")

        if lines and lines[0].startswith(">"):
            header_full = lines[0][1:].strip()
            
            if not manual_sample:
                # Protein의 경우 A|B|Species 형태가 많으므로 마지막 섹션을 종으로 취급
                if is_protein and "|" in header_full:
                    manual_sample = header_full.split("|")[-1].strip()
                    # 자동으로 tag를 채우지 않음 - 사용자가 입력한 값만 사용
                else:
                    # DNA 등 일반적인 파싱
                    parts = header_full.split()
                    if len(parts) >= 3:
                        if len(parts) >= 4:
                            raw_strain = parts[4].replace(",","") if parts[3].lower()=="strain" and len(parts)>4 else parts[3].replace(",","")
                            manual_sample = f"{parts[1]} {parts[2]} {raw_strain}"
                            manual_tag    = raw_strain.lower()
                        else:
                            manual_sample = f"{parts[1]} {parts[2]}"
            
        final_sample = manual_sample or "Unknown Sample"
        final_tag    = manual_tag    or "user"
        SPECIES      = final_sample.split()[0].lower() if final_sample else "protein"
        
        # --- 헤더 없을 때 보정 (Text Input 지원) ---
        if not lines or not lines[0].startswith(">"):
            raw_content = f">User_Input_{SPECIES}\n" + raw_content
            lines = raw_content.splitlines()

        self.config['radar']['sample_name'] = final_sample
        self.config['radar']['strain_tag']  = final_tag
        
        print(f"📌 Sample Name : {final_sample}")
        print(f"📌 Strain Tag  : {final_tag}")
        print(f"📌 Cutoff      : {cutoff}")
        
        # --- 디렉토리 구성 ---
        fna_dir  = os.path.join(RADAR_PATH, "pipeline", "genome", "sequence", "1.fna",  SPECIES)
        anno_dir = os.path.join(RADAR_PATH, "pipeline", "genome", "annotation", "2.anno", SPECIES)
        b6_dir   = os.path.join(RADAR_PATH, "pipeline", "genome", "blastp", "3.align",  SPECIES)
        report_dir = os.path.join(RADAR_PATH, "pipeline", "result", "wgs", "table", SPECIES)
        
        for d in [fna_dir, anno_dir, b6_dir, report_dir]:
            os.makedirs(d, exist_ok=True)
        
        # --- 서열별 타입 매핑 사전 생성 (멀티 파스타/보정 헤더 지원) ---
        contig_type_map = {}
        curr_id = None
        for line in lines:
            if line.startswith(">"):
                curr_id = line[1:].strip().split()[0]
                is_p = "plasmid" in line.lower()
                if is_protein:
                    contig_type_map[curr_id] = "Protein"
                else:
                    contig_type_map[curr_id] = "Plasmid" if is_p else "Chromosome"

        if is_protein:
            faa_path_dest = os.path.join(anno_dir, "test_genome.faa")
            with open(faa_path_dest, "w") as f:
                f.write(raw_content)
            query_faa = faa_path_dest
            print(f"🧪 아미노산 서열 감지됨. (총 {len(contig_type_map)}개 서열)")
        else:
            fna_path_dest = os.path.join(fna_dir, "test_genome.fna")
            with open(fna_path_dest, "w") as f:
                f.write(raw_content)
            query_faa = os.path.join(anno_dir, "test_genome.faa")
            print(f"🧬 DNA 염기서열 감지됨. (총 {len(contig_type_map)}개 서열, Prodigal 실행...)")
            
            # Prodigal 실행 (DNA → 단백질)
            prodigal_start = datetime.now()
            print(f"[DNA → Protein] Prodigal 시작 시간: {prodigal_start.strftime('%Y-%m-%d %H:%M:%S')}")
            prodigal_cmd = [PRODIGAL, "-i", fna_path_dest, "-a", query_faa, "-p", "meta", "-q"]
            result = subprocess.run(prodigal_cmd, capture_output=True, text=True)
            prodigal_end = datetime.now()
            self._log_task("Prodigal", prodigal_start, prodigal_end)
            if result.returncode != 0:
                raise RuntimeError(f"Prodigal 실패:\n{result.stderr}")
            print(f"✅ Prodigal 완료. (종료: {prodigal_end.strftime('%H:%M:%S')}, 경과: {prodigal_end - prodigal_start})")
        
        # --- DIAMOND BlastP 실행 ---
        b6_path = os.path.join(b6_dir, "test_genome.b6")
        diamond_start = datetime.now()
        print(f"[DIAMOND] BlastP 정렬 시작 시간: {diamond_start.strftime('%Y-%m-%d %H:%M:%S')}")
        diamond_cmd = [
            DIAMOND, "blastp",
            "-q", query_faa,
            "-d", BOARDS_DB,
            "-o", b6_path,
            "--outfmt", "6",
            "--evalue", "1e-9",
            "--max-target-seqs", "1",
            "-p", "4",
            "--quiet"
        ]
        result = subprocess.run(diamond_cmd, capture_output=True, text=True)
        diamond_end = datetime.now()
        self._log_task("BLASTp", diamond_start, diamond_end)
        if result.returncode != 0:
            raise RuntimeError(f"DIAMOND 실패:\n{result.stderr}")
        print(f"✅ DIAMOND 완료. 결과: {b6_path}")
        print(f"[DIAMOND] 종료 시간: {diamond_end.strftime('%Y-%m-%d %H:%M:%S')}, 경과 시간: {diamond_end - diamond_start}")
        
        # --- Cell 2: .b6 → TSV 리포트 파싱 ---
        if not os.path.exists(b6_path) or os.path.getsize(b6_path) == 0:
            raise RuntimeError("DIAMOND 결과(.b6) 파일이 비어있거나 생성되지 않았습니다. 서열이 BOARDS DB와 전혀 매칭되지 않았습니다.")
        
        df = pd.read_csv(b6_path, sep='\t',
                         names=['q_id','s_id','identity','aln_len','mismatch','gap',
                                'q_start','q_end','s_start','s_end','evalue','bitscore'])
        df['p.identity']    = df['identity'] / 100.0
        df['contig_id']     = df['q_id'].str.rsplit('_', n=1).str[0]
        
        # 개별 컨티그 타입 매핑
        if is_protein:
            df['seq_type'] = df['q_id'].map(contig_type_map).fillna("Protein")
        else:
            df['seq_type'] = df['contig_id'].map(contig_type_map).fillna("Chromosome")
            
        df['snp_marker']    = "-"
        # 추출: 무조건 첫 두 단어(학명의 이명법 형태 - Genus species)만 남기기
        clean_species = " ".join(str(final_sample).split()[:2]) if final_sample else "Unknown Sample"
        df['sample_species']= clean_species
        df['strain_tag']    = user_provided_strain if user_provided_strain else "-"
        df['temp_id']       = df['s_id'].str.split('|').str[0]
        
        # BOARDS.faa 메타데이터 매핑 (약물 정보)
        drug_info_db = {
            'van':('Glycopeptides','Vancomycin'),   'mec':('Beta-lactams','Methicillin'),
            'mcr':('Polymyxins','Colistin'),         'tem':('Beta-lactams','Ampicillin'),
            'shv':('Beta-lactams','Amoxicillin'),    'kpc':('Beta-lactams','Carbapenem'),
            'ndm':('Beta-lactams','Carbapenem'),     'vim':('Beta-lactams','Carbapenem'),
            'oxa':('Beta-lactams','Oxacillin'),      'erm':('Macrolides','Erythromycin'),
            'tet':('Tetracyclines','Tetracycline'),  'aac':('Aminoglycosides','Gentamicin'),
            'aph':('Aminoglycosides','Kanamycin'),   'dfr':('Trimethoprim','Trimethoprim'),
            'sul':('Sulfonamides','Sulfamethoxazole'),'cat':('Phenicols','Chloramphenicol'),
            'gyr':('Fluoroquinolones','Ciprofloxacin'),'par':('Fluoroquinolones','Levofloxacin')
        }
        boards_faa = os.path.join(RADAR_PATH, "database", "BOARDS.faa")
        gene_info_map = {}
        if os.path.exists(boards_faa):
            with open(boards_faa, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith(">"):
                        parts2 = line.replace(">","").strip().split("|")
                        sid = parts2[0]
                        orig_gene = "Unknown"
                        for p in parts2:
                            if "_[" in p:
                                orig_gene = p.split("_[")[0].split("|")[-1].split(":")[-1]; break
                        tc, tm = "Others", "Others"
                        for key,(cls,major) in drug_info_db.items():
                            if key in orig_gene.lower(): tc, tm = cls, major; break
                        gene_info_map[sid] = {"amr_gene_detail": orig_gene, "target_antibiotics_class": tc, "target_antibiotics_major": tm}
        
        df['amr_gene_detail']           = df['temp_id'].map(lambda x: gene_info_map.get(x,{}).get('amr_gene_detail','Unknown'))
        df['amr_gene_class']            = df['amr_gene_detail'].apply(lambda g: next((k for k in ['van','mec','mcr'] if k in g.lower()), 'others'))
        df['target_antibiotics_class']  = df['temp_id'].map(lambda x: gene_info_map.get(x,{}).get('target_antibiotics_class','Others'))
        df['target_antibiotics_major']  = df['temp_id'].map(lambda x: gene_info_map.get(x,{}).get('target_antibiotics_major','Others'))
        
        file_prefix = f"{SPECIES}_{final_tag}" if final_tag else SPECIES
        
        # 리포트 1: 통계 (q_id 정제 및 노출)
        df_t1 = df[df['bitscore']>=100].drop_duplicates(subset=['q_id','s_id']).sort_values(['bitscore','evalue'],ascending=[False,True])
        if is_protein:
            # 아미노산 서열 입력시: Sample Name의 괄호 안 번호(NCBI Taxonomy ID)를 추출
            import re as _re
            _m = _re.search(r'\((\d+)\)', str(final_sample))
            ncbi_id = _m.group(1) if _m else "-"
            df_t1['q_id'] = ncbi_id
            cols1 = ['sample_species','strain_tag','seq_type','q_id','temp_id','q_start','q_end','snp_marker','p.identity','bitscore','evalue']
        else:
            cols1 = ['sample_species','strain_tag','seq_type','q_id','temp_id','q_start','q_end','snp_marker','p.identity','bitscore','evalue']
        
        df_t1 = df_t1[[c for c in cols1 if c in df_t1.columns]]
        df_t1.to_csv(os.path.join(report_dir, f"{file_prefix}_statistical_report.tsv"), sep='\t', index=False)
        
        # 리포트 2: 임상 (q_id 정제 및 노출)
        df_t2 = df[(df['bitscore']>=100)&(df['p.identity']>=cutoff)].drop_duplicates(subset=['q_id','s_id']).sort_values(['p.identity','bitscore','evalue'],ascending=[False,False,True])
        if is_protein:
            # 아미노산 서열 입력시: Sample Name의 괄호 안 번호(NCBI Taxonomy ID)를 추출
            import re as _re
            _m = _re.search(r'\((\d+)\)', str(final_sample))
            ncbi_id = _m.group(1) if _m else "-"
            df_t2['q_id'] = ncbi_id
            cols2 = ['sample_species','strain_tag','seq_type','q_id','temp_id','snp_marker','amr_gene_class','amr_gene_detail','target_antibiotics_class','target_antibiotics_major','p.identity','bitscore','evalue']
        else:
            cols2 = ['sample_species','strain_tag','seq_type','q_id','temp_id','snp_marker','amr_gene_class','amr_gene_detail','target_antibiotics_class','target_antibiotics_major','p.identity','bitscore','evalue']
        
        df_t2 = df_t2[[c for c in cols2 if c in df_t2.columns]]
        df_t2.to_csv(os.path.join(report_dir, f"{file_prefix}_clinical_report.tsv"), sep='\t', index=False)
        
        print(f"✅ 리포트 생성 완료 → {report_dir}")
        print(f"   통계 리포트: {len(df_t1)}건 | 임상 리포트: {len(df_t2)}건")
        
        # --- [추가] 가장 유력한 타겟 단백질 서열(Top Hit) 하나만 추출 및 '*' 제거 ---
        if not df_t2.empty:
            # 원본 df에서 정렬 후 첫 번째 유효 q_id 추출
            top_original_q_id = df[(df['bitscore']>=100)&(df['p.identity']>=cutoff)].drop_duplicates(subset=['q_id','s_id']).sort_values(['p.identity','bitscore','evalue'],ascending=[False,False,True])['q_id'].iloc[0]
            
            extracted_faa_path = os.path.join(anno_dir, "target_protein.faa")
            with open(query_faa, 'r') as f_in, open(extracted_faa_path, 'w') as f_out:
                write_flag = False
                for line in f_in:
                    if line.startswith(">"):
                        curr_header_id = line[1:].strip().split()[0]
                        if curr_header_id == top_original_q_id:
                            f_out.write(line)
                            write_flag = True
                        else:
                            write_flag = False
                    elif write_flag:
                        f_out.write(line.replace("*", ""))
            
            query_faa = extracted_faa_path
            print(f"🎯 Top Hit 단백질 서열 추출 및 별표(*) 정제 완료: {extracted_faa_path} (ID: {top_original_q_id})")
        else:
            # 매칭 결과가 없으면 전체 서열에서 '*' 기호만 제거
            cleaned_faa_path = os.path.join(anno_dir, "cleaned_test_genome.faa")
            with open(query_faa, 'r') as f_in, open(cleaned_faa_path, 'w') as f_out:
                for line in f_in:
                    if not line.startswith(">"):
                        f_out.write(line.replace("*", ""))
                    else:
                        f_out.write(line)
            query_faa = cleaned_faa_path
            print(f"⚠️ 매칭 결과 없음. 입력된 모든 서열 대상 특수문자(*) 정제 완료: {cleaned_faa_path}")
            
        # BOARDS PDB 구조 파일 확보 (상위 타겟)
        if not df_t2.empty:
            top_id = df_t2['temp_id'].iloc[0]
            struct_dir = os.path.join(RADAR_PATH, "database", "structures")
            os.makedirs(struct_dir, exist_ok=True)
            for rank in range(5):
                local_pdb = os.path.join(struct_dir, f"{top_id}_ref_rank{rank}.pdb")
                if not os.path.exists(local_pdb):
                    url = f"https://sbml.unist.ac.kr/psp/{top_id}_BOARDS/ref_model/ranked_{rank}.pdb"
                    try:
                        import requests as _req
                        r = _req.get(url, timeout=30)
                        if r.status_code == 200:
                            with open(local_pdb, 'wb') as pf: pf.write(r.content)
                            print(f"⬇️  PDB rank{rank} 다운로드 완료")
                    except Exception:
                        pass
            # PAE 이미지
            pae_local = os.path.join(struct_dir, f"{top_id}_PAE.png")
            if not os.path.exists(pae_local):
                try:
                    import requests as _req
                    r = _req.get(f"https://sbml.unist.ac.kr/psp/{top_id}_BOARDS/ref_model/result_PAE.png", timeout=30)
                    if r.status_code == 200:
                        with open(pae_local, 'wb') as pf: pf.write(r.content)
                except Exception:
                    pass
        
        # 최종 화면 데이터 로드
        res_hit, res_data = self.load_radar_results(final_sample, final_tag)
        if res_hit and res_data is not None:
            res_data["query_faa"] = query_faa
            res_data["is_protein"] = is_protein
            
        engine_end = datetime.now()
        self._log_task("BOARDS DB 연동", engine_start, engine_end)
        return res_hit, res_data




    def _check_boards_id_exists(self, index):
        """Checks if a BOARDS index exists by pinging the info.html page."""
        url = f"{BOARDS_BASE_URL}/{index}_BOARDS/{index}_psp_info.html"
        try:
            response = requests.head(url, timeout=5)
            return response.status_code == 200
        except:
            return False

        # [mock code removed: run_radar_and_check_boards completely migrated to load_radar_results]

    def save_radar_boards(self, run_id, radar_output, boards_id):
        rb = RadarBoardsDB(run_id=run_id, radar_output=radar_output, boards_entry_id=boards_id)
        self.session.add(rb)
        self.session.commit()
        return rb

    # ========================== STAGE 1 ==========================
    def run_alphafold(self, run_id, fasta_path, output_dir):
        af_start = datetime.now()
        conf = self.config['alphafold']
        
        cmd = [self.tool_paths['colabfold']]
        
        if conf.get('model_type'): cmd.extend(["--model-type", conf['model_type']])
        if conf.get('num_models'): cmd.extend(["--num-models", str(conf['num_models'])])
        if conf.get('num_recycle'): cmd.extend(["--num-recycle", str(conf['num_recycle'])])
        if conf.get('recycle_early_stop_tolerance'): cmd.extend(["--recycle-early-stop-tolerance", str(conf['recycle_early_stop_tolerance'])])
        if conf.get('num_ensemble'): cmd.extend(["--num-ensemble", str(conf['num_ensemble'])])
        if conf.get('num_seeds'): cmd.extend(["--num-seeds", str(conf['num_seeds'])])
        if conf.get('random_seed') is not None: cmd.extend(["--random-seed", str(conf['random_seed'])])
        if conf.get('use_dropout'): cmd.append("--use-dropout")
        if conf.get('num_relax', 0) > 0: cmd.extend(["--num-relax", str(conf['num_relax'])])

        # New Advanced MSA/Pairing arguments
        if conf.get('msa_mode'): cmd.extend(["--msa-mode", conf['msa_mode']])
        if conf.get('pair_mode'): cmd.extend(["--pair-mode", conf['pair_mode']])
        if conf.get('template_mode') and conf.get('template_mode') != "none": cmd.append("--custom-template-path") # Simplified representation
        if conf.get('relax_max_iterations', 0) > 0: cmd.extend(["--relax-max-iterations", str(conf['relax_max_iterations'])])
        if conf.get('pairing_strategy'): cmd.extend(["--pairing-strategy", conf['pairing_strategy']])
        if conf.get('max_msa') and conf['max_msa'] != "auto": cmd.extend(["--max-msa", conf['max_msa']])

        cmd.extend([fasta_path, output_dir])
        
        self.run_command(cmd, "AlphaFold2_ColabFold_Batch", skip_on_dev=False)
        
        af_out_dir = output_dir
        
        models_data = []
        import json
        import glob

        # 실제 ColabFold 결과 파일 파싱
        for i in range(conf.get('num_models', 5)):
            pdb_pattern = os.path.join(af_out_dir, f"*_ranked_{i}.pdb")
            pdb_files = glob.glob(pdb_pattern)
            
            if not pdb_files:
                print(f"⚠️ [AlphaFold2] Rank {i} PDB 파일을 찾을 수 없습니다. (Pattern: {pdb_pattern})")
                continue
                
            pdb_file = pdb_files[0]
            abs_pdb_path = os.path.abspath(pdb_file)
            
            # pLDDT 추출
            plddt_list = self.get_plddt_stats_fixed(abs_pdb_path)
            avg_plddt = np.mean(plddt_list) if plddt_list else 0.0
            
            # PAE 데이터 로드 (JSON)
            pae_data = []
            json_pattern = os.path.join(af_out_dir, f"*_predicted_aligned_error_v1_{i}.json")
            json_files = glob.glob(json_pattern)
            if json_files:
                try:
                    with open(json_files[0], 'r') as jf:
                        pae_json = json.load(jf)
                        if isinstance(pae_json, list) and len(pae_json) > 0:
                            pae_data = pae_json[0].get('predicted_aligned_error', [])
                        elif isinstance(pae_json, dict):
                            pae_data = pae_json.get('predicted_aligned_error', [])
                except Exception as e:
                    print(f"⚠️ [AlphaFold2] PAE JSON 파싱 실패: {e}")

            with open(abs_pdb_path, "r") as f:
                pdb_content = f.read()

            models_data.append({
                "rank": i,
                "file_name": os.path.basename(pdb_file),
                "pdb_path": abs_pdb_path,
                "plddt": avg_plddt,
                "model_name": f"ranked_{i}",
                "plddt_list": plddt_list,
                "pae_data": pae_data,
                "pdb_data": pdb_content
            })
            
            model_id = f"{run_id}_rank_{i}"
            db_model = TargetStructureDB(model_id=model_id, run_id=run_id, pdb_path=abs_pdb_path, 
                                         plddt_score=avg_plddt, pae_data_path="")
            self.session.merge(db_model)
            
        self.session.commit()
        af_end = datetime.now()
        self._log_task("AlphaFold2 Structure Prediction (Total)", af_start, af_end)
        return models_data, ""
