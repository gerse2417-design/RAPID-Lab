import streamlit as st
import re
import os
import sys
import glob
import pandas as pd
import numpy as np
import streamlit.components.v1 as components
from datetime import datetime
import json
import html

# === Docker/Subfolder Path Configuration ===
# app/streamlit_app.py에서 실행되므로 BASE_DIR은 프로젝트 루트(/app)를 가리켜야 함
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.dirname(CURRENT_DIR) # /app root
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

# 환경 변수에서 경로 가져오기 (Dockerfile의 ENV 설정 활용) 혹은 기본값 설정
INPUT_DIR  = os.environ.get("INPUT_DIR", os.path.join(BASE_DIR, "inputs"))
RESULT_DIR = os.environ.get("RESULT_DIR", os.path.join(BASE_DIR, "results"))
MPBIND_DIR = os.environ.get("MPBIND_DIR", os.path.join(BASE_DIR, "scripts", "MPBind"))

os.makedirs(INPUT_DIR,  exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)

# V2 백엔드 유틸리티 및 엔진 임포트
from scripts.mcsa_utils import get_mcsa_residues_by_uniprot, get_mcsa_residues_by_pdb_id, get_mcsa_residues_by_keyword
from scripts.pdb_utils import extract_pdb_info, clean_pdb_file
from scripts.pocket_utils import run_p2rank, run_fpocket, get_pocket_overlap
from scripts.apbs_utils import score_pocket_electrostatics, run_apbs_pipeline, DXParser
from scripts.ranking_utils import rank_pockets, calculate_motif_concordance
from scripts.predict_site import MPBindWrapper
from scripts.consensus_analysis import ConsensusEngine
import importlib
importlib.reload(sys.modules['scripts.consensus_analysis'])
from scripts.consensus_analysis import ConsensusEngine
from scripts.validate_binding import VinaValidator

# === 페이지 설정 및 테마 초기화 ===
st.set_page_config(
    page_title="ActSiteScope",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# === 디자인 변수 (라이트 모드 고정) ===
vars = {
    "main_bg": "#f8fafc",
    "sidebar_bg": "#f1f5f9",
    "text_primary": "#0f172a",
    "text_secondary": "#475569",
    "glass_bg": "rgba(255, 255, 255, 0.8)",
    "glass_border": "rgba(203, 213, 225, 0.8)",
    "card_bg": "#ffffff",
    "header_color": "#0f172a",
    "primary": "#2563eb"
}

st.markdown(f"""
<link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;400;500;600;700&family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<link href="https://fonts.googleapis.com/css2?family=Material+Symbols+Outlined:wght,FILL@100..700,0..1" rel="stylesheet">
<style>
    :root {{
        --main-bg: {vars['main_bg']};
        --sidebar-bg: {vars['sidebar_bg']};
        --text-primary: {vars['text_primary']};
        --text-secondary: {vars['text_secondary']};
        --glass-bg: {vars['glass_bg']};
        --glass-border: {vars['glass_border']};
        --card-bg: {vars['card_bg']};
        --header-color: {vars['header_color']};
        --primary: {vars['primary']};
    }}
    header[data-testid="stHeader"] {{ background: transparent !important; }}
    .main {{ background-color: var(--main-bg) !important; color: var(--text-primary) !important; font-family: 'Inter', sans-serif !important; }}
    .stApp, .stAppViewContainer {{ background-color: var(--main-bg) !important; color: var(--text-primary) !important; }}
    .main .block-container {{ max-width: 95% !important; padding-top: 2rem !important; padding-bottom: 0rem !important; padding-left: 5rem !important; padding-right: 5rem !important; }}
    #stSidebar, section[data-testid="stSidebar"], [data-testid="stSidebarContent"] {{ background-color: var(--sidebar-bg) !important; border-right: 1px solid var(--glass-border) !important; padding-top: 2rem !important; }}
    button[data-testid="stSidebarCollapse"] svg, button[aria-label="Expand sidebar"] svg, button[aria-label="Collapse sidebar"] svg, header[data-testid="stHeader"] button svg {{ color: #000000 !important; fill: #000000 !important; opacity: 1 !important; }}
    .material-symbols-outlined {{ font-variation-settings: 'FILL' 0, 'wght' 300, 'GRAD' 0, 'opsz' 24; vertical-align: middle; }}
    .glass-panel {{ background: var(--glass-bg) !important; backdrop-filter: blur(12px); -webkit-backdrop-filter: blur(12px); border: 1px solid var(--glass-border); }}
    .stButton>button {{ border-radius: 0.375rem !important; border: 1px solid var(--glass-border) !important; background-color: var(--card-bg) !important; color: var(--text-primary) !important; font-weight: 500 !important; font-size: 16px !important; transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1) !important; }}
    .stButton>button:hover {{ background-color: rgba(59, 130, 246, 0.1) !important; border-color: var(--primary) !important; color: var(--primary) !important; }}
    .gold-btn button {{ background: linear-gradient(135deg, rgba(234, 179, 8, 0.1) 0%, rgba(234, 179, 8, 0.05) 100%) !important; border: 1px solid rgba(234, 179, 8, 0.3) !important; color: #eab308 !important; }}
    .gold-btn button:hover {{ background: rgba(234, 179, 8, 0.2) !important; border-color: #eab308 !important; }}
    .final-results-table {{ width: 100%; table-layout: fixed; border-collapse: collapse; font-size: 18px !important; color: var(--text-primary) !important; }}
    .final-results-table th {{ background-color: var(--glass-bg); color: var(--header-color) !important; padding: 12px; text-align: left; font-weight: 700; text-transform: uppercase; border-bottom: 2px solid var(--glass-border); font-size: 20px !important; }}
    .final-results-table td {{ padding: 12px; border-bottom: 1px solid var(--glass-border); color: var(--text-primary) !important; }}
    .fixed-docking-table {{ table-layout: fixed !important; width: 100% !important; border-collapse: collapse; }}
    .fixed-docking-table td {{ word-break: break-all !important; white-space: normal !important; line-height: 1.4 !important; vertical-align: top !important; }}
    [data-testid="stWidgetLabel"], [data-testid="stMarkdownContainer"] p, section[data-testid="stSidebar"] label {{ color: var(--text-primary) !important; opacity: 1 !important; }}
    [data-testid="stFileUploader"] {{ background-color: var(--card-bg) !important; border: 2px dashed var(--glass-border) !important; border-radius: 1rem !important; padding: 10px !important; }}
    [data-testid="stFileUploader"] section {{ background-color: transparent !important; }}
    [data-testid="stFileUploader"] div[data-testid="stMarkdownContainer"] p {{ color: var(--text-primary) !important; opacity: 0.7 !important; }}
    div[data-baseweb="select"] > div {{ background-color: var(--card-bg) !important; color: var(--text-primary) !important; border-color: var(--glass-border) !important; }}
    input {{ color: var(--text-primary) !important; background-color: var(--card-bg) !important; }}
    .stDownloadButton > button {{ color: var(--text-primary) !important; background-color: var(--card-bg) !important; border: 1px solid var(--glass-border) !important; }}
    .stDownloadButton > button:hover {{ border-color: var(--primary) !important; color: var(--primary) !important; }}
    [data-testid="stStatusWidget"], [data-testid="stStatusWidget"] div, [data-testid="stStatusWidget"] p, [data-testid="stStatusWidget"] summary {{ color: var(--text-primary) !important; }}
    [data-testid="stStatusWidget"] {{ background-color: var(--card-bg) !important; border: 1px solid var(--glass-border) !important; }}
    [data-testid="stStatusWidget"] summary {{ background-color: #FFFFFF !important; color: #000000 !important; font-weight: 700 !important; padding: 5px 10px !important; border-radius: 0.5rem !important; }}
    div[data-testid="stToggle"] > div {{ background-color: var(--primary) !important; border: 1px solid #000000 !important; }}
    div[data-testid="stToggle"] {{ border: 1px solid #000000 !important; border-radius: 2rem !important; padding: 2px !important; }}
</style>
""", unsafe_allow_html=True)

st.markdown(f"""
<script src="https://cdn.tailwindcss.com?plugins=forms,container-queries"></script>
<script>
    tailwind.config = {{
      darkMode: "class",
      theme: {{
        extend: {{
          "colors": {{
                  "primary": "{vars['primary']}",
                  "surface-container-lowest": "{vars['main_bg']}",
                  "outline-variant": "{vars['glass_border']}",
                  "on-surface": "{vars['text_primary']}",
                  "secondary": "#4edea3",
                  "surface": "{vars['main_bg']}",
                  "background": "{vars['main_bg']}",
                  "accent-gold": "#eab308"
          }}
        }}
      }}
    }}
</script>
""", unsafe_allow_html=True)

# Top Navigation Bar
st.markdown(f'<header class="fixed top-0 left-0 w-full z-[1000] bg-[{vars["main_bg"]}]/80 backdrop-blur-2xl px-6 h-16 border-b border-black/5 flex justify-between items-center shadow-sm"><div class="flex items-center gap-3"><div class="w-8 h-8 bg-primary/20 rounded-lg flex items-center justify-center border border-primary/30"><span class="material-symbols-outlined text-primary text-xl">science</span></div><h1 class="text-[45px] font-bold tracking-widest text-[{vars["header_color"]}] uppercase opacity-90">ActSiteScope</h1></div></header><div style="height: 64px;"></div>', unsafe_allow_html=True)

# === 세션 상태 초기화 ===
for k in ["mcsa_residues", "mp_hotspots", "p2rank_results", "fpocket_results", "apbs_results",
          "hotspots", "consensus_hotspots", "consensus_formats", "vina_results_mcsa", "vina_results_consensus"]:
    if k not in st.session_state: st.session_state[k] = []
for k in ["site_affinities", "site_ligands"]:
    if k not in st.session_state: st.session_state[k] = {}
for k in ["v_mcsa", "v_mpbind", "v_p2rank", "v_fpocket", "v_apbs"]:
    if k not in st.session_state: st.session_state[k] = False

# 시각화 기본 상태: 초기 업로드 시에는 False여야 함
if "show_stability" not in st.session_state: st.session_state["show_stability"] = False
if "show_residues" not in st.session_state: st.session_state["show_residues"] = False
if "extracted_uniprot_id" not in st.session_state: st.session_state["extracted_uniprot_id"] = None
if "extracted_pdb_id" not in st.session_state: st.session_state["extracted_pdb_id"] = None
if "extracted_keyword" not in st.session_state: st.session_state["extracted_keyword"] = None

# ============================================================
# 헬퍼 함수 정의 (필수 유틸리티)
# ============================================================
def get_system_status():
    status_list = []
    prot_t5 = os.path.join(MPBIND_DIR, "src", "ProtT5", "prot_t5_xl_uniref50", "pytorch_model.bin")
    status_list.append(("MPBind 모델", "✅ 준비 완료" if os.path.exists(prot_t5) else "⚠️ 모델 가중치가 없습니다.", "success" if os.path.exists(prot_t5) else "warning"))
    vina = os.environ.get("VINA_PATH", os.path.join(BASE_DIR, "scripts", "bin", "vina"))
    status_list.append(("Vina 도킹 엔진", "✅ 준비 완료" if os.path.exists(vina) else "❌ 미설치", "success" if os.path.exists(vina) else "error"))
    p2rank = os.environ.get("P2RANK_PATH", os.path.join(BASE_DIR, "scripts", "bin", "p2rank", "prank"))
    jre = os.path.join(BASE_DIR, "scripts", "bin", "jre", "bin", "java")
    status_list.append(("P2Rank (Hotspot)", "✅ 준비 완료" if (os.path.exists(p2rank) and os.path.exists(jre)) else "❌ 미설치", "success" if (os.path.exists(p2rank) and os.path.exists(jre)) else "error"))
    fpocket = os.path.join(BASE_DIR, "scripts", "bin", "fpocket", "fpocket")
    status_list.append(("Fpocket (Cavity)", "✅ 준비 완료" if os.path.exists(fpocket) else "❌ 미설치", "success" if os.path.exists(fpocket) else "error"))
    apbs_bin = os.path.join(BASE_DIR, "bin", "apbs_bin", "bin", "apbs")
    status_list.append(("APBS (Electrostatics)", "✅ 준비 완료" if os.path.exists(apbs_bin) else "❌ 미설치", "success" if os.path.exists(apbs_bin) else "error"))
    return status_list

def get_pdb_residue_map(pdb_str):
    mapping = {}
    if not pdb_str: return mapping
    for line in pdb_str.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21].strip() or 'A'
            res_seq = line[22:26].strip()
            res_name = line[17:20].strip().upper()
            mapping[(chain, res_seq)] = res_name
    return mapping

def generate_styled_table(data, cols, res_map, col_map=None, widths=None):
    if not col_map: col_map = {}
    html = '<div class="overflow-x-auto"><table class="final-results-table"><thead><tr>'
    
    # Calculate widths indexing if provided
    w_idx = 0
    for c in cols:
        header = col_map.get(c, c.upper())
        w_attr = f' style="width: {widths[w_idx]}%;"' if widths and w_idx < len(widths) else ""
        html += f'<th{w_attr}>{header}</th>'
        w_idx += 1
        
        if c == "residues":
            w_attr_res = f' style="width: {widths[w_idx]}%;"' if widths and w_idx < len(widths) else ""
            html += f'<th{w_attr_res}>아미노산 정보</th>'
            w_idx += 1
            
    html += '</tr></thead><tbody>'
    for item in data:
        html += '<tr>'
        for c in cols:
            val = item.get(c, "")
            if c == "residues":
                ident, names = format_residue_info_columns(str(val), res_map)
                html += f'<td>{ident}</td><td>{names}</td>'
            elif isinstance(val, (float, np.float32, np.float64)):
                html += f'<td class="font-mono text-primary">{val:.3f}</td>'
            else:
                html += f'<td>{val}</td>'
        html += '</tr>'
    return html + '</tbody></table></div>'

def format_residue_info_columns(residue_str, res_map):
    if not residue_str: return "", ""
    import re
    # 1. 쉼표와 공백을 모두 처리하기 위해 정규식으로 분리
    raw_items = re.split(r'[,\s]+', str(residue_str).strip())
    
    # 2. 토큰을 (체인, 번호, 원본이름) 형태로 정규화
    targets = []
    idx = 0
    while idx < len(raw_items):
        item = raw_items[idx]
        if not item: 
            idx += 1
            continue
            
        if '_' in item:
            # P2Rank 형식: A_70
            ch_n = item.split('_')
            ch = ch_n[0] if ch_n[0] else 'A'
            n = ch_n[1] if len(ch_n) > 1 else ''
            targets.append((ch, n, None))
            idx += 1
        elif item.isalpha() and len(item) == 3 and idx + 1 < len(raw_items):
            # MPBind 형식: THR 70
            targets.append(('A', raw_items[idx+1], item))
            idx += 2
        else:
            # 기타 (번호만 있거나 알 수 없는 형식)
            targets.append(('A', item, None))
            idx += 1

    ids = []
    names = []
    
    for ch_in, n, orig_name in targets:
        if not n: continue
        
        # 1. 주어진 체인 및 번호로 우선 조회
        r_name = res_map.get((ch_in, n))
        target_id = f"{ch_in}{n}"
        
        # 2. 실패 시 번호만으로 전체 체인 검색
        if r_name is None:
            found = False
            for (m_ch, m_seq), name in res_map.items():
                if m_seq == n:
                    r_name = name
                    target_id = f"{m_ch}{n}"
                    found = True
                    break
            if not found:
                r_name = orig_name if orig_name else "UNK"
                target_id = f"{ch_in}{n}"
        
        ids.append(target_id)
        names.append(f"{r_name} {n}")
        
    return ", ".join(ids), ", ".join(names)

def reset_parameter_defaults():
    rc = st.session_state.get("reset_count", 0)
    st.session_state[f"cfg_mp_{rc}"] = 0.5
    st.session_state[f"cfg_fpc_{rc}"] = 0.5
    st.session_state[f"cfg_p2r_{rc}"] = 30
    st.session_state[f"cfg_apbs_{rc}"] = 30

def render_3d_view(pdb_str, hotspots=None, consensus=None, dx_str=None, ligand_str=None, mcsa_data=None,
                   ligand_name=None, ligand_format="pdb", show_legend=True, mode="individual", show_stability=False,
                   show_base_residues=False, vis_mcsa=False, vis_mpbind=False, vis_p2rank=False,
                   vis_fpocket=False, vis_apbs=False, limit_atoms=None):
    import json

    def _get_top_pose(pdb_str, max_atoms=None):
        if not pdb_str: return pdb_str
        lines = []
        atom_count = 0
        for line in pdb_str.split('\n'):
            if line.startswith("MODEL") and atom_count > 0: break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_count += 1
                if max_atoms and atom_count > max_atoms: break
                try:
                    serial = int(line[6:11].strip())
                    if serial == 1 and atom_count > 1: break
                except: pass
            lines.append(line)
            if atom_count > 0:
                if line.startswith("ENDMDL") or line.startswith("TER") or line.startswith("END"): break
        return "\n".join(lines)

    processed_ligand = _get_top_pose(ligand_str, max_atoms=limit_atoms) if ligand_str and mode == "docking" else (ligand_str if ligand_str else "")

    data_json = json.dumps({
        "pdb": pdb_str if pdb_str else "",
        "hotspots": hotspots if hotspots else [],
        "ligand": processed_ligand,
        "ligand_format": ligand_format if ligand_format else "pdb",
        "mcsa": mcsa_data if mcsa_data else [],
        "consensus": consensus if consensus else [],
        "dx": dx_str if dx_str else "",
        "mode": mode,
        "ligand_name": ligand_name,
        "show_stability": show_stability,
        "show_base_residues": show_base_residues,
        "vis_mcsa": vis_mcsa,
        "vis_mpbind": vis_mpbind,
        "vis_p2rank": vis_p2rank,
        "vis_fpocket": vis_fpocket,
        "vis_apbs": vis_apbs,
        "theme": "Light"
    })

    viewer_bg  = "#f8fafc"
    panel_bg   = "rgba(255, 255, 255, 0.8)"
    panel_bdr  = "rgba(0,0,0,0.1)"
    text_color = "#0f172a"
    shadow     = "rgba(0,0,0,0.1)"
    legend_display = "block" if show_legend else "none"

    html_code = """
    <!DOCTYPE html><html><head><meta charset="utf-8">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.3/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.1/3Dmol-min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pako/2.1.0/pako.min.js"></script>
    <style>
        body, html { margin:0; padding:0; width:100%; height:100%; overflow:hidden; background:VAR_VIEWER_BG; color:VAR_TEXT_COLOR; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }
        #viewer3d { width:100%; height:500px; min-height:500px; position:relative; background:VAR_VIEWER_BG; }
        .glass-panel { background:VAR_PANEL_BG; backdrop-filter:blur(12px); -webkit-backdrop-filter:blur(12px); border:1px solid VAR_PANEL_BDR; }
        #diag-overlay { position:absolute; bottom:20px; left:24px; display:flex; align-items:center; gap:12px; padding:8px 16px; border-radius:10px; z-index:1000; box-shadow:0 10px 30px rgba(0,0,0,0.4); }
        #premium-legend { position:absolute; top:24px; right:24px; width:160px; border-radius:14px; padding:16px; z-index:1000; box-shadow:0 20px 50px VAR_SHADOW; display:VAR_LEGEND_DISPLAY; }
        .pulse { width:8px; height:8px; border-radius:50%; background:#4edea3; animation:pulse-anim 2s infinite; }
        @keyframes pulse-anim { 0%{transform:scale(0.95);box-shadow:0 0 0 0 rgba(78,222,163,0.7);} 70%{transform:scale(1);box-shadow:0 0 0 8px rgba(78,222,163,0);} 100%{transform:scale(0.95);box-shadow:0 0 0 0 rgba(78,222,163,0);} }
        .dot { display:inline-block; width:8px; height:8px; border-radius:50%; margin-right:12px; }
        .header-label { font-size:14px; font-weight:bold; text-transform:uppercase; letter-spacing:0.15em; color:VAR_TEXT_COLOR; margin-bottom:15px; border-bottom:1px solid VAR_PANEL_BDR; padding-bottom:8px; }
        .loading { position:absolute; top:50%; left:50%; transform:translate(-50%,-50%); color:VAR_TEXT_COLOR; opacity:0.4; font-size:12px; letter-spacing:0.2em; text-transform:uppercase; font-weight:700; }
        ul { list-style:none; padding:0; margin:0; }
        li { display:flex; align-items:center; font-size:14px; color:VAR_TEXT_COLOR; margin-bottom:10px; opacity:0.8; }
    </style></head><body>
    <div id="viewer3d">
        <div class="loading">Initializing Engine...</div>
        <div id="diag-overlay" class="glass-panel">
            <span style="font-size:12px;font-family:monospace;font-weight:800;letter-spacing:0.1em;color:VAR_TEXT_COLOR;display:flex;align-items:center;gap:8px;">
                <div class="pulse" style="background:#10b981;"></div>SYSTEM: <span id="sys-mode">READY</span>
            </span>
            <span style="color:rgba(0,0,0,0.15);">|</span>
            <span id="atom-stat" style="font-size:12px;font-family:monospace;color:VAR_TEXT_COLOR;font-weight:600;">0 atoms</span>
            <span style="color:rgba(0,0,0,0.15);">|</span>
            <span id="ligand-stat" style="font-size:12px;font-family:monospace;color:VAR_TEXT_COLOR;font-weight:800;">LIGAND: NONE</span>
        </div>
        <div id="premium-legend" class="glass-panel">
            <div class="header-label">Analysis Legend</div>
            <ul>
                <li><span class="dot" style="background:#000;border:0.5px solid rgba(0,0,0,0.2);"></span> M-CSA Catalytic Motif</li>
                <li><span class="dot" style="background:#00008F;border:0.5px solid rgba(0,0,0,0.2);"></span> AI-Based Binding Site (MPBind)</li>
                <li><span class="dot" style="background:#FF5722;border:0.5px solid rgba(0,0,0,0.2);"></span> Geometry (P2Rank)</li>
                <li><span class="dot" style="background:#2ECC71;border:0.5px solid rgba(0,0,0,0.2);"></span> Druggable (fPocket)</li>
                <li><span class="dot" style="background:#eab308;border:0.5px solid rgba(0,0,0,0.2);"></span> Consensus Hotspot</li>
                <li><span class="dot" style="background:#800080;border:0.5px solid rgba(0,0,0,0.2);"></span> Ligand</li>
            </ul>
        </div>
    </div>
    <script>
        var appData = VAR_DATA_JSON;
        var viewer = null;

        function updateTelemetry(ac, gl, mode, ligName) { 
            document.getElementById('atom-stat').innerText = ac + " atoms"; 
            document.getElementById('sys-mode').innerText = "READY (" + mode + ")"; 
            var ligElem = document.getElementById('ligand-stat');
            if(ligElem) {
                ligElem.innerText = (ligName && ligName.length > 0) ? "LIGAND: " + ligName : "LIGAND: NONE";
            }
        }

        function getResiIds(resStr) { if(!resStr) return []; var m=String(resStr).match(/\\d+/g); return m?m.map(Number):[]; }

        function addSimpleLabel(v, text, center, bgColor, fontColor) {
            if(!center) return;
            v.addLabel(text, {position:{x:center[0],y:center[1],z:center[2]}, backgroundColor:bgColor||"#222", fontColor:fontColor||"white", fontSize:13, fontFamily:"sans-serif", backgroundOpacity:0.9, padding:4});
        }

        function addEPLabel(v, val, center, resInfo) {
            if(!center||val===undefined) return;
            var isZero = Math.abs(val) <= 0.001;
            var bgCol = isZero ? "#FFFFFF" : (val > 0 ? "#2563eb" : "#ef4444");
            var txtCol = isZero ? "black" : "white";
            var text = resInfo + " | EP " + val.toFixed(2) + " kT/e";
            v.addLabel(text, {position:{x:center[0]+1.5,y:center[1]+1.5,z:center[2]+1.5}, backgroundColor:bgCol, fontColor:txtCol, fontSize:14, backgroundOpacity:0.95});
        }

        var getAdaptiveColorFunc = function(atoms) {
            var maxB=0; for(var i=0;i<atoms.length;i++) { if(atoms[i].b>maxB) maxB=atoms[i].b; }
            var isExp = maxB < 65;
            return function(atom) {
                if(!appData.show_stability) return '#CCCCCC';
                var v=atom.b;
                if(isExp) { if(v<=15) return '#0053D6'; if(v<=30) return '#65CBF3'; if(v<=45) return '#FFDB13'; return '#FF7D45'; }
                else { if(v>=90) return '#0053D6'; if(v>=70) return '#65CBF3'; if(v>=50) return '#FFDB13'; return '#FF7D45'; }
            };
        };

        function initViewer() {
            setTimeout(function() {
                try {
                    if(typeof $3Dmol === 'undefined') return;
                    var c=document.createElement("canvas"); 
                    var gl=c.getContext("webgl")||c.getContext("experimental-webgl");
                    var glStatus = gl ? "WebGL OK" : "WebGL Fail";
                    
                    var viewerEle = $("#viewer3d");
                    viewer = $3Dmol.createViewer(viewerEle, {backgroundColor:"VAR_VIEWER_BG"});
                    viewer.resize();

                    if(!appData.pdb || appData.pdb.length < 10) { $(".loading").text("No PDB Data"); return; }

                    var m = viewer.addModel(appData.pdb, "pdb", {keepH: true});
                    var atomCount = m.selectedAtoms().length;
                    
                    var adaptiveColorFunc = getAdaptiveColorFunc(m.selectedAtoms());
                    var baseStyle = {cartoon:{colorfunc:adaptiveColorFunc, opacity:0.85}};
                    if(appData.show_base_residues) { baseStyle.stick = {colorfunc:adaptiveColorFunc, radius:0.05, opacity:0.4}; }
                    viewer.setStyle({model:0}, baseStyle);


                    // APBS Surface (Individual mode ONLY)
                    if(appData.mode === "individual" && appData.vis_apbs && appData.dx && appData.dx.length > 100) {
                        try {
                            var binaryStr = atob(appData.dx);
                            var bytes = new Uint8Array(binaryStr.length);
                            for(var i=0; i<binaryStr.length; i++) { bytes[i] = binaryStr.charCodeAt(i); }
                            var rawDxStr = pako.ungzip(bytes, {to: 'string'});
                            var voldata = new $3Dmol.VolumeData(rawDxStr, "dx");
                            viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.85, volscheme: new $3Dmol.Gradient.RWB(-5, 5), voldata: voldata}, {model: 0});
                        } catch(eApbs) { console.warn("APBS error:", eApbs); }
                    }

                    // Mode Specific Rendering
                    if(appData.mode === "individual") {
                        if(appData.vis_mcsa && appData.mcsa && appData.mcsa.length > 0) {
                            var m_ids = appData.mcsa.map(r => parseInt(r.res_num)).filter(id => !isNaN(id));
                            viewer.addStyle({model:0, resi:m_ids}, {stick:{color:"#000", radius:0.2}, sphere:{color:"#000", radius:0.35}});
                            appData.mcsa.forEach(mc => {
                                try {
                                    var resNum = parseInt(mc.res_num);
                                    var center = mc.center;
                                    if(!center) {
                                        var sa = viewer.selectedAtoms({model:0, resi:[resNum]});
                                        if(sa && sa.length > 0) {
                                            var sx=0, sy=0, sz=0; 
                                            sa.forEach(a => { sx+=a.x; sy+=a.y; sz+=a.z; });
                                            center = [sx/sa.length, sy/sa.length, sz/sa.length];
                                        }
                                    }
                                    if(center) addSimpleLabel(viewer, mc.res_name + " " + mc.res_num, center, "#222", "white");
                                } catch(eMcsa) { console.warn("M-CSA individual label err:", eMcsa); }
                            });
                        }
                        if(appData.hotspots) {
                            appData.hotspots.forEach(h => {
                                try {
                                    var isMPB=h.label&&(h.label.includes("MPBind")||h.label.includes("AI-Based Binding")), isP2R=h.label&&h.label.includes("P2Rank"), isFP=h.label&&h.label.includes("fPocket");
                                    if(isMPB&&!appData.vis_mpbind) return;
                                    if(isP2R&&!appData.vis_p2rank) return;
                                    if(isFP&&!appData.vis_fpocket) return;
                                    var ids=getResiIds(h.residues), col=isP2R?"#FF5722":(isFP?"#2ECC71":"#00008F");
                                    viewer.addStyle({model:0, resi:ids}, {stick:{color:col, radius:0.2}, sphere:{color:col, radius:0.35}});
                                    if(h.center) viewer.addSphere({center:{x:h.center[0], y:h.center[1], z:h.center[2]}, radius:5.0, color:col, opacity:0.1, wireframe:true});
                                    if(h.detailed_residues) { 
                                        h.detailed_residues.forEach(dr => { if(dr.center) addSimpleLabel(viewer, dr.name+" "+dr.num, dr.center, col, "white"); }); 
                                    }
                                } catch(eInd) { console.warn("Ind mode err:", eInd); }
                            });
                        }
                    } else if(appData.mode === "consensus") {
                        if(appData.consensus) {
                            appData.consensus.forEach(h => {
                                try {
                                    var ids = [parseInt(h.ResNo)];
                                    viewer.addStyle({model:0, resi:ids}, {stick:{color:"#FFD700", radius:0.2}, sphere:{color:"#FFD700", radius:0.35}});
                                    if(h.center) {
                                        if(appData.vis_apbs && h.ep_value !== undefined) addEPLabel(viewer, h.ep_value, h.center, h.ResName+" "+h.ResNo);
                                        else addSimpleLabel(viewer, h.ResName+" "+h.ResNo, h.center, "#FFD700", "black");
                                    }
                                } catch(eCons) { console.warn("Cons mode err:", eCons); }
                            });
                        }
                    } else if(appData.mode === "docking") {
                        if(appData.vis_mcsa) {
                            if(appData.mcsa) {
                                appData.mcsa.forEach(mc => {
                                    var resNo = parseInt(mc.res_num);
                                    viewer.addStyle({model:0, resi:[resNo]}, {stick:{color:"#000", radius:0.2}, sphere:{color:"#000", radius:0.35}});
                                    var center = mc.center;
                                    if(!center) {
                                        var atoms = viewer.selectedAtoms({model:0, resi:[resNo]});
                                        if(atoms.length > 0) {
                                            var x=0,y=0,z=0; atoms.forEach(a=>{x+=a.x;y+=a.y;z+=a.z;});
                                            center = [x/atoms.length, y/atoms.length, z/atoms.length];
                                        }
                                    }
                                    if(center) {
                                        if(appData.vis_apbs && mc.ep_value !== undefined) addEPLabel(viewer, mc.ep_value, center, mc.res_name+" "+mc.res_num);
                                        else addSimpleLabel(viewer, mc.res_name+" "+mc.res_num, center, "#222", "white");
                                    }
                                });
                            }
                        } else {
                            if(appData.consensus) {
                                appData.consensus.forEach(h => {
                                    var resNo = parseInt(h.ResNo);
                                    viewer.addStyle({model:0, resi:[resNo]}, {stick:{color:"#FFD700", radius:0.2}, sphere:{color:"#FFD700", radius:0.35}});
                                    if(h.center) {
                                        if(appData.vis_apbs && h.ep_value !== undefined) addEPLabel(viewer, h.ep_value, h.center, h.ResName+" "+h.ResNo);
                                        else addSimpleLabel(viewer, h.ResName+" "+h.ResNo, h.center, "#FFD700", "black");
                                    }
                                });
                            }
                        }
                        if(appData.ligand && appData.ligand.length > 20) {
                            var ligM = viewer.addModel(appData.ligand, appData.ligand_format || "pdb");
                            viewer.setStyle({model:1}, {stick:{color:"#800080", radius:0.22}});
                            if(appData.ligand_name) {
                                var la = ligM.selectedAtoms(); if(la.length > 0) {
                                    var x=0,y=0,z=0; la.forEach(a=>{x+=a.x;y+=a.y;z+=a.z;});
                                    addSimpleLabel(viewer, appData.ligand_name, [x/la.length, y/la.length, z/la.length], "#800080", "white");
                                }
                            }
                        }
                    }

                    viewer.zoomTo({model:0});
                    viewer.render();
                    setTimeout(function() { viewer.zoomTo({model:0}); viewer.render(); }, 500);
                    updateTelemetry(atomCount, glStatus, appData.mode, appData.ligand_name);
                    $(".loading").fadeOut();

                } catch(eGlobal) { 
                    console.error("Init Error:", eGlobal); 
                    $(".loading").text("ERROR: " + eGlobal.message);
                }
            }, 150);
        }
        $(document).ready(function() { setTimeout(initViewer, 500); });
    </script></body></html>
    """
    replacements = {
        "VAR_DATA_JSON": data_json,
        "VAR_VIEWER_BG": viewer_bg,
        "VAR_TEXT_COLOR": text_color,
        "VAR_PANEL_BG": panel_bg,
        "VAR_PANEL_BDR": panel_bdr,
        "VAR_SHADOW": shadow,
        "VAR_LEGEND_DISPLAY": legend_display
    }
    for k, v in replacements.items():
        html_code = html_code.replace(k, str(v))
    
    components.html(html_code, height=600)

def get_real_mpbind_predictions(pdb_path):
    wrapper = MPBindWrapper(BASE_DIR)
    results = wrapper.predict(pdb_path, os.path.join(BASE_DIR, "results"))
    if not results: return []
    return wrapper.extract_hotspots(results[0])

def get_residue_center(pdb_str, residue_query):
    import numpy as np
    coords = []
    parts = residue_query.split()
    target_chain = parts[0] if len(parts) > 1 else None
    target_res = parts[1] if len(parts) > 1 else parts[0]
    for line in pdb_str.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            r_num = line[22:26].strip()
            chain = line[21].strip()
            if r_num == target_res:
                if target_chain is None or target_chain == chain:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not coords: return None
    return np.mean(coords, axis=0).tolist()

def get_detailed_residues(res_str, info_map):
    if not res_str: return []
    import re
    nums = re.findall(r'(\d+)', str(res_str))
    results = []
    seen = set()
    for r_num in nums:
        if r_num in info_map and r_num not in seen:
            results.append({"name": info_map[r_num]["name"], "num": r_num, "center": info_map[r_num]["center"]})
            seen.add(r_num)
    return results

def run_batch_docking_real(receptor_pdb, ligand_sdf, hotspots):
    validator = VinaValidator(BASE_DIR)
    return validator.run_batch_docking_parallel(receptor_pdb, ligand_sdf, hotspots, max_workers=4)

# === 메인 레이아웃 및 업로드 상태 ===
pdb_loaded = "current_pdb_path" in st.session_state
pdb_data = st.session_state.get("pdb_data", "")

# === 사이드바 구성 ===
with st.sidebar:
    # 초기화 카운터 가져오기 (위젯 강제 리셋용)
    rc = st.session_state.get("reset_count", 0)

    st.markdown('<div class="h-2"></div>', unsafe_allow_html=True)
    st.markdown("""<div class="flex items-center gap-2 mb-5"><span class="material-symbols-outlined text-primary text-xl">visibility</span><h3 class="text-[18px] font-bold uppercase tracking-[0.2em] opacity-40" style="color: var(--text-primary);">시각화 설정</h3></div>""", unsafe_allow_html=True)
    st.session_state.show_stability = st.toggle("구조 안정성 색상 표시", value=st.session_state.get('show_stability', False), key=f"toggle_stability_{rc}")
    st.session_state.show_residues  = st.toggle("구조 전체 잔기 표시", value=st.session_state.get('show_residues', False), key=f"toggle_residues_{rc}")
    
    st.markdown('<div class="h-8"></div>', unsafe_allow_html=True)
    st.markdown("""<div class="flex items-center gap-2 mb-5"><span class="material-symbols-outlined text-primary text-xl">folder_open</span><h3 class="text-[18px] font-bold uppercase tracking-[0.2em] opacity-40" style="color: var(--text-primary);">데이터 업로드</h3></div>""", unsafe_allow_html=True)
    pdb_file = st.file_uploader("항생제 내성 유발 단백질 구조 파일(PDB) [필수]", type=None, label_visibility="visible", key=f"pdb_uploader_{rc}")
    
    st.text_input("UniProt ID 직접 입력 (e.g. P62593) [선택]", key=f"user_uniprot_id_{rc}")
    st.markdown('<div style="margin-top: -15px; margin-bottom: 15px;"><a href="https://www.ebi.ac.uk/thornton-srv/m-csa/" target="_blank" style="font-size: 14px; color: var(--primary-color); text-decoration: none; display: flex; items-center: center; gap: 4px;"><span class="material-symbols-outlined" style="font-size: 14px;">open_in_new</span> Mechanism and Catalytic Site Atlas (M-CSA) </a></div>', unsafe_allow_html=True)
    
    ligand_files = st.file_uploader("리간드(항생제) 구조 파일 (SDF/MOL2/PDB) [선택]", type=["sdf", "mol2", "pdb"], accept_multiple_files=True, label_visibility="visible", key=f"ligand_uploader_{rc}")

    if st.button("🔄 세션 초기화", use_container_width=True):
        next_rc = st.session_state.get("reset_count", 0) + 1
        st.session_state.clear()
        st.session_state["reset_count"] = next_rc
        st.rerun()

    if st.button(":house: 서비스 홈", type="primary", use_container_width=True, key="back_to_main_pselector"):
        # This will be integrated with the main selection screen later
        pass

    if ligand_files:
        paths = []
        for lf in ligand_files: # 입력 순서대로 처리
            l_path = os.path.join(INPUT_DIR, lf.name)
            with open(l_path, "wb") as f: f.write(lf.getbuffer())
            paths.append(l_path)
        
        # 첫 번째 리간드 데이터 로드 (시각화용)
        try:
            st.session_state["ligand_data"] = ligand_files[0].getvalue().decode("utf-8")
            st.session_state["ligand_name"] = ligand_files[0].name
        except: pass

        if st.session_state.get("current_ligand_paths") != paths:
            st.session_state["current_ligand_paths"] = paths
            # 리간드 순서 맵 생성 (파일명 기준)
            st.session_state["ligand_order_map"] = {os.path.basename(p): i for i, p in enumerate(paths)}
            st.rerun()
    else:
        if st.session_state.get("current_ligand_paths") or st.session_state.get("ligand_data"):
            st.session_state["current_ligand_paths"] = []
            st.session_state["ligand_data"] = ""
            st.session_state["ligand_name"] = ""
            st.rerun()

    if pdb_file:
        orig_path = os.path.join(INPUT_DIR, f"orig_{pdb_file.name}")
        clean_path = os.path.join(INPUT_DIR, f"clean_{pdb_file.name}")
        if not os.path.exists(orig_path):
            with open(orig_path, "wb") as f: f.write(pdb_file.getbuffer())
            if not clean_pdb_file(orig_path, clean_path):
                import shutil; shutil.copy2(orig_path, clean_path)
        if st.session_state.get("current_pdb_path") != clean_path:
            # 새로운 파일 업로드 시 시각화 옵션 초기화 (회색 기본 구조만 표시)
            st.session_state["show_stability"] = False
            st.session_state["show_residues"]  = False
            
            st.session_state["current_pdb_path"] = clean_path
            with open(clean_path, "r") as f: st.session_state["pdb_data"] = f.read()
            info = extract_pdb_info(st.session_state.get("pdb_data", ""), pdb_filename=pdb_file.name)
            
            # 메타데이터 추출 및 저장
            st.session_state["extracted_pdb_id"] = info.get("pdb_id")
            st.session_state["extracted_uniprot_id"] = info.get("uniprot_id")
            st.session_state["extracted_keyword"] = info.get("title_keyword")
            
            # 분석 결과 초기화 (삭제 시가 아니라 '다른 파일 업로드' 시에만 초기화)
            for k in ["mcsa_residues", "mp_hotspots", "p2rank_results", "fpocket_results",
                      "consensus_hotspots", "apbs_results", "vina_results_mcsa", "vina_results_consensus"]:
                st.session_state[k] = []
            st.rerun()
    else:
        # PDB 파일 삭제 시 분석 데이터는 유지하고 이미지만 제거
        if st.session_state.get("current_pdb_path"):
            del st.session_state["current_pdb_path"]
            st.session_state["pdb_data"] = ""
            st.rerun()

# === 메인 레이아웃 영역 ===

res_info_map = {}
if pdb_data:
    temp_map = {}
    for line in pdb_data.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            r_num, r_name = line[22:26].strip(), line[17:20].strip()
            c = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if r_num not in temp_map: temp_map[r_num] = {"name": r_name, "coords": []}
            temp_map[r_num]["coords"].append(c)
    for r_num, info in temp_map.items():
        res_info_map[r_num] = {"name": info["name"], "center": np.mean(info["coords"], axis=0).tolist()}

all_hotspots = []
for h in st.session_state.get("mp_hotspots", []):
    c = h.copy(); c["label"] = f"AI-Based Binding {h.get('site')}"; c["detailed_residues"] = get_detailed_residues(h.get("residues"), res_info_map); all_hotspots.append(c)
for p in st.session_state.get("p2rank_results", []):
    c = p.copy(); c["label"] = f"P2Rank {p.get('site')}"; c["detailed_residues"] = get_detailed_residues(p.get("residues"), res_info_map); all_hotspots.append(c)
for f in st.session_state.get("fpocket_results", []):
    c = f.copy(); c["label"] = f"fPocket {f.get('site')}"; c["detailed_residues"] = get_detailed_residues(f.get("residues"), res_info_map); all_hotspots.append(c)

r1_c1, r1_c2 = st.columns(2)
with r1_c2:
    st.markdown("""<div class="flex items-center gap-2 mb-4"><span class="material-symbols-outlined text-primary text-xl">tune</span><h3 class="text-[18px] font-bold uppercase tracking-widest opacity-80">임계값 설정 파라미터</h3></div>""", unsafe_allow_html=True)
    with st.container(border=True):
        ic1, ic2 = st.columns(2)
        with ic1:
            th_mp = st.number_input("MPBind: 결합 확률", 0.0, 1.0, 0.5, 0.05, key=f"cfg_mp_{rc}")
            th_fpc = st.number_input("fPocket: 약물 적합도", 0.0, 1.0, 0.5, 0.05, key=f"cfg_fpc_{rc}")
        with ic2:
            th_p2r = st.number_input("P2Rank: 상위 % (포켓 점수 기반)", 0, 100, 30, 5, key=f"cfg_p2r_{rc}")
            th_apbs = st.number_input("APBS: 상위 % (절대값 EP 기반)", 0, 100, 30, 5, key=f"cfg_apbs_{rc}")
        st.button("🔄 임계값 초기화", use_container_width=True, on_click=reset_parameter_defaults)

with r1_c1:
    st.markdown("""<div class="flex items-center gap-2 mb-4"><span class="material-symbols-outlined text-primary text-xl">bolt</span><h3 class="text-[18px] font-bold uppercase tracking-widest opacity-80">타깃 식별 및 예측 모듈 콘솔</h3></div>""", unsafe_allow_html=True)
    if st.button("🧬 M-CSA 표준 모티프(활성 부위) 식별", use_container_width=True):
        if not pdb_loaded: st.warning("PDB 구조 파일을 먼저 업로드해주세요.")
        else:
            with st.spinner("M-CSA 데이터베이스 탐색 중..."):
                u_id = st.session_state.get(f"user_uniprot_id_{rc}", "").strip()
                ex_u_id = st.session_state.get("extracted_uniprot_id")
                ex_p_id = st.session_state.get("extracted_pdb_id")
                ex_key  = st.session_state.get("extracted_keyword")
                
                res = []
                if u_id:
                    st.write(f"검색 중: UniProt({u_id})...")
                    res = get_mcsa_residues_by_uniprot(u_id, target_pdb_id=ex_p_id)
                
                if not res and ex_u_id:
                    st.write(f"검색 중: UniProt(헤더 추출 {ex_u_id})...")
                    res = get_mcsa_residues_by_uniprot(ex_u_id, target_pdb_id=ex_p_id)
                
                if not res and ex_p_id:
                    st.write(f"🔍 PDB ID({ex_p_id})로 검색 중...")
                    res = get_mcsa_residues_by_pdb_id(ex_p_id)
                
                if not res and u_id:
                    st.write(f"🔍 입력값({u_id})을 PDB ID로 간주하여 검색 중...")
                    res = get_mcsa_residues_by_pdb_id(u_id)
                
                if not res and ex_key:
                    st.write(f"🔍 키워드({ex_key})로 검색 중...")
                    res = get_mcsa_residues_by_keyword(ex_key)
                
                if not res and u_id:
                    st.write(f"🔍 입력값({u_id})을 키워드로 간주하여 검색 중...")
                    res = get_mcsa_residues_by_keyword(u_id)
                
                if res:
                    st.success(f"✅ {len(res)}개의 촉매 잔기를 확인했습니다.")
                    for r in res: r["center"] = get_residue_center(pdb_data, str(r.get("res_num")))
                    st.session_state["mcsa_residues"] = res
                else:
                    st.error("❌ M-CSA 분석 결과가 없습니다. UniProt ID를 직접 입력해 보세요.")
                st.rerun()

    if st.button("🧠 MPBind", use_container_width=True):
        if not pdb_loaded: st.warning("PDB 필요")
        else:
            with st.spinner("AI 분석 중..."):
                mp_res = get_real_mpbind_predictions(st.session_state["current_pdb_path"])
                mp_filtered = sorted([h for h in mp_res if h.get("prob",0)>=th_mp], key=lambda x: x.get("prob",0), reverse=True)
                for i, r in enumerate(mp_filtered): r["site"] = f"Site {i+1}"
                st.session_state["mp_hotspots"] = mp_filtered; st.rerun()

    if st.button("🔍 P2Rank", use_container_width=True):
        if not pdb_loaded: st.warning("PDB 필요")
        else:
            with st.spinner("P2Rank 실행..."):
                p2r_pockets, p2r_res_scores = run_p2rank(st.session_state["current_pdb_path"], RESULT_DIR)
                import logging
                logging.warning(f"STREAMLIT P2RANK RETURN: pockets={len(p2r_pockets) if isinstance(p2r_pockets, list) else type(p2r_pockets)}, scores={p2r_res_scores}")
                if isinstance(p2r_pockets, list) and len(p2r_pockets) > 0:
                    # p가 dict인지 보장 (AttributeError 방지)
                    p2r_pockets = [p for p in p2r_pockets if isinstance(p, dict)]
                    if p2r_pockets:
                        mx = max(p.get("score",0) for p in p2r_pockets)
                        p2r_filtered = sorted([p for p in p2r_pockets if p.get("score",0) >= mx * (1.0 - (th_p2r/100.0))], key=lambda x: x.get("score",0), reverse=True)
                        for i, r in enumerate(p2r_filtered): r["site"] = f"Site {i+1}"
                        st.session_state["p2rank_results"] = p2r_filtered
                        st.session_state["p2rank_residue_scores"] = p2r_res_scores
                        st.rerun()
                else:
                    err_msg = p2r_res_scores.get("error", "분석 결과가 없거나 실행 중 오류가 발생했습니다.")
                    st.error(f"P2Rank 오류: {err_msg}")

    if st.button("🧪 fPocket", use_container_width=True):
        if not pdb_loaded: st.warning("PDB 필요")
        else:
            with st.spinner("fPocket 실행..."):
                fpc = run_fpocket(st.session_state["current_pdb_path"], RESULT_DIR)
                res = sorted([f for f in fpc if f.get("score",0)>=th_fpc], key=lambda x: x.get("score",0), reverse=True)
                for i, r in enumerate(res): r["site"] = f"Site {i+1}"
                st.session_state["fpocket_results"] = res; st.rerun()

    if st.button("⚡ APBS", use_container_width=True):
        if not pdb_loaded: st.warning("PDB 필요")
        else:
            with st.spinner("APBS 계산..."):
                dx = run_apbs_pipeline(st.session_state["current_pdb_path"], RESULT_DIR)
                if dx: st.session_state["apbs_results"] = [dx]; st.session_state["v_apbs"] = True
                st.rerun()

    st.markdown('<div class="gold-btn mt-4">', unsafe_allow_html=True)
    if st.button("🏆 다차원 통합 모듈 엔진", use_container_width=True):
        if not pdb_loaded: st.warning("PDB 필요")
        else:
            with st.status("🚀 통합 분석 파이프라인 가동..."):
                # 1. M-CSA
                st.write("🧬 M-CSA 검색...")
                # 1. M-CSA
                st.write("🧬 M-CSA 표준 모티프 검색 중...")
                u_id = st.session_state.get(f"user_uniprot_id_{rc}", "").strip()
                ex_u_id = st.session_state.get("extracted_uniprot_id")
                ex_p_id = st.session_state.get("extracted_pdb_id")
                ex_key  = st.session_state.get("extracted_keyword")
                
                m = []
                if u_id: m = get_mcsa_residues_by_uniprot(u_id, target_pdb_id=ex_p_id)
                if not m and ex_u_id: m = get_mcsa_residues_by_uniprot(ex_u_id, target_pdb_id=ex_p_id)
                if not m and ex_p_id: m = get_mcsa_residues_by_pdb_id(ex_p_id)
                if not m and u_id:    m = get_mcsa_residues_by_pdb_id(u_id)
                if not m and ex_key:  m = get_mcsa_residues_by_keyword(ex_key)
                if not m and u_id:    m = get_mcsa_residues_by_keyword(u_id)
                
                for r in m: r["center"] = get_residue_center(pdb_data, str(r.get("res_num")))
                st.session_state["mcsa_residues"] = m
                
                # 2. MPBind
                st.write("🧠 AI 추론...")
                mp_res = get_real_mpbind_predictions(st.session_state["current_pdb_path"])
                mp_filtered = sorted([h for h in mp_res if h.get("prob",0)>=th_mp], key=lambda x: x.get("prob",0), reverse=True)
                for i, r in enumerate(mp_filtered): r["site"] = f"Site {i+1}"
                st.session_state["mp_hotspots"] = mp_filtered
                
                # 3. P2Rank
                st.write("🔍 P2Rank 포켓 검색...")
                p2r_pockets, p2r_res_scores = run_p2rank(st.session_state["current_pdb_path"], RESULT_DIR)
                if isinstance(p2r_pockets, list) and len(p2r_pockets) > 0:
                    # p가 dict인지 보장
                    p2r_pockets = [p for p in p2r_pockets if isinstance(p, dict)]
                    if p2r_pockets:
                        mx = max(p.get("score",0) for p in p2r_pockets)
                        p2r_filtered = sorted([p for p in p2r_pockets if p.get("score",0) >= mx * (1.0 - (th_p2r/100.0))], key=lambda x: x.get("score",0), reverse=True)
                        for i, r in enumerate(p2r_filtered): r["site"] = f"Site {i+1}"
                        st.session_state["p2rank_results"] = p2r_filtered
                        st.session_state["p2rank_residue_scores"] = p2r_res_scores
                else:
                    err_msg = p2r_res_scores.get("error", "분석 실패")
                    st.error(f"P2Rank 에러: {err_msg}")
                    st.stop()
                
                # 4. fPocket
                st.write("🧪 fPocket 분석...")
                fpc_res = run_fpocket(st.session_state["current_pdb_path"], RESULT_DIR)
                fpc_filtered = sorted([f for f in fpc_res if f.get("score",0)>=th_fpc], key=lambda x: x.get("score",0), reverse=True)
                for i, r in enumerate(fpc_filtered): r["site"] = f"Site {i+1}"
                st.session_state["fpocket_results"] = fpc_filtered
                
                # 5. APBS
                st.write("⚡ APBS 정전기장 계산...")
                dx = run_apbs_pipeline(st.session_state["current_pdb_path"], RESULT_DIR)
                if dx:
                    st.session_state["apbs_results"] = [dx]
                    for key in ["mp_hotspots", "p2rank_results", "fpocket_results", "mcsa_residues"]:
                        upd = []
                        for p in st.session_state.get(key, []):
                            p_c = p.copy()
                            res_str = p.get("residues", "")
                            if key == "mcsa_residues": res_str = f"{p['res_name']} {p['res_num']}"
                            ns, rs = score_pocket_electrostatics(res_str, dx, st.session_state["current_pdb_path"])
                            p_c["apbs_score"], p_c["raw_apbs_score"] = ns, rs
                            if key == "mcsa_residues": p_c["ep_value"] = rs
                            upd.append(p_c)
                        st.session_state[key] = upd
                    st.session_state["v_apbs"] = True

                # 6. Consensus
                st.write("🏆 통합 결합 부위 도출 중...")
                engine = ConsensusEngine(st.session_state["current_pdb_path"], st.session_state.get("mp_hotspots",[]), st.session_state.get("p2rank_results",[]), st.session_state.get("fpocket_results",[]), st.session_state.get("mcsa_residues",[]), 
                                         p2r_res_scores_dict=st.session_state.get("p2rank_residue_scores",{}), 
                                         p2r_res_scores=st.session_state.get("p2rank_residue_scores",{}), 
                                         mp_th=th_mp, p2r_pct=th_p2r, fpc_th=th_fpc, apbs_pct=th_apbs)
                final = engine.run()
                
                rfd, res_f = engine.generate_report(final, os.path.join(RESULT_DIR, "Final_Hotspots.txt"))
                st.session_state["consensus_hotspots"] = final; st.session_state["consensus_formats"] = {"rfd": rfd, "res": res_f}
            st.rerun()
    st.markdown('</div>', unsafe_allow_html=True)

st.markdown("---")
st.markdown("""<div class="flex items-center gap-2 mb-4"><span class="material-symbols-outlined text-primary text-xl">view_in_ar</span><h2 class="text-[25px] font-bold opacity-90">단백질/리간드 3D 구조 뷰어</h2></div>""", unsafe_allow_html=True)

vt1, vt2, vt3 = st.tabs(["🧬 모듈별 분석 결과", "🏆 통합 모듈 분석 결과", "⚓ 도킹 검증"])
with vt1:
    if pdb_loaded:
        c1, c2, c3, c4, c5 = st.columns(5)
        v_mcsa = c1.checkbox("M-CSA", key=f"vis_mcsa_toggle_{rc}")
        v_mp   = c2.checkbox("MPBind", key=f"vis_mp_toggle_{rc}")
        v_p2   = c3.checkbox("P2Rank", key=f"vis_p2r_toggle_{rc}")
        v_fp   = c4.checkbox("fPocket", key=f"vis_fp_toggle_{rc}")
        v_ap   = c5.checkbox("APBS", key=f"vis_apbs_toggle_{rc}")
        
        _dx_local = ""
        v_ap_state = st.session_state.get("v_apbs", False) or v_ap
        if v_ap_state and st.session_state.get("apbs_results"):
            try:
                import gzip, base64
                with open(st.session_state["apbs_results"][0], "rb") as f:
                    _dx_compressed = gzip.compress(f.read())
                    _dx_local = base64.b64encode(_dx_compressed).decode("utf-8")
            except Exception as e: 
                st.error(f"APBS 로드 에러: {e}")

        render_3d_view(pdb_data, hotspots=all_hotspots, mcsa_data=st.session_state.get("mcsa_residues",[]), 
                       mode="individual", vis_mcsa=v_mcsa, vis_mpbind=v_mp, vis_p2rank=v_p2, 
                       vis_fpocket=v_fp, vis_apbs=v_ap, dx_str=_dx_local,
                       ligand_str=st.session_state.get("ligand_data", ""),
                       ligand_name=st.session_state.get("ligand_name", ""),
                       show_stability=st.session_state.show_stability, 
                       show_base_residues=st.session_state.show_residues)
with vt2:
    if pdb_loaded:
        render_3d_view(pdb_data, consensus=st.session_state.get("consensus_hotspots",[]), mode="consensus",
                       vis_apbs=st.session_state.get("v_apbs", False), dx_str=_dx_local,
                       ligand_str=st.session_state.get("ligand_data", ""),
                       ligand_name=st.session_state.get("ligand_name", ""),
                       show_stability=st.session_state.show_stability, 
                       show_base_residues=st.session_state.show_residues)
with vt3:
    if pdb_loaded:
        res_m = st.session_state.get("vina_results_mcsa", [])
        res_c = st.session_state.get("vina_results_consensus", [])
        
        if not res_m and not res_c:
            st.info("도킹 분석 결과가 없습니다. 상단의 하위 모듈에서 결합력 평가(도킹)를 먼저 수행해주세요.")
        else:
            cb1, cb2 = st.columns(2)
            with cb1:
                base_type = st.radio("분석 기반 선택", ["M-CSA 표준 모티프 기반", "최종 통합 타겟 핫스팟 기반"], horizontal=True, key=f"base_type_{rc}")
            
            target_results = res_m if base_type == "M-CSA 표준 모티프 기반" else res_c
            
            if not target_results:
                st.warning(f"선택하신 '{base_type}'에 대한 도킹 결과가 존재하지 않습니다.")
            else:
                with cb2:
                    # 입력 순서(order_idx)대로 정렬하여 제공
                    target_results.sort(key=lambda x: x.get("order_idx", 999))
                    lig_names = [r.get("리간드 이름", r.get("리간드", "Unknown")) for r in target_results]
                    selected_lig_name = st.radio("리간드(항생제) 선택", lig_names, horizontal=True, key=f"sel_lig_dock_{rc}")
                
                # 선택된 리간드 결과 추출
                sel_res = next((r for r in target_results if r.get("리간드 이름", r.get("리간드", "Unknown")) == selected_lig_name), None)
                
                if sel_res and os.path.exists(sel_res["docked_pdb"]):
                    with open(sel_res["docked_pdb"], "r") as f:
                        _lig_str = f.read()
                    
                    # 뷰어 렌더링
                    render_3d_view(pdb_data, 
                                   mcsa_data=st.session_state.get("mcsa_residues",[]), 
                                   consensus=st.session_state.get("consensus_hotspots",[]), 
                                   mode="docking", 
                                   vis_mcsa=(base_type == "M-CSA 표준 모티프 기반"), 
                                   vis_apbs=st.session_state.get("v_apbs", False), 
                                   dx_str="", ligand_str=_lig_str, 
                                   ligand_name=selected_lig_name,
                                   show_stability=st.session_state.show_stability, 
                                   show_base_residues=st.session_state.show_residues)

st.markdown("---")
r3_c1, r3_c2 = st.columns([0.55, 0.45])
with r3_c1:
    st.markdown('<h3 class="text-[25px] font-bold uppercase tracking-widest opacity-80 mb-4">단위 모듈 분석 결과 패널</h3>', unsafe_allow_html=True)
    tabs = st.tabs(["M-CSA", "MPBind", "P2Rank", "fPocket"])
    rm = get_pdb_residue_map(pdb_data)
    with tabs[0]: 
        if st.session_state.get("mcsa_residues"):
            mc_html = '<div class="overflow-x-auto"><table class="final-results-table"><thead><tr><th>잔기 ID</th><th>아미노산 정보</th><th>역할</th></tr></thead><tbody>'
            for m in st.session_state["mcsa_residues"]: mc_html += f'<tr><td>A{m["res_num"]}</td><td>{m["res_name"]} {m["res_num"]}</td><td>{m.get("role","Catalytic")}</td></tr>'
            st.markdown(mc_html+"</tbody></table></div>", unsafe_allow_html=True)
    with tabs[1]: st.markdown(generate_styled_table(st.session_state.get("mp_hotspots",[]), ["site","prob","residues"], rm, 
                                                  col_map={"site":"사이트", "prob":"결합 확률", "residues":"잔기 ID"}, widths=[15,20,35,30]), unsafe_allow_html=True)
    with tabs[2]: st.markdown(generate_styled_table(st.session_state.get("p2rank_results",[]), ["site","score","residues"], rm, 
                                                  col_map={"site":"사이트", "score":"포켓 점수", "residues":"잔기 ID"}, widths=[15,25,35,30]), unsafe_allow_html=True)
    with tabs[3]: st.markdown(generate_styled_table(st.session_state.get("fpocket_results",[]), ["site","score","residues"], rm, 
                                                  col_map={"site":"사이트", "score":"약물 적합도", "residues":"잔기 ID"}, widths=[15,25,35,30]), unsafe_allow_html=True)

with r3_c2:
    st.markdown('<h3 class="text-[25px] font-bold uppercase tracking-widest opacity-80 mb-4">다차원 통합 모듈 분석 결과</h3>', unsafe_allow_html=True)
    if st.session_state.get("consensus_hotspots"):
        table_html = '<div class="overflow-x-auto"><table class="final-results-table"><thead><tr><th>잔기 ID</th><th>아미노산 정보</th><th>분석 근거(모듈)</th><th>통합 신뢰 점수</th></tr></thead><tbody>'
        rf_list = []
        seq_list = []
        for item in st.session_state["consensus_hotspots"]: 
            is_anchor = "M-CSA" in item.get("Tools", "")
            anchor_tag = '<span style="color: #ef4444; margin-left: 6px; font-size: 0.85em; font-weight: bold;">(Anchor)</span>' if is_anchor else ""
            table_html += f'<tr><td>A{item["ResNo"]}</td><td>{item["ResName"]} {item["ResNo"]}</td><td>{item["Tools"]}</td><td class="text-primary font-bold">{item["Score"]}{anchor_tag}</td></tr>'
            rf_list.append(f"A{item['ResNo']}")
            seq_list.append(f"{item['ResName']} {item['ResNo']}")
        st.markdown(table_html+"</tbody></table></div>", unsafe_allow_html=True)
        
        st.markdown('<h4 class="text-[18px] font-bold uppercase tracking-widest opacity-80 mt-6 mb-3">최종 통합 타겟 핫스팟 정보</h4>', unsafe_allow_html=True)
        rf_str = ",".join(rf_list)
        seq_str = ", ".join(seq_list)
        
        info_l, info_r = st.columns(2)
        with info_l:
            st.markdown('<div class="text-[18px] font-bold text-gray-500 mb-1">RFdiffusion 포맷</div>', unsafe_allow_html=True)
            st.code(rf_str, language="text")
            st.download_button("⬇️ RFdiffusion 포맷 다운로드 (.txt)", rf_str, "rfdiffusion_format.txt", use_container_width=True)
        with info_r:
            st.markdown('<div class="text-[18px] font-bold text-gray-500 mb-1">연구자용 서열 정보</div>', unsafe_allow_html=True)
            st.code(seq_str, language="text")
            st.download_button("⬇️ 서열 정보 다운로드 (.txt)", seq_str, "sequence_info.txt", use_container_width=True)

st.markdown("---")
st.markdown('<h3 class="text-[25px] font-bold uppercase tracking-widest opacity-80 mb-4">⚓ AutoDock Vina: 결합력 평가</h3>', unsafe_allow_html=True)
v1, v2 = st.columns(2)
with v1:
    if st.button("⚓ M-CSA 표준 모티프 기반", use_container_width=True):
        if not pdb_loaded or not st.session_state.get("current_ligand_paths"): st.warning("파일 필요")
        else:
            with st.spinner("도킹 중..."):
                coords = [get_residue_center(pdb_data, str(m['res_num'])) for m in st.session_state.get("mcsa_residues",[])]
                center = np.mean([c for c in coords if c], axis=0).tolist() if any(coords) else None
                if center:
                    res = []
                    mcsa_hotspots = [f"{h['ResName']} {h['ResNo']}" for h in st.session_state.get("consensus_hotspots",[]) if "M-CSA" in h.get("Tools", "")]
                    mcsa_hotspot_str = ", ".join(mcsa_hotspots)
                    for lp in st.session_state["current_ligand_paths"]:
                        d = run_batch_docking_real(st.session_state["current_pdb_path"], lp, [{"center": center}])
                        for s, r in d.items(): 
                            l_name = os.path.splitext(os.path.basename(lp))[0]
                            res.append({
                                "리간드 이름": l_name, 
                                "M-CSA 표준 모티프": mcsa_hotspot_str, 
                                "결합 에너지": f"{r.get('affinity',0):.3f} kcal/mol", 
                                "docked_pdb": r.get('output_pdb'),
                                "order_idx": st.session_state.get("ligand_order_map", {}).get(os.path.basename(lp), 999)
                            })
                    st.session_state["vina_results_mcsa"] = res; st.rerun()
with v2:
    if st.button("⚓ 최종 통합 타겟 핫스팟 기반", use_container_width=True):
        if not pdb_loaded or not st.session_state.get("current_ligand_paths"): st.warning("파일 필요")
        else:
            with st.spinner("도킹 중..."):
                coords = [get_residue_center(pdb_data, str(h['ResNo'])) for h in st.session_state.get("consensus_hotspots",[])]
                center = np.mean([c for c in coords if c], axis=0).tolist() if any(coords) else None
                if center:
                    res = []
                    cons_hotspots = [f"{h['ResName']} {h['ResNo']}" for h in st.session_state.get("consensus_hotspots",[])]
                    cons_hotspot_str = ", ".join(cons_hotspots)
                    for lp in st.session_state["current_ligand_paths"]:
                        d = run_batch_docking_real(st.session_state["current_pdb_path"], lp, [{"center": center}])
                        for s, r in d.items(): 
                            res.append({
                                "리간드 이름": os.path.splitext(os.path.basename(lp))[0], 
                                "타겟 핫스팟": cons_hotspot_str, 
                                "결합 에너지": f"{r.get('affinity',0):.3f} kcal/mol", 
                                "docked_pdb": r.get('output_pdb')
                            })
                    st.session_state["vina_results_consensus"] = res; st.rerun()

if st.session_state.get("vina_results_mcsa") or st.session_state.get("vina_results_consensus"):
    c1, c2 = st.columns(2)
    with c1: 
        if st.session_state.get("vina_results_mcsa"): 
            st.markdown('<div class="text-sm font-bold opacity-80 mb-2">M-CSA 표준 모티프 기반 결합력 평가 결과</div>', unsafe_allow_html=True)
            df_mcsa = pd.DataFrame(st.session_state["vina_results_mcsa"])
            if "order_idx" in df_mcsa.columns:
                df_mcsa = df_mcsa.sort_values("order_idx").drop(columns=["order_idx"])
            st.table(df_mcsa.drop(columns=["docked_pdb"]))
    with c2: 
        if st.session_state.get("vina_results_consensus"): 
            st.markdown('<div class="text-sm font-bold opacity-80 mb-2">최종 통합 타겟 핫스팟 기반 결합력 평가 결과</div>', unsafe_allow_html=True)
            df_cons = pd.DataFrame(st.session_state["vina_results_consensus"])
            if "order_idx" in df_cons.columns:
                df_cons = df_cons.sort_values("order_idx").drop(columns=["order_idx"])
            st.table(df_cons.drop(columns=["docked_pdb"]))
