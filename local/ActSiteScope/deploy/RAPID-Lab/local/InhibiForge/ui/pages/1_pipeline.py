"""
pages/1_pipeline.py
파이프라인 설정 및 실행 페이지 — Re-Off_Action_UI_code_Light_v2.html 대응.
레이아웃: 설계 파라미터(흰 박스) + 내성 유발 단백질 PDB 업로드 2단 배치
"""
import re
import sys
from pathlib import Path
import streamlit as st

# ── 경로 부트스트랩 ─────────────────────────────────────────────────────────
_app_dir = Path(__file__).resolve().parents[1]
for _p in (str(_app_dir), str(_app_dir.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from lib.session import (
    init_session, RFDIFF_INPUT_DIR, RFDIFF_OUTPUT_BASE,
    get_rf_output_dir, get_ld_output_dir,
)
from lib.html_loader import load_and_inject_css
from lib.viewer import render_pdb_viewer, VIEWER_STYLES, parse_hotspots
from lib.pipeline_runner import run_pipeline
from lib.page_header import render_page_header
from lib.utils import parse_hotspot_input

init_session()
load_and_inject_css(_app_dir / "styles" / "theme.css")

render_page_header("InhibiForge", "InhibiForge 설정 및 실행")

# ── 첫 번째 페이지 글씨 1.5배 확대 CSS ───────────────────────────────────────
st.markdown("""
<style>
/* 페이지 헤더 */
.ifg-page-title    { font-size: 2.25em !important; }
.ifg-page-subtitle { font-size: 2.25em !important; }

/* Streamlit 위젯 라벨 */
[data-testid="stMainBlockContainer"] label,
[data-testid="stMainBlockContainer"] .stMarkdown p,
[data-testid="stMainBlockContainer"] .stMarkdown li,
[data-testid="stMainBlockContainer"] [data-testid="stWidgetLabel"] p {
  font-size: 1.3125rem !important;   /* 0.875rem × 1.5 */
}

/* h3 / h4 태그 */
[data-testid="stMainBlockContainer"] h3 {
  font-size: 2.1rem !important;      /* 1.4rem × 1.5 */
}
[data-testid="stMainBlockContainer"] h4 {
  font-size: 1.8rem !important;      /* 1.2rem × 1.5 */
}

/* 입력 위젯 내부 텍스트 */
[data-testid="stMainBlockContainer"] .stTextInput input,
[data-testid="stMainBlockContainer"] .stNumberInput input,
[data-testid="stMainBlockContainer"] .stTextArea textarea,
[data-testid="stMainBlockContainer"] .stSelectbox [data-baseweb="select"] {
  font-size: 1.3125rem !important;
}

/* 슬라이더 라벨 / 값 */
[data-testid="stMainBlockContainer"] [data-testid="stSlider"] label {
  font-size: 1.275rem !important;    /* 0.85rem × 1.5 */
}

/* 버튼 텍스트 */
[data-testid="stMainBlockContainer"] .stButton button {
  font-size: 1.3125rem !important;
}

/* caption */
[data-testid="stMainBlockContainer"] [data-testid="stCaptionContainer"] p {
  font-size: 1.2rem !important;      /* 0.8rem × 1.5 */
}

/* success / info / error 알림 */
[data-testid="stMainBlockContainer"] [data-testid="stAlert"] p {
  font-size: 1.3125rem !important;
}

/* file uploader 텍스트 */
[data-testid="stMainBlockContainer"] [data-testid="stFileUploader"] span,
[data-testid="stMainBlockContainer"] [data-testid="stFileUploader"] p,
[data-testid="stMainBlockContainer"] [data-testid="stFileUploader"] small {
  font-size: 1.2rem !important;
}

/* selectbox 옵션 텍스트 */
[data-testid="stMainBlockContainer"] [data-baseweb="select"] span {
  font-size: 1.3125rem !important;
}

/* popover / expander 내부 */
[data-testid="stMainBlockContainer"] [data-testid="stExpander"] summary {
  font-size: 1.3125rem !important;
}

/* number input 증감 버튼 옆 값 */
[data-testid="stMainBlockContainer"] .stNumberInput [data-baseweb="input"] {
  font-size: 1.3125rem !important;
}

/* 텍스트 입력 placeholder 크기 축소 (0.5배) */
[data-testid="stMainBlockContainer"] .stTextInput input::placeholder,
[data-testid="stMainBlockContainer"] .stTextArea textarea::placeholder {
  font-size: 1rem !important;
}
</style>
""", unsafe_allow_html=True)

st.markdown("<div style='margin-top:3rem;'></div>", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# 도움말 상수
# ─────────────────────────────────────────────────────────────────────────────
HOTSPOT_HELP_TEXT = """\
💡 Hotspot 잔기 입력 방법

입력 형식: A70, A73, A166  (체인ID + 잔기번호, 쉼표 구분)

🔄 아미노산 3글자 이름도 자동 변환!
  Ser70         →  A70
  Lys73, Glu166 →  A73, A166
  SER70 / ser70 →  A70   (대소문자 무관)
  Ser 70        →  A70   (공백 있어도 OK)
  A70, A73      →  그대로 유지

⚠️ 권장 개수: 3~6개
  → 너무 많으면 오히려 품질이 떨어질 수 있어요
  → RFdiffusion은 지정한 것보다 더 많은 접촉을
     스스로 만들도록 학습되어 있어요

✅ 좋은 hotspot 고르는 기준
  1. 바인더가 동시에 닿을 수 있는 인근 잔기
  2. 소수성(Hydrophobic) 잔기 우선
     (Ala, Val, Ile, Leu, Phe, Trp 등)
  3. Active site 중심부 잔기 2~3개

🔍 잔기 번호 찾는 방법
  → RCSB PDB: rcsb.org/structure/{PDB_ID}
     "Annotations" 탭 → Active Site 확인
  → PDBsum: ebi.ac.uk/pdbsum/{PDB_ID}
     "Ligands" → 접촉 잔기 다이어그램

예시 (TEM-1 베타락타마제, 1M40):
  활성부위 잔기: Ser70, Lys73, Glu166, Lys234
  입력: Ser70, Lys73, Glu166  또는  A70, A73, A166
  → 둘 다 자동으로 A70, A73, A166 으로 변환됩니다.

비워두면: RFdiffusion이 전체 표면을 탐색
  (시간 더 걸리고 결과 품질 낮아질 수 있음)
"""

# ─────────────────────────────────────────────────────────────────────────────
# PDB 파싱 헬퍼
# ─────────────────────────────────────────────────────────────────────────────

def _parse_pdb_chains(pdb_path: str) -> dict:
    """ATOM 레코드에서 체인별 (첫 잔기 번호, 마지막 잔기 번호) 반환."""
    chains: dict[str, list[int]] = {}
    try:
        with open(pdb_path, "r", errors="replace") as _f:
            for _line in _f:
                if not _line.startswith("ATOM"):
                    continue
                _chain = _line[21]
                try:
                    _res = int(_line[22:26])
                except ValueError:
                    continue
                if _chain not in chains:
                    chains[_chain] = [_res, _res]
                else:
                    if _res < chains[_chain][0]:
                        chains[_chain][0] = _res
                    if _res > chains[_chain][1]:
                        chains[_chain][1] = _res
    except Exception:
        return {}
    return {c: (v[0], v[1]) for c, v in sorted(chains.items())}


def _build_contig(chains: dict, old_contig: str) -> str:
    """체인 범위로 RFdiffusion contig 문자열 생성. 기존 바인더 길이 유지."""
    _m = re.search(r"/0\s+([\d]+-[\d]+)", old_contig)
    _binder = _m.group(1) if _m else "80-80"
    _parts = [f"{c}{first}-{last}" for c, (first, last) in chains.items()]
    return "/".join(_parts) + f"/0 {_binder}"


# ─────────────────────────────────────────────────────────────────────────────
# PDB 업로드 사전 처리 (columns 렌더링 전 — session_state로 위젯 값 미리 접근)
# Streamlit은 위젯 렌더링 전에도 session_state[key]로 현재 업로드 파일 접근 가능
# ─────────────────────────────────────────────────────────────────────────────

_upload_obj   = st.session_state.get("pdb_uploader")      # 현재 업로드된 파일 객체
saved_pdb_path = st.session_state.get("ifg_input_pdb_path") or ""

if _upload_obj is not None:
    RFDIFF_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    _sp   = RFDIFF_INPUT_DIR / _upload_obj.name
    _buf  = _upload_obj.getbuffer()
    if not _sp.exists() or _sp.stat().st_size != len(_buf):
        with open(_sp, "wb") as _wf:
            _wf.write(bytes(_buf))
    saved_pdb_path = str(_sp)
    st.session_state["ifg_input_pdb_path"] = saved_pdb_path

# 새 PDB일 때만 체인 파싱 + Contigs 자동 생성
if saved_pdb_path and saved_pdb_path != st.session_state.get("ifg_last_parsed_pdb", ""):
    if Path(saved_pdb_path).exists():
        _chains = _parse_pdb_chains(saved_pdb_path)
        if _chains:
            _new_contig = _build_contig(_chains, st.session_state.get("ifg_contigs", ""))
            st.session_state["ifg_contigs"]       = _new_contig
            st.session_state["input_contigs"]     = _new_contig   # 위젯 값 갱신
            st.session_state["ifg_pdb_chain_info"] = _chains
        st.session_state["ifg_last_parsed_pdb"] = saved_pdb_path

# ─────────────────────────────────────────────────────────────────────────────
# 섹션 A+B: 타겟 PDB 업로드(좌) + 설계 파라미터(우, 흰 박스) — 2단 배치
# ─────────────────────────────────────────────────────────────────────────────
main_left, main_right = st.columns([3, 4], gap="large")


def _normalize_hotspot(raw: str) -> str:
    """핫스팟 입력을 표준형으로 정규화.

    - 구분자는 ',' 만 사용 (공백/세미콜론 등은 ','로 치환)
    - 각 토큰 내부 공백 제거 후 체인 문자는 대문자로 (예: 'a 450' → 'A450')
    """
    if not raw:
        return ""
    # 콤마 외의 구분자(세미콜론, 줄바꿈)는 콤마로 치환
    _raw = re.sub(r"[;\n\r\t]+", ",", raw)
    tokens = [t for t in _raw.split(",")]
    out: list[str] = []
    for tok in tokens:
        tok = re.sub(r"\s+", "", tok)  # 토큰 내부 공백 제거
        if not tok:
            continue
        m = re.match(r"^([A-Za-z])(\d+)$", tok)
        if m:
            out.append(f"{m.group(1).upper()}{m.group(2)}")
        else:
            out.append(tok.upper())
    return ", ".join(out)


def _on_hotspot_change() -> None:
    """텍스트 입력이 바뀔 때 핫스팟을 A450 형태로 정규화."""
    st.session_state["input_hotspot"] = _normalize_hotspot(
        st.session_state.get("input_hotspot", "")
    )


def _on_job_name_change() -> None:
    """새 Job Name 입력 시 파이프라인 상태 및 타이밍 초기화."""
    st.session_state["ifg_pipeline_status"] = {
        "step_rfdiff":    "waiting",
        "step_split":     "waiting",
        "step_lightdock": "waiting",
        "step_1":         "waiting",
        "step_2":         "waiting",
        "step_3":         "waiting",
    }
    st.session_state["ifg_pipeline_elapsed"]     = {}
    st.session_state["ifg_pipeline_start_times"] = {}
    st.session_state["ifg_pipeline_end_times"]   = {}
    st.session_state["ifg_overall_status"]       = "waiting"
    st.session_state["ifg_job_name"] = st.session_state["input_job_name"]

# ── 왼쪽: 타겟 PDB 업로드 (흰 박스) ───────────────────────────────────────
with main_left:
    st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">upload_file</span>내성 유발 단백질 PDB 파일 업로드</h3>', unsafe_allow_html=True)

    with st.container(border=True):
        st.file_uploader(
            "내성 단백질 구조 파일 (타겟 리셉터, .pdb, 최대 5MB)",
            type=["pdb"],
            help="드래그 앤 드롭 또는 클릭하여 파일을 선택하세요.",
            key="pdb_uploader",
        )

        # 파일 저장은 columns 렌더링 전에 완료됨 — 여기서는 상태만 표시
        if _upload_obj is not None:
            st.success(f"✅ `{_upload_obj.name}` → `{RFDIFF_INPUT_DIR / _upload_obj.name}` 저장 완료")

        if saved_pdb_path and Path(saved_pdb_path).exists():
            st.info(f"현재 입력 PDB: `{Path(saved_pdb_path).name}`")
        else:
            # 기존 파일 선택
            existing = sorted(RFDIFF_INPUT_DIR.glob("*.pdb"))
            if existing:
                sel = st.selectbox(
                    "또는 기존 업로드 파일 선택",
                    options=["(없음)"] + [p.name for p in existing],
                    key="existing_pdb_sel",
                )
                if sel != "(없음)":
                    _new_path = str(RFDIFF_INPUT_DIR / sel)
                    if _new_path != saved_pdb_path:
                        st.session_state["ifg_input_pdb_path"] = _new_path
                        st.rerun()  # 다음 렌더에서 Contigs 자동 생성

        # 체인 정보 표시 (PDB 파싱 성공 시)
        _chain_info = st.session_state.get("ifg_pdb_chain_info", {})
        if _chain_info and saved_pdb_path:
            _info_parts = [
                f"Chain {c}: {first}–{last} ({last - first + 1}잔기)"
                for c, (first, last) in sorted(_chain_info.items())
            ]
            #st.caption("🔍 잔기 범위: " + " | ".join(_info_parts))

# ── 오른쪽: 설계 파라미터 (흰 박스) ────────────────────────────────────────────
with main_right:
    st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">tune</span> 설계 파라미터</h3>', unsafe_allow_html=True)
    
    with st.container(border=True):
        col_l, col_r = st.columns(2, gap="medium")

        with col_l:
            job_name = st.text_input(
                "프로젝트 이름",
                value=st.session_state.ifg_job_name,
                help="결과 폴더명으로 사용됩니다. 영문·숫자·밑줄만 사용하세요.",
                key="input_job_name",
                on_change=_on_job_name_change,
            )

            contigs = st.text_input(
                "컨티그(Contigs)",
                value=st.session_state.ifg_contigs,
                placeholder="내성 유발 단백질 PDB 파일 업로드 시 자동 생성",
                help="확산(diffusion) 과정에서 표적 단백질의 어느 부분을 고정할 것인지, 그리고 생성할 결합 단백질(binder)의 길이를 얼마로 할 것인지 지정하는 역할",
                key="input_contigs",
            )

            hs_label_col, hs_help_col = st.columns([2, 1], gap="small")
            with hs_label_col:
                st.markdown("**핫스팟 잔기**")
            with hs_help_col:
                with st.popover("💡 도움말", use_container_width=True):
                    st.markdown("### Hotspot 잔기 입력 방법 도움말")
                    st.code(HOTSPOT_HELP_TEXT, language="text")
                    st.download_button(
                        label="⬇️ hotspot_help.txt 다운로드",
                        data=HOTSPOT_HELP_TEXT,
                        file_name="hotspot_help.txt",
                        mime="text/plain",
                        use_container_width=True,
                        key="dl_hotspot_help",
                    )

            hotspot_raw = st.text_input(
                "핫스팟 잔기",
                value=st.session_state.ifg_hotspot,
                placeholder="예: A70, A73, A166 또는 Ser70, Lys73",
                help="','로 구분하여 입력하세요. 'A70', 'a70', 'Ser70', 'SER 70' 모두 자동으로 'A70' 형태로 변환됩니다.",
                key="input_hotspot",
                on_change=_on_hotspot_change,
                label_visibility="collapsed",
            )
            # 3글자 아미노산 표기('Ser70' 등)까지 포함해 표준형('A70')으로 변환
            hotspot = parse_hotspot_input(hotspot_raw)
            # 변환 전/후가 다를 때만 실시간 변환 결과를 캡션으로 표시
            if hotspot_raw.strip() and hotspot != hotspot_raw.strip():
                st.caption(f"✅ 변환 결과: {hotspot}")

        with col_r:
            num_designs = st.number_input(
                "억제 단백질 후보군 갯수",
                min_value=1,
                value=st.session_state.ifg_num_designs,
                step=1,
                key="input_num_designs",
            )

            iterations = st.selectbox(
                "확산 반복 횟수(Diffusion Iterations)",
                options=[50, 100, 150, 200],
                index=[50, 100, 150, 200].index(st.session_state.ifg_iterations)
                      if st.session_state.ifg_iterations in [50, 100, 150, 200] else 0,
                help="RFdiffusion이 노이즈에서 단백질 구조를 만들어내는 단계 수로, 숫자가 증가할수록 품질은 좋아지지만 속도가 느려짐",
                key="input_iterations",
            )

            hotspot_radius = st.slider(
                "핫스팟 트리밍 반경 (Å)",
                min_value=5.0, max_value=30.0,
                value=float(st.session_state.ifg_hotspot_radius),
                step=0.5,
                key="input_hotspot_radius",
                help="핫스팟 잔기로부터 해당 반경(Å) 이내의 잔기만 사용합니다.",
            )
        
run_disabled = not (saved_pdb_path and Path(saved_pdb_path).exists())



# ─────────────────────────────────────────────────────────────────────────────
# 섹션 C: 내성 유발 단백질 3D 구조 미리보기
# ─────────────────────────────────────────────────────────────────────────────
st.markdown('<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">view_in_ar</span> 내성 유발 단백질 3D 구조 미리보기</h4>', unsafe_allow_html=True)
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
if saved_pdb_path and Path(saved_pdb_path).exists():
    _hs_input = st.session_state.get("input_hotspot", "") or st.session_state.get("ifg_hotspot", "")
    _hs_list  = parse_hotspots(_hs_input)
    render_pdb_viewer(
        saved_pdb_path,
        style="체인별 색상",
        height=600,
        key="preview_viewer",
        hotspots=_hs_list,
    )
    if _hs_list:
        _hs_txt = ", ".join(f'{h["chain"]}{h["resi"]}' for h in _hs_list)
        st.caption(f"🔴 핫스팟 {len(_hs_list)}개 강조 표시: {_hs_txt} | 💡 마우스 드래그: 회전 | 스크롤: 확대·축소")
    else:
        st.caption("💡 마우스 드래그: 회전 | 스크롤: 확대·축소 (핫스팟 잔기를 입력하면 3D에 빨간색으로 표시됩니다)")
else:
    st.markdown(
        '<div style="height:600px;background:#f0f4f8;border-radius:12px;'
        'display:flex;align-items:center;justify-content:center;color:#aaa;'
        'font-size:1.35rem;">PDB 파일을 업로드하면 미리보기가 표시됩니다</div>',
        unsafe_allow_html=True,
    )

# ─────────────────────────────────────────────────────────────────────────────
# 섹션 D: 워크플로우 스텝 카드
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)
st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">rocket</span>InhibiForge 실행</h3>', unsafe_allow_html=True)
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)

_step_cards = [
    ("Step 1", "RFdiffusion",        "내성 억제 단백질 3D 구조 설계", "primary",   "step_1"),
    ("Step 2", "ProteinMPNN",        "내성 억제 단백질 아미노산 서열 생성",          "secondary", "step_2"),
    ("Step 3", "AlphaFold-multimer", "내성 유발 단백질-억제 단백질<br> 결합체 구조 예측",      "primary",   "step_3"),
    ("Step 4", "LightDock",          "내성 유발 단백질-억제 단백질 결합체 구조 평가", "secondary", "step_lightdock"),
]

_STATUS_LABEL = {"waiting": "Waiting", "running": "Running", "done": "Success", "error": "Error"}
_STATUS_COLOR = {"waiting": "#888888", "running": "#005ac1", "done": "#006c4c", "error": "#b02638"}
_icons        = {"waiting": "⏳", "running": "🔄", "done": "✅", "error": "❌"}


def _fmt_elapsed(sec) -> str:
    if sec is None or sec == 0:
        return "--:--:--"
    h, rem = divmod(int(sec), 3600)
    m, s   = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


# 5 컬럼: 4개 단계 카드 + 마지막 "전체 실행" 컬럼
cols = st.columns(5, gap="small")

_status_slots: dict[int, object] = {}
_run_clicks:   dict[int, bool]   = {}

for i, (step, name, desc, color, sk) in enumerate(_step_cards):
    with cols[i]:
        border_c = "#005ac1" if color == "primary" else "#006c4c"
        sts  = st.session_state.get("ifg_pipeline_status", {}).get(sk, "waiting")
        icon = _icons.get(sts, "⏳")
        st.markdown(
            f"""
            <div style="background:#fff;border-left:4px solid {border_c};
                        border-radius:0 12px 12px 0;padding:1rem;
                        box-shadow:0 2px 8px rgba(0,0,0,0.06);height:160px;">
              <div style="font-size:1.08rem;color:{border_c};font-weight:700;
                          letter-spacing:.06em;">{step} {icon}</div>
              <div style="font-family:'Space Grotesk',sans-serif;font-weight:600;
                          font-size:1.425rem;margin:4px 0;">{name}</div>
              <div style="font-size:1.17rem;color:#666;line-height:1.4;">{desc}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
        _run_clicks[i] = st.button(
            "▶ 실행",
            key=f"btn_step_{i}",
            use_container_width=True,
            disabled=run_disabled,
            help="PDB 파일을 먼저 업로드하세요." if run_disabled else f"{name} 단계만 실행합니다.",
        )
        _status_slots[i] = st.empty()

# ── 마지막 컬럼: 전체 실행 ─────────────────────────────────────────────
with cols[4]:
    st.markdown(
        """
        <div style="background:#fff;border-left:4px solid #b02638;
                    border-radius:0 12px 12px 0;padding:1rem;
                    box-shadow:0 2px 8px rgba(0,0,0,0.06);height:160px;">
          <div style="font-size:1.08rem;color:#b02638;font-weight:700;
                      letter-spacing:.06em;">ALL</div>
          <div style="font-family:'Space Grotesk',sans-serif;font-weight:600;
                      font-size:1.425rem;margin:4px 0;">전체 파이프라인</div>
          <div style="font-size:1.17rem;color:#666;line-height:1.4;">
            모든 단계를 순차적으로 실행
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
    run_all_clicked = st.button(
        "▶ 전체 실행",
        key="btn_run_all",
        type="primary",
        use_container_width=True,
        disabled=run_disabled,
        help="PDB 파일을 먼저 업로드하세요." if run_disabled else "모든 단계를 순차 실행합니다.",
    )
    _total_slot = st.empty()


def _render_step_slot(i: int, sk: str) -> None:
    status_map = st.session_state.get("ifg_pipeline_status", {})
    elapsed    = st.session_state.get("ifg_pipeline_elapsed", {})
    end_map    = st.session_state.get("ifg_pipeline_end_times", {})
    sts = status_map.get(sk, "waiting")
    el  = elapsed.get(sk)
    end_ts = end_map.get(sk)
    # "YYYY-MM-DD HH:MM:SS" 에서 시각 부분만 표시
    end_disp = end_ts.split(" ")[-1] if end_ts else "--:--:--"
    lbl = _STATUS_LABEL.get(sts, sts)
    color = _STATUS_COLOR.get(sts, "#888")
    _status_slots[i].markdown(
        f"""
        <div style="padding:0.55rem 0.25rem;font-size:1.17rem;color:#444;
                    line-height:1.55;">
          <div>상태&nbsp;:&nbsp;<b style="color:{color};">{lbl}</b></div>
          <div>소요 시간&nbsp;:&nbsp;<b>{_fmt_elapsed(el)}</b></div>
          <div>종료 시각&nbsp;:&nbsp;<b>{end_disp}</b></div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_total_slot() -> None:
    # 전체 파이프라인 상태는 "전체 실행" 버튼을 눌렀을 때만 갱신된다.
    # 개별 Step 버튼 실행은 전체 상태에 영향을 주지 않는다.
    elapsed = st.session_state.get("ifg_pipeline_elapsed", {})
    overall = st.session_state.get("ifg_overall_status", "waiting")
    total_sec = sum(
        (elapsed.get(sk) or 0)
        for sk in ("step_1", "step_2", "step_3", "step_split", "step_lightdock")
    )
    lbl   = _STATUS_LABEL.get(overall, overall)
    color = _STATUS_COLOR.get(overall, "#888")
    _total_slot.markdown(
        f"""
        <div style="padding:0.55rem 0.25rem;font-size:1.17rem;color:#444;
                    line-height:1.55;">
          <div>상태&nbsp;:&nbsp;<b style="color:{color};">{lbl}</b></div>
          <div>총 소요 시간&nbsp;:&nbsp;<b>{_fmt_elapsed(total_sec if total_sec > 0 else None)}</b></div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_all_status() -> None:
    for i, (_, _, _, _, sk) in enumerate(_step_cards):
        _render_step_slot(i, sk)
    _render_total_slot()


# 실행 버튼을 누르지 않아도 디폴트 상태로 항상 표시
_render_all_status()

st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)
st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)


# ── 파라미터를 세션에 저장 ─────────────────────────────────────────────────
st.session_state.ifg_job_name       = job_name
st.session_state.ifg_contigs        = contigs
st.session_state.ifg_num_designs    = num_designs
st.session_state.ifg_iterations     = iterations
st.session_state.ifg_hotspot        = hotspot
st.session_state.ifg_hotspot_radius = hotspot_radius

st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)

# ── 실행 ──────────────────────────────────────────────────────────────────────
from lib.pipeline_runner import run_rfdiff_step, run_lightdock_step


def _show_error_log() -> None:
    s = st.session_state.get("ifg_pipeline_status", {})
    for step, lbl in [
        ("step_rfdiff",   "RFdiffusion"),
        ("step_1",        "RFdiffusion"),
        ("step_2",        "ProteinMPNN"),
        ("step_3",        "AlphaFold-Multimer"),
        ("step_split",    "체인 분리"),
        ("step_lightdock","LightDock"),
    ]:
        if s.get(step) == "error":
            st.error(f"❌ **{lbl}** 단계에서 오류가 발생했습니다.")

    _log_path   = RFDIFF_OUTPUT_BASE / job_name / "pipeline.log"
    _stderr_msg = st.session_state.get("ifg_pipeline_stderr", "")
    if _log_path.exists() or _stderr_msg:
        with st.expander("🔍 오류 상세 로그 보기"):
            if _log_path.exists():
                _log_lines = _log_path.read_text(encoding="utf-8", errors="replace").splitlines()
                st.code("\n".join(_log_lines[-60:]), language="")
            elif _stderr_msg:
                st.code(_stderr_msg[-3000:], language="")


if not run_disabled:
    cfg = {
        "job_name":        job_name,
        "input_pdb_path":  saved_pdb_path,
        "contigs":         contigs,
        "num_designs":     num_designs,
        "iterations":      iterations,
        "hotspot":         hotspot,
        "hotspot_radius":  hotspot_radius,
    }

    # Step 1~3 (RFdiffusion/ProteinMPNN/AlphaFold-Multimer) 은 동일 백엔드 호출로 묶여 있지만
    # 각 버튼은 자신의 카드 상태 키(step_1/step_2/step_3)에만 영향을 준다.
    # 전체 파이프라인 상태(ifg_overall_status)는 건드리지 않는다.
    _individual_rfdiff = {
        0: ("step_1", "RFdiffusion"),
        1: ("step_2", "ProteinMPNN"),
        2: ("step_3", "AlphaFold-Multimer"),
    }
    _clicked_rfdiff_idx = next(
        (i for i in (0, 1, 2) if _run_clicks.get(i)), None
    )

    if _clicked_rfdiff_idx is not None:
        _sk, _lbl = _individual_rfdiff[_clicked_rfdiff_idx]
        ok = run_rfdiff_step(cfg, status_cb=_render_all_status, step_key=_sk)
        _render_all_status()
        if ok:
            st.success(f"✅ {_lbl} 단계 완료")
        else:
            _show_error_log()

    elif _run_clicks.get(3):
        # LightDock 단독 실행 — 전체 상태는 갱신하지 않는다.
        ok = run_lightdock_step(cfg, status_cb=_render_all_status)
        _render_all_status()
        if ok:
            st.success("✅ LightDock 단계 완료")
        else:
            _show_error_log()

    elif run_all_clicked:
        # 전체 실행: run_pipeline 이 step_1/step_2/step_3 를 stdout 파싱으로
        # 각각 독립 갱신한다. 별도 미러링 불필요.
        for _k in ("step_1", "step_2", "step_3", "step_split", "step_lightdock"):
            st.session_state["ifg_pipeline_status"][_k] = "waiting"
        st.session_state["ifg_pipeline_start_times"] = {}
        st.session_state["ifg_pipeline_end_times"]   = {}
        st.session_state["ifg_pipeline_elapsed"]     = {}
        st.session_state["ifg_overall_status"] = "running"
        _render_all_status()
        ok = run_pipeline(cfg, status_cb=_render_all_status)
        st.session_state["ifg_overall_status"] = "done" if ok else "error"
        _render_all_status()
        if ok:
            st.success("🎉 전체 파이프라인 완료! 결과 보기 탭에서 결과를 확인하세요.")
            st.balloons()
        else:
            _show_error_log()
