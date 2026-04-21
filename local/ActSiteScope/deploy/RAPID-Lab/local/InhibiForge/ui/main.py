"""
InhibiForge_app — main.py
st.navigation 기반 다중 페이지 앱 진입점.

실행: streamlit run /home/sooyeon/amr/LightdockDockQ/InhibiForge_app/main.py
"""
import os
import sys
import copy
from pathlib import Path

# ── 허브 URL (RAPID Lab Hub 으로 돌아가는 홈 버튼용) ────────────────────────
# 도커 통합 실행 시 환경변수로 덮어씀. 기본값은 로컬 WSL에서 hub:8501.
HUB_URL = os.getenv("HUB_URL", "http://localhost:8501")

# ── 경로 부트스트랩 ─────────────────────────────────────────────────────────
# 1) InhibiForge_app/ 를 sys.path에 추가 → `from lib.xxx import ...` 가능
_app_dir = Path(__file__).resolve().parent
if str(_app_dir) not in sys.path:
    sys.path.insert(0, str(_app_dir))

# 2) LightdockDockQ/ 추가 → `from amr_dock_backend import ...` 가능
_ld_base = str(_app_dir.parent)
if _ld_base not in sys.path:
    sys.path.insert(0, _ld_base)

# ── Streamlit import ────────────────────────────────────────────────────────
import streamlit as st
import streamlit.components.v1 as components

# ── 세션 초기화 ─────────────────────────────────────────────────────────────
from lib.session import init_session, _DEFAULTS
init_session()

# ── CSS 로드 ────────────────────────────────────────────────────────────────
from lib.html_loader import load_and_inject_css
_css_path = _app_dir / "styles" / "theme.css"

# ── 페이지 설정 (main.py에서만 호출) ───────────────────────────────────────
st.set_page_config(
    page_title="InhibiForge",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

load_and_inject_css(_css_path)

# ── 사이드바 로고 (HTML 원본과 동일: Re-Off + 알파폴드 팀) ────────────────
st.markdown(
    """
    <div class="sidebar-logo">
      <div class="logo-title">InhibiForge</div>
      <div class="logo-desc">RFdiffusion 및 ProteinMPNN 기반<br/>맞춤형 내성 억제제 설계 및 검증</div>
      <div class="logo-sub">By ResistBreakers</div>
    </div>
    """,
    unsafe_allow_html=True,
)

# 사이드바 너비를 --sidebar-width CSS 변수에 실시간 동기화
components.html(
    """
    <script>
    (function() {
        var sidebar = window.parent.document.querySelector('[data-testid="stSidebar"]');
        if (!sidebar) return;
        var root = window.parent.document.documentElement;
        function sync() {
            var w = sidebar.getBoundingClientRect().width;
            if (w > 50) root.style.setProperty('--sidebar-width', w + 'px');
        }
        sync();
        new ResizeObserver(sync).observe(sidebar);
    })();
    </script>
    """,
    height=0,
)

# ── 사이드바 하단 "새 프로젝트" 버튼 ────────────────────────────────────────
with st.sidebar:
    if st.button("새 프로젝트", key="btn_new_project", use_container_width=True):
        for key, val in _DEFAULTS.items():
            st.session_state[key] = copy.deepcopy(val)
        st.switch_page("pages/1_pipeline.py")

    # --- Home Return Button → RAPID Lab Hub ---
    st.markdown('<div style="margin-top:2rem"></div>', unsafe_allow_html=True)
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"] .st-key-back_to_main_pselector button,
        [data-testid="stSidebar"] .stButton.st-key-back_to_main_pselector > button,
        [data-testid="stSidebar"] .st-key-back_to_main_pselector a,
        [data-testid="stSidebar"] .stLinkButton.st-key-back_to_main_pselector > a {
            background: rgb(239, 68, 68) !important;
            border-color: rgb(239, 68, 68) !important;
            color: #ffffff !important;
            box-shadow: 0 2px 6px rgba(239,68,68,0.25) !important;
        }
        [data-testid="stSidebar"] .st-key-back_to_main_pselector button:hover,
        [data-testid="stSidebar"] .stButton.st-key-back_to_main_pselector > button:hover,
        [data-testid="stSidebar"] .st-key-back_to_main_pselector a:hover,
        [data-testid="stSidebar"] .stLinkButton.st-key-back_to_main_pselector > a:hover {
            background: rgb(220, 38, 38) !important;
            border-color: rgb(220, 38, 38) !important;
            color: #ffffff !important;
            box-shadow: 0 2px 8px rgba(220,38,38,0.35) !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.link_button(
        ":house: 서비스 홈",
        HUB_URL,
        type="primary",
        use_container_width=True,
        key="back_to_main_pselector",
    )

# ── 네비게이션 정의 ─────────────────────────────────────────────────────────
pg_pipeline   = st.Page("pages/1_pipeline.py",        title="InhibiForge 설정 & 실행", icon="⚗️")
pg_validation = st.Page("pages/2_validation_score.py",title="결합 신뢰도 분석 결과",              icon="📈")
pg_lightdock  = st.Page("pages/3_lightdock.py",       title="단백질 도킹 분석 결과",             icon="⚓")
pg_dnaworks   = st.Page("pages/4_dnaworks.py",        title="DNA 설계",              icon="🧬")
pg_storage    = st.Page("pages/5_file_storage.py",    title="파일 보관함",            icon="🗂️")

nav = st.navigation(
    {
        "메인":     [pg_pipeline],
        "결과 보기": [pg_validation, pg_lightdock, pg_dnaworks],
        "보관함":   [pg_storage],
    }
)
nav.run()
