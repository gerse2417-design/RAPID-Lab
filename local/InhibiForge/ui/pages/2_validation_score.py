"""
pages/3_validation_score.py
AlphaFold-Multimer 검증 점수 페이지 — validationscore_code.html 대응.
"""
import io
import sys
import zipfile
from pathlib import Path
import streamlit as st

_app_dir = Path(__file__).resolve().parents[1]
for _p in (str(_app_dir), str(_app_dir.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from lib.session import init_session, get_rf_output_dir, get_ld_output_dir
from lib.html_loader import load_and_inject_css
from lib.page_header import render_page_header
from lib.viewer import get_pdb_stats, render_pdb_viewer, VIEWER_STYLES, PLDDT_4TIER

init_session()
load_and_inject_css(_app_dir / "styles" / "theme.css")

import pandas as pd

render_page_header("내성 유발 단백질과 억제 단백질의 결합 신뢰도 분석 결과")

# ── 이 페이지 글씨 1.5배 확대 CSS ────────────────────────────────────────────
st.markdown("""
<style>
/* 페이지 헤더 */
.ifg-page-title    { font-size: 2.25em !important; }
.ifg-page-subtitle { font-size: 2.25em !important; }

/* Streamlit 위젯 라벨 및 본문 텍스트 */
[data-testid="stMainBlockContainer"] label,
[data-testid="stMainBlockContainer"] .stMarkdown p,
[data-testid="stMainBlockContainer"] .stMarkdown li,
[data-testid="stMainBlockContainer"] [data-testid="stWidgetLabel"] p {
  font-size: 1.3125rem !important;
}

/* h3 / h4 태그 */
[data-testid="stMainBlockContainer"] h3 {
  font-size: 2.1rem !important;
}
[data-testid="stMainBlockContainer"] h4 {
  font-size: 1.8rem !important;
}

/* 입력 위젯 내부 텍스트 */
[data-testid="stMainBlockContainer"] .stTextInput input,
[data-testid="stMainBlockContainer"] .stNumberInput input,
[data-testid="stMainBlockContainer"] .stTextArea textarea,
[data-testid="stMainBlockContainer"] .stSelectbox [data-baseweb="select"] {
  font-size: 1.3125rem !important;
}

/* 버튼 텍스트 */
[data-testid="stMainBlockContainer"] .stButton button {
  font-size: 1.3125rem !important;
}

/* caption */
[data-testid="stMainBlockContainer"] [data-testid="stCaptionContainer"] p {
  font-size: 1.2rem !important;
}

/* success / info / error 알림 */
[data-testid="stMainBlockContainer"] [data-testid="stAlert"] p {
  font-size: 1.3125rem !important;
}

/* metric 값 */
[data-testid="stMainBlockContainer"] [data-testid="stMetricValue"] {
  font-size: 2.625rem !important;
}
[data-testid="stMainBlockContainer"] [data-testid="stMetricLabel"] {
  font-size: 1.2rem !important;
}

/* selectbox 옵션 텍스트 */
[data-testid="stMainBlockContainer"] [data-baseweb="select"] span {
  font-size: 1.3125rem !important;
}

/* 다운로드 버튼 */
[data-testid="stMainBlockContainer"] [data-testid="stDownloadButton"] button {
  font-size: 1.3125rem !important;
}

/* 데이터프레임 */
[data-testid="stMainBlockContainer"] [data-testid="stDataFrame"] {
  font-size: 1.2rem !important;
}

/* 체크박스 */
[data-testid="stMainBlockContainer"] .stCheckbox label {
  font-size: 1.35rem !important;
}
</style>
""", unsafe_allow_html=True)

st.markdown("<div style='margin-top:1.5rem;'></div>", unsafe_allow_html=True)

rf_dir   = get_rf_output_dir()
rank_csv = rf_dir / "mpnn_results.csv"
if not rank_csv.exists():
    rank_csv = rf_dir / "rank.csv"  # fallback
if not rank_csv.exists():
    st.info("💡 파이프라인을 먼저 실행하면 검증 점수가 여기에 표시됩니다.")
    st.stop()

# ─────────────────────────────────────────────────────────────────────────────
# 지표 설명 테이블 (정적 — fragment 밖)
# ─────────────────────────────────────────────────────────────────────────────
_TH = ("background:#b3d9f0;color:#1a1c1e;font-weight:600;padding:0.5rem 0.75rem;"
       "text-align:left;border-bottom:1px solid #8ec8e8;")
_TD = "padding:0.45rem 0.75rem;border-bottom:1px solid #d0e8f5;"
st.markdown(f"""
<div style="background:#e8f4fd;border-radius:12px;padding:1.25rem 1.5rem;margin-bottom:1rem;">
  <div style="font-family:'Space Grotesk',sans-serif;font-size:1.5rem;font-weight:700;
              color:#005ac1;margin-bottom:0.75rem;">결합 신뢰도 평가지표 설명 및 합격 기준</div>
  <table style="width:100%;border-collapse:collapse;font-size:1.275rem;">
    <thead>
      <tr>
        <th style="{_TH}">결합 신뢰도 <br>평가지표</th>
        <th style="{_TH}">의미</th>
        <th style="{_TH}text-align:center;">신뢰 가능</th>
        <th style="{_TH}text-align:center;">매우 신뢰 가능</th>
        <th style="{_TH}">비고</th>
      </tr>
    </thead>
    <tbody>
      <tr style="background:white;"><td style="{_TD}"><b>종합 신뢰도 <br> 점수</b></td><td style="{_TD}">단백질 복합체 전체의 전역적(Global) 및<br> 결합면(Interface) 신뢰도</td><td style="{_TD}text-align:center;">0.5 ~ 0.75</td><td style="{_TD}text-align:center;">≥ 0.75</td><td style="{_TD}">0.8 * ipTM + 0.2 * pTM</td></tr>
      <tr style="background:white;"><td style="{_TD}"><b>pLDDT</b></td><td style="{_TD}">각 아미노산 잔기(Residue)<br> 단위의 신뢰도 (0–100)</td><td style="{_TD}text-align:center;">≥ 70</td><td style="{_TD}text-align:center;">≥ 90</td><td style="{_TD}">생성된 억제 단백질(Inhibitor) 자체가 안정적인 구조를 형성하고 있는지 확인</td></tr>
      <tr style="background:white;"><td style="{_TD}"><b>ipTM</b></td><td style="{_TD}">결합 인터페이스(Interface) 부위의 <br> 예측 정확도에 특화된 점수 (0–1)</td><td style="{_TD}text-align:center;">≥ 0.8</td><td style="{_TD}text-align:center;">-</td><td style="{_TD}">억제 단백질이 타겟 단백질의 작용 부위에 얼마나 정확하게 결합했는지를 나타냄<br> • 0.6 미만: 예측 실패 가능성이 높음<br>
 • 0.6~0.8 사이: 예측이 정확할 수도 있고 아닐 수도 있는 불확실한 영역</td></tr>
      <tr style="background:white;border-bottom:none;"><td style="{_TD}border-bottom:none;"><b>pTM</b></td><td style="{_TD}border-bottom:none;">전체적인 구조적 유사성(Global fold) <br> 예측 점수 (0–1)</td><td style="{_TD}border-bottom:none;text-align:center;">≥ 0.5</td><td style="{_TD}border-bottom:none;text-align:center;">-</td><td style="{_TD}border-bottom:none;">복합체 전체의 형태가 얼마나 정확하게 모델링 되었는지 나타냄 <br> • pTM 점수가 0.5 이상이면, 예측된 복합체(complex)의 전반적인 폴딩(fold)이 <br> 실제 구조와 유사할 가능성이 있음.</td></tr>
    </tbody>
  </table>
</div>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# 내성 유발 단백질과 억제 단백질 결합체 3D 구조 뷰어
# ─────────────────────────────────────────────────────────────────────────────
ld_dir   = get_ld_output_dir()
best_pdb = rf_dir / "best.pdb"

st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)
st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">view_in_ar</span>내성 유발 단백질과 억제 단백질 결합체 3D 구조 뷰어</h3>', unsafe_allow_html=True)
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)

# ── PDB 파일 목록 수집 ─────────────────────────────────────────────────────
pdb_files: dict[str, Path] = {}

if best_pdb.exists():
    pdb_files["🏆 best.pdb (최고 점수 바인더)"] = best_pdb

for p in sorted(rf_dir.glob("*.pdb")):
    if p.name != "best.pdb":
        pdb_files[f"🔷 {p.name} (RFdiffusion 백본)"] = p

# LightDock 포즈 상위 5개
rank_path = ld_dir / "rank_by_scoring.list"
if rank_path.exists():
    try:
        sys.path.insert(0, str(_app_dir.parent))
        from amr_dock_backend import parse_rank_list
        df_ld_top = parse_rank_list(rank_path).sort_values("Scoring", ascending=False).head(5)
        for _, row in df_ld_top.iterrows():
            pose_p = ld_dir / f"swarm_{int(row['Swarm'])}" / row["PDB"]
            if pose_p.exists():
                pdb_files[f"⚓ {row['PDB']} (LightDock 포즈, 점수={row['Scoring']:.3f})"] = pose_p
    except Exception:
        pass

if pdb_files:
    # 컨트롤 패널
    st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
    st.markdown("""
    <style>
    div[data-testid="stVerticalBlockBorderWrapper"],
    div[data-testid="stVerticalBlockBorderWrapper"] > div[data-testid="stVerticalBlock"],
    div[data-testid="stVerticalBlock"][style*="border"] {
        background-color: #ffffff !important;
    }
    </style>
    """, unsafe_allow_html=True)
    with st.container(border=True):
        ctrl1, ctrl3 = st.columns([7, 1])
        with ctrl1:
            sel_label = st.selectbox("PDB 파일 선택", list(pdb_files.keys()), key="pdb_sel_3d")
        with ctrl3:
            spin = st.checkbox("자동 회전", value=False, key="spin_3d")
        viewer_style = "pLDDT 신뢰도"

    selected_pdb = pdb_files[sel_label]

    # 체인 하이라이트 session state 초기화
    if "highlight_chain" not in st.session_state:
        st.session_state.highlight_chain = None
    # PDB 파일 변경 시 하이라이트 초기화
    if st.session_state.get("_prev_pdb_sel") != sel_label:
        st.session_state.highlight_chain = None
        st.session_state._prev_pdb_sel = sel_label

    # 메인 뷰어 + 우측 패널
    st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
    with st.container(border=True):
        col_main, col_side = st.columns([9, 3], gap="medium")

        with col_main:
            render_pdb_viewer(
                selected_pdb, style=viewer_style, spin=spin, height=520,
                key="main_3d_viewer",
                highlight_chain=st.session_state.highlight_chain,
            )

            # pLDDT 신뢰도 범례 (4단계)
            if viewer_style == "pLDDT 신뢰도":
                _sw = ("display:inline-block;width:14px;height:14px;"
                       "border-radius:3px;vertical-align:middle;margin-right:6px;")
                items_html = "".join(
                    f'<span style="display:flex;align-items:center;gap:6px;white-space:nowrap;">'
                    f'<span style="{_sw}background:{c};"></span>'
                    f'<span style="font-weight:600;">{label}</span>'
                    f'<span style="color:#888;">({rng})</span>'
                    f'</span>'
                    for c, label, rng in PLDDT_4TIER
                )
                st.markdown(
                    f"""
                    <div style="background:#ffffff;border:1px solid #e0e4ea;border-radius:10px;
                                padding:0.75rem 1rem;margin-top:0.5rem;">
                      <div style="font-weight:600;font-size:1.275rem;color:#1a1c1e;margin-bottom:0.6rem;">
                        각 아미노산 잔기 단위 신뢰도 분류(pLDDT)
                      </div>
                      <div style="display:grid;grid-template-columns:repeat(2,1fr);gap:0.4rem 1rem;font-size:1.2rem;color:#333;">
                        {items_html}
                      </div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
            else:
                # 일반 조작 안내
                st.markdown(
                    """
                    <div style="display:flex;flex-wrap:wrap;gap:0.5rem 1.5rem;
                                padding:0.75rem 1rem;background:#ffffff;
                                border-radius:10px;font-size:1.17rem;color:#555;
                                margin-top:0.5rem;width:100%;box-sizing:border-box;">
                      <span>🖱️ 왼쪽 드래그: 회전</span>
                      <span>🖱️ 스크롤: 확대·축소</span>
                      <span>🖱️ 오른쪽 드래그: 이동</span>
                      <span>🖱️ 더블클릭: 초점</span>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

        with col_side:
            with st.container(border=True):
                st.markdown(
                    "<style>"
                    "div[data-testid='stVerticalBlockBorderWrapper']:has(h4) {"
                    "  background-color: #ffffff !important;"
                    "}"
                    "div[data-testid='column'] div[data-testid='stVerticalBlockBorderWrapper'] {"
                    "  background-color: #ffffff !important;"
                    "}"
                    "</style>",
                    unsafe_allow_html=True,
                )
                st.markdown('<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">bar_chart</span> 구조 정보</h4>', unsafe_allow_html=True)
                st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
                stats_3d = get_pdb_stats(selected_pdb)

                # 파일 메타
                import os, datetime
                mtime = os.path.getmtime(selected_pdb)
                st.markdown(
                    f"""
                    <div style="font-size:1.2rem;color:#555;margin-bottom:1rem;">
                      <div>📁 <b>{selected_pdb.name}</b></div>
                      <div>🗓️ {datetime.datetime.fromtimestamp(mtime).strftime('%Y.%m.%d')}</div>
                      <div>⚖️ {os.path.getsize(selected_pdb) / 1024:.1f} KB</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

                st.metric("총 잔기 수", stats_3d.get("total_res", "—"))
                st.metric("총 원자 수", stats_3d.get("total_atoms", "—"))

                # 체인별 정보 (클릭하면 3D 뷰어에서 해당 체인 하이라이트)
                if stats_3d.get("chains"):
                    st.markdown("**체인 분석** <span style='font-size:0.85rem;color:#888;'>(클릭하여 하이라이트)</span>", unsafe_allow_html=True)

                    for ch in stats_3d["chains"]:
                        plddt_pct = min(100, ch["avg_plddt"])
                        color     = "#006c4c" if plddt_pct >= 70 else ("#e88000" if plddt_pct >= 50 else "#b02638")
                        is_active = st.session_state.highlight_chain == ch["id"]
                        border_style = f"border:2px solid {color};" if is_active else "border:2px solid transparent;"
                        bg_style = "background:#f0f7ff;" if is_active else "background:transparent;"

                        if st.button(
                            f"Chain {ch['id']}  ·  {ch['n_res']} 잔기  ·  pLDDT {ch['avg_plddt']:.0f}",
                            key=f"chain_btn_{ch['id']}",
                            use_container_width=True,
                        ):
                            # 토글: 같은 체인 다시 클릭하면 하이라이트 해제
                            if st.session_state.highlight_chain == ch["id"]:
                                st.session_state.highlight_chain = None
                            else:
                                st.session_state.highlight_chain = ch["id"]
                            st.rerun()

                        # 프로그레스 바
                        st.markdown(
                            f"""
                            <div style="height:6px;background:#eee;border-radius:3px;margin:-8px 0 8px;">
                              <div style="height:100%;width:{plddt_pct}%;background:{color};border-radius:3px;"></div>
                            </div>
                            """,
                            unsafe_allow_html=True,
                        )

                # 다운로드
                st.download_button(
                    "📥 PDB 다운로드",
                    data=selected_pdb.read_bytes(),
                    file_name=selected_pdb.name,
                    mime="chemical/x-pdb",
                    use_container_width=True,
                    key="dl_3d_pdb",
                )

        st.markdown("<div style='padding-bottom:3rem;'></div>", unsafe_allow_html=True)
else:
    st.info("💡 파이프라인을 먼저 실행하면 PDB 파일을 여기서 확인할 수 있습니다.")

st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# CSV 로드 (fragment 밖 — 한 번만 읽음)
# ─────────────────────────────────────────────────────────────────────────────
try:
    _df_raw = pd.read_csv(rank_csv)
except Exception as e:
    st.error(f"CSV 읽기 실패: {e}")
    st.stop()

if "design" in _df_raw.columns:
    _df_raw = _df_raw.rename(columns={"design": "pdb_file"})
elif "pdb_file" not in _df_raw.columns:
    _df_raw.insert(0, "pdb_file", "unknown")

# pdb_file 값을 파일명만 남기고 .pdb 확장자 보장
_df_raw["pdb_file"] = _df_raw["pdb_file"].astype(str).apply(
    lambda x: Path(x).stem + ".pdb" if not x.endswith(".pdb") else Path(x).name
)

if "plddt" in _df_raw.columns:
    # pLDDT 값이 0~1 범위인 경우 100을 곱하여 0~100 범위로 변환
    if _df_raw["plddt"].max() <= 1.0:
        _df_raw["plddt"] = _df_raw["plddt"] * 100
    _df_raw = _df_raw.sort_values("plddt", ascending=False).reset_index(drop=True)
_df_raw.insert(0, "RANK", range(1, len(_df_raw) + 1))


# ─────────────────────────────────────────────────────────────────────────────
# Fragment: 필터 + 결과표 + 차트 + 인사이트 + 다운로드
# (필터 적용 시 이 영역만 재실행됨)
# ─────────────────────────────────────────────────────────────────────────────
@st.fragment
def _results_panel():
    df = _df_raw.copy()

    # ── 필터 컨트롤 ──────────────────────────────────────────────────────────
    st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">filter_list</span>내성 유발 단백질과 억제 단백질 결합 신뢰도 검증 결과</h3>', unsafe_allow_html=True)
    st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
    st.markdown('<div style="font-family:\'Space Grotesk\',sans-serif;font-size:1.575rem;font-weight:700;color:var(--on-surface);margin-bottom:0.75rem;">검증 필터</div>', unsafe_allow_html=True)
    with st.container(border=True):
        fc1, fc2, fc3, fc4, fc5 = st.columns([2, 2, 2, 2, 1])

        thr_overall = fc1.number_input(
            "종합 신뢰도 점수 최솟값", min_value=0.0, max_value=1.0, step=0.05,
            value=float(st.session_state.ifg_thr_overall), key="thr_overall_in",
        )
        thr_plddt = fc2.number_input(
            "pLDDT 최솟값", min_value=0, max_value=100, step=1,
            value=st.session_state.ifg_thr_plddt, key="thr_plddt_in",
        )
        thr_iptm = fc3.number_input(
            "ipTM 최솟값", min_value=0.0, max_value=1.0, step=0.05,
            value=float(st.session_state.ifg_thr_iptm), key="thr_iptm_in",
        )
        thr_ptm = fc4.number_input(
            "pTM 최솟값", min_value=0.0, max_value=1.0, step=0.05,
            value=float(st.session_state.ifg_thr_ptm), key="thr_ptm_in",
        )
        fc5.markdown("<div style='padding-top:1.7rem'>", unsafe_allow_html=True)
        if fc5.button("✅ 적용", key="apply_filter", use_container_width=True):
            st.session_state.ifg_thr_overall = thr_overall
            st.session_state.ifg_thr_plddt   = thr_plddt
            st.session_state.ifg_thr_iptm    = thr_iptm
            st.session_state.ifg_thr_ptm     = thr_ptm
        fc5.markdown("</div>", unsafe_allow_html=True)

        eff_overall = st.session_state.ifg_thr_overall
        eff_plddt   = st.session_state.ifg_thr_plddt
        eff_iptm    = st.session_state.ifg_thr_iptm
        eff_ptm     = st.session_state.ifg_thr_ptm
        st.caption(f"적용 중: 종합 ≥ **{eff_overall}** | pLDDT ≥ **{eff_plddt}** | i_ptm ≥ **{eff_iptm}** | ptm ≥ **{eff_ptm}**")

    # ── 종합 신뢰도 점수 계산 (0.8 * ipTM + 0.2 * pTM) ──────────────────────
    if "i_ptm" in df.columns and "ptm" in df.columns:
        df["overall_score"] = 0.8 * df["i_ptm"] + 0.2 * df["ptm"]

    # ── STATUS 계산 ──────────────────────────────────────────────────────────
    def _status(row):
        checks = []
        if "overall_score" in row.index: checks.append(row["overall_score"] >= eff_overall)
        if "plddt" in row.index: checks.append(row["plddt"] >= eff_plddt)
        if "i_ptm" in row.index: checks.append(row["i_ptm"] >= eff_iptm)
        elif "ptm" in row.index: checks.append(row["ptm"] >= eff_iptm)
        if "ptm" in row.index:   checks.append(row["ptm"] >= eff_ptm)
        if not checks:
            return "PASSED"
        if all(checks):
            return "PASSED"
        if sum(checks) >= len(checks) - 1:
            return "MARGINAL"
        return "FAILED"

    df["STATUS"] = df.apply(_status, axis=1)

    numeric_cols = [c for c in df.columns if c in ("overall_score", "plddt", "i_ptm", "ptm")]
    col_cfg = {
        "RANK":          st.column_config.NumberColumn("순위", width=5),
        "STATUS":        st.column_config.TextColumn("검증 필터 통과 여부"),
        "pdb_file":      st.column_config.TextColumn("PDB 파일"),
        "overall_score": st.column_config.NumberColumn("종합 신뢰도 점수", format="%.3f"),
        "plddt":         st.column_config.NumberColumn("pLDDT", format="%.3f"),
        "i_ptm":         st.column_config.NumberColumn("인터페이스 pTM(i_pTM)", format="%.3f"),
        "ptm":           st.column_config.NumberColumn("전체 pTM(ptm)", format="%.3f"),
    }

    # ── 결과 테이블 ──────────────────────────────────────────────────────────
    st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
    st.markdown('<div style="font-family:\'Space Grotesk\',sans-serif;font-size:1.575rem;font-weight:700;color:var(--on-surface);margin-bottom:0.75rem;">검증 결과</div>', unsafe_allow_html=True)
    show_cols = ["RANK", "STATUS", "overall_score", "plddt", "i_ptm", "ptm", "pdb_file"]
    show_cols = [c for c in show_cols if c in df.columns]
    # 상위 5개만 표시 (종합 신뢰도 점수 / pLDDT 기준 기본 정렬)
    display_df = df[show_cols].head(5)

    # 필터 기준에 연동된 셀 색상: 통과=초록, 미달=빨강, 배경=흰색
    _PASS = "background-color: #d4edda"  # 연한 초록
    _FAIL = "background-color: #f8d7da"  # 연한 빨강
    _DFLT = "background-color: white"

    def _highlight_by_filter(row):
        styles = [_DFLT] * len(row)
        for i, col in enumerate(row.index):
            if col == "overall_score" and col in display_df.columns:
                styles[i] = _PASS if row[col] >= eff_overall else _FAIL
            elif col == "plddt" and col in display_df.columns:
                styles[i] = _PASS if row[col] >= eff_plddt else _FAIL
            elif col == "i_ptm" and col in display_df.columns:
                styles[i] = _PASS if row[col] >= eff_iptm else _FAIL
            elif col == "ptm" and col in display_df.columns:
                styles[i] = _PASS if row[col] >= eff_ptm else _FAIL
        return styles

    styled = display_df.style.apply(_highlight_by_filter, axis=1)
    if numeric_cols:
        styled = styled.format("{:.3f}", subset=numeric_cols)

    st.dataframe(styled, use_container_width=True, hide_index=True, column_config=col_cfg)

    # ── Quick Insight ────────────────────────────────────────────────────────
    n_pass  = int((df["STATUS"] == "PASSED").sum())
    n_total = len(df)
    pct = n_pass / n_total * 100 if n_total else 0
    st.markdown(
        f"""
        <div style="background:linear-gradient(135deg,#f0f4f8,#e8f0fd);
                    border-radius:14px;padding:1.25rem 1.5rem;margin-top:1rem;">
          <div style="font-family:'Space Grotesk',sans-serif;font-size:1.65rem;
                      font-weight:700;color:#1a1a2e;">💡 검증 결과 해석</div>
          <div style="font-size:1.35rem;color:#444;margin-top:0.5rem;">
            데이터 세트의 <b>{pct:.0f}%</b>가 높은 신뢰성을 보이고 있습니다.
            ({n_pass} / {n_total} 디자인 통과)
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ── 다운로드 ─────────────────────────────────────────────────────────────
    st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
    st.markdown('<div style="font-family:\'Space Grotesk\',sans-serif;font-size:1.575rem;font-weight:700;color:var(--on-surface);margin-bottom:0.75rem;">검증 결과 표 및 Top 5 우수 복합체 다운로드</div>', unsafe_allow_html=True)
    st.markdown("<div style='margin-top:0rem;'></div>", unsafe_allow_html=True)
    dl_col1, dl_col2 = st.columns(2)

    with dl_col1:
        # 화면에 보이는 Top 5 검증결과표 그대로 CSV 로 다운로드
        top5_csv_bytes = display_df.to_csv(index=False).encode("utf-8-sig")
        st.download_button(
            "📥 Top 5 검증결과표 다운로드 (CSV)",
            data=top5_csv_bytes,
            file_name="validation_top5.csv",
            mime="text/csv",
            use_container_width=True,
        )

    with dl_col2:
        if "plddt" in df.columns and "pdb_file" in df.columns:
            top5 = df.nlargest(5, "plddt")
            pdb_zip_buf = io.BytesIO()
            all_pdb_dir = rf_dir / "all_pdb"
            written = 0
            with zipfile.ZipFile(pdb_zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for rank_idx, (_, row) in enumerate(top5.iterrows(), start=1):
                    pdb_fname = str(row["pdb_file"]).replace("🏆 ", "").strip()
                    design_val = Path(pdb_fname).stem if pdb_fname else ""
                    n_val = row.get("n", None)

                    candidates: list[Path] = []
                    # 1) CSV의 pdb_file 값이 실제 파일명인 경우
                    if pdb_fname:
                        candidates += [
                            rf_dir / pdb_fname,
                            rf_dir / (pdb_fname + ".pdb"),
                        ]
                        if all_pdb_dir.exists():
                            candidates.append(all_pdb_dir / pdb_fname)
                    # 2) design/n 인덱스로 실제 파일 경로 재구성 (all_pdb/designX_nY.pdb)
                    if design_val != "" and pd.notna(n_val):
                        try:
                            pattern_name = f"design{int(float(design_val))}_n{int(n_val)}.pdb"
                            if all_pdb_dir.exists():
                                candidates.append(all_pdb_dir / pattern_name)
                            candidates.append(rf_dir / pattern_name)
                        except (TypeError, ValueError):
                            pass

                    for p in candidates:
                        if p.exists():
                            # ZIP 내부 파일명에 pLDDT 순위 prefix 부여
                            plddt_val = row.get("plddt", None)
                            if pd.notna(plddt_val):
                                arcname = f"rank{rank_idx:02d}_plddt{float(plddt_val):.1f}_{p.name}"
                            else:
                                arcname = f"rank{rank_idx:02d}_{p.name}"
                            zf.write(p, arcname)
                            written += 1
                            break
            pdb_zip_buf.seek(0)

            st.download_button(
                f"📥 Top 5 PDB 다운로드 (ZIP) · {written}개 포함",
                data=pdb_zip_buf.getvalue(),
                file_name="top5_pdbs.zip",
                mime="application/zip",
                use_container_width=True,
                disabled=(written == 0),
            )
            if written == 0:
                st.warning("⚠️ 상위 5개 디자인에 해당하는 PDB 파일을 찾을 수 없습니다.")


_results_panel()
