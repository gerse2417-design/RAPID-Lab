"""
pages/4_lightdock.py
LightDock 도킹 결과 페이지 — lightdock_code.html 대응.
"""
import io
import sys
import zipfile
from pathlib import Path
import pandas as pd
import streamlit as st

_app_dir = Path(__file__).resolve().parents[1]
for _p in (str(_app_dir), str(_app_dir.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from lib.session import init_session, get_ld_output_dir
from lib.html_loader import load_and_inject_css
from lib.viewer import render_pdb_viewer, PLDDT_4TIER
from lib.page_header import render_page_header

init_session()
load_and_inject_css(_app_dir / "styles" / "theme.css")

from amr_dock_backend import (
    parse_rank_list, parse_rank_filtered, collect_all_swarm_clusters,
)

render_page_header("내성 유발 단백질과 억제 단백질의 도킹 분석 결과")

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

/* success / info / error / warning 알림 */
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

/* 체크박스 / 토글 */
[data-testid="stMainBlockContainer"] .stCheckbox label,
[data-testid="stMainBlockContainer"] [data-testid="stToggle"] label {
  font-size: 1.35rem !important;
}

/* custom-metric 카드 (LightDock 양성 포즈 + fnat 요약 공용) */
.custom-metric {
  background: var(--secondary-background-color, #f0f2f6);
  border-radius: 0.5rem;
  padding: 1rem;
  text-align: center;
  position: relative;
}
.custom-metric .cm-label {
  font-size: 1.3125rem !important;
  color: var(--text-color, #555);
  line-height: 1.4;
}
.custom-metric .cm-value {
  font-size: 2.625rem !important;
  font-weight: 700;
  color: var(--text-color, #1a1c1e);
  margin-top: 0.25rem;
}
.custom-metric .cm-help {
  position: absolute;
  top: 0.5rem;
  right: 0.5rem;
  border-radius: 50%;
  background: rgba(0,0,0,0.08);
  color: var(--text-color, #888);
  text-align: center;
  cursor: default;
  font-size: 1.125rem !important;
  width: 1.8rem !important;
  height: 1.8rem !important;
  line-height: 1.8rem !important;
}
.cm-tooltip {
  display: none;
  position: absolute;
  top: 2.4rem;
  right: 0.5rem;
  background: #333;
  color: #fff;
  padding: 0.5rem 0.75rem;
  border-radius: 0.4rem;
  line-height: 1.4;
  text-align: left;
  z-index: 999;
  box-shadow: 0 2px 8px rgba(0,0,0,0.2);
  font-size: 1.125rem !important;
  width: 330px !important;
}
.cm-help:hover + .cm-tooltip {
  display: block;
}
</style>
""", unsafe_allow_html=True)

st.markdown("<div style='margin-top:1.5rem;'></div>", unsafe_allow_html=True)

ld_dir = get_ld_output_dir()
rank_path       = ld_dir / "rank_by_scoring.list"
fnat_path       = ld_dir / "filtered" / "rank_filtered.list"
restraints_path = ld_dir / "restraints.list"

has_data = rank_path.exists()
df_rank = pd.DataFrame()
if has_data:
    try:
        df_rank = parse_rank_list(rank_path).sort_values("Scoring", ascending=False).reset_index(drop=True)
    except Exception as e:
        st.error(f"랭킹 파일 파싱 실패: {e}")
        has_data = False

# ─────────────────────────────────────────────────────────────────────────────
# 지표 설명 테이블 (정적 — fragment 밖)
# ─────────────────────────────────────────────────────────────────────────────
_TH = ("background:#b3d9f0;color:#1a1c1e;font-weight:600;padding:0.5rem 0.75rem;"
       "text-align:left;border-bottom:1px solid #8ec8e8;")
_TD = "padding:0.45rem 0.75rem;border-bottom:1px solid #d0e8f5;"
st.markdown(f"""
<div style="background:#e8f4fd;border-radius:12px;padding:1.25rem 1.5rem;margin-bottom:1rem;">
  <div style="font-family:'Space Grotesk',sans-serif;font-size:1.5rem;font-weight:700;
              color:#005ac1;margin-bottom:0.75rem;">단백질간 도킹 평가지표 설명 및 합격 기준</div>
  <table style="width:100%;border-collapse:collapse;font-size:1.275rem;">
    <thead>
      <tr>
        <th style="{_TH}">단백질간 도킹 평가지표</th>
        <th style="{_TH}">의미</th>
        <th style="{_TH}text-align:center;">양호 (Pass)</th>
        <th style="{_TH}text-align:center;">우수 (High)</th>
        <th style="{_TH}">비고</th>
      </tr>
    </thead>
    <tbody>
      <tr style="background:white;"><td style="{_TD}"><b>LightDock Score</b></td><td style="{_TD}">GSO 기반 도킹 점수 (높을수록 안정적)</td><td style="{_TD}text-align:center;">&gt; 0</td><td style="{_TD}text-align:center;">—</td><td style="{_TD}">결합 자유 에너지 <br> 기반 점수</td></tr>
      <tr style="background:white;"><td style="{_TD}"><b>LightDock 양성 점수 <br> 포즈 수</b></td><td style="{_TD}">LightDock Scoring &gt; 0인 유효 포즈 개수</td><td style="{_TD}text-align:center;">≥ 1</td><td style="{_TD}text-align:center;">—</td><td style="{_TD}">유효 결합 포즈 존재 여부</td></tr>
      <tr style="background:white;border-bottom:none;"><td style="{_TD}border-bottom:none;"><b>핫스팟 결합 정확도(fnat)</b></td><td style="{_TD}border-bottom:none;">설계된 억제 단백질이 내성 유발 단백질의 핵심 부위(핫스팟)에 결합한 비율 <br> (0–1)</td><td style="{_TD}border-bottom:none;text-align:center;">≥ 0.5</td><td style="{_TD}border-bottom:none;text-align:center;">≥ 0.75</td><td style="{_TD}border-bottom:none;">지정 결합 잔기 정확도</td></tr>
    </tbody>
  </table>
</div>
""", unsafe_allow_html=True)
st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# 상단: 3D 뷰어 (Top-1 포즈) + 에너지 요약
# ─────────────────────────────────────────────────────────────────────────────
with st.container(border=True):
    col_3d, col_info = st.columns([7, 5], gap="large")

    with col_3d:
        st.markdown('<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">biotech</span>내성 유발 단백질-억제제 결합체(Complex) 3D 구조 뷰어</h4>', unsafe_allow_html=True)
        st.caption("리셉터: 파란색 Cartoon · 리간드: 카툰 · 옵션으로 pLDDT 색상 및 핫스팟 Sphere 표시 가능")

        if has_data and not df_rank.empty:
            pose_options = [
                f"Rank {i+1} | Swarm {int(r['Swarm'])} Glowworm {int(r['Glowworm'])} "
                f"(score={r['Scoring']:.3f})"
                for i, (_, r) in enumerate(df_rank.head(5).iterrows())
            ]
            sel_idx = st.selectbox("내성 유발 단백질-억제제 결합체 선택", range(len(pose_options)),
                                   format_func=lambda i: pose_options[i], key="ld_pose_sel_main")
            sel_row  = df_rank.iloc[sel_idx]
            pdb_path = ld_dir / f"swarm_{int(sel_row['Swarm'])}" / sel_row["PDB"]

            if pdb_path.exists():
                # ── 뷰어 표시 옵션 토글 ────────────────────────────────────
                col_t1, col_t2 = st.columns(2)
                with col_t1:
                    show_plddt = st.toggle(
                        "pLDDT 신뢰도 색상",
                        value=False,
                        key="ld_toggle_plddt",
                        help="켜면 B-factor(pLDDT) 기준 4단계 색상으로 표시됩니다.",
                    )
                with col_t2:
                    show_hotspot = st.toggle(
                        "핫스팟 잔기 표시 (빨간 Sphere)",
                        value=False,
                        key="ld_toggle_hotspot",
                        help="파이프라인 설정의 핫스팟 잔기를 빨간색 Sphere로 강조 표시합니다.",
                    )

                viewer_style = "pLDDT 신뢰도" if show_plddt else "리셉터 파란색 + 리간드 카툰"
                hotspots_val = st.session_state.get("ifg_hotspot", "") if show_hotspot else None

                render_pdb_viewer(
                    str(pdb_path),
                    style=viewer_style,
                    height=452,
                    key="ld_main_viewer",
                    hotspots=hotspots_val,
                )

                # pLDDT 4단계 범례 — pLDDT 색상 토글이 켜져 있을 때만 표시
                if show_plddt:
                    _sw = ("display:inline-block;width:14px;height:14px;"
                           "border-radius:3px;vertical-align:middle;margin-right:6px;")
                    _items_html = "".join(
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
                            pLDDT 신뢰도 분류
                            <span style="font-weight:400;font-size:1.125rem;color:#888;margin-left:6px;">pLDDT 신뢰도 분류 (B-factor 기준)</span>
                          </div>
                          <div style="display:grid;grid-template-columns:repeat(2,1fr);gap:0.4rem 1rem;font-size:1.2rem;color:#333;">
                            {_items_html}
                          </div>
                        </div>
                        """,
                        unsafe_allow_html=True,
                    )
                    st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)
            else:
                st.warning(f"PDB 파일 없음: `{pdb_path}`")
        else:
            st.markdown(
                """
                <div style="height:452px;background:#f0f4f8;border-radius:12px;
                            display:flex;align-items:center;justify-content:center;
                            color:#9aa0ac;font-size:1.425rem;">
                  파이프라인 실행 후 Top 포즈가 표시됩니다
                </div>
                """,
                unsafe_allow_html=True,
            )

    with col_info:
        st.markdown('<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">bar_chart</span>결합 신뢰도 및 점수 요약</h4>', unsafe_allow_html=True)
        
        if has_data and not df_rank.empty:
            # 1. LightDock Score 파트 (2열 배치로 깔끔하게)
            c1, c2 = st.columns(2)
            top1 = df_rank.iloc[0]
            
            with c1:
                st.metric(
                    label="최고 LightDock 점수",
                    value=f"{float(top1['Scoring']):.3f}",
                    help="값이 클수록 결합이 안정적임을 의미합니다."
                )
            with c2:
                n_passed = len(df_rank[df_rank["Scoring"] > 0]) if "Scoring" in df_rank.columns else 0
                c2.markdown(
                    f'<div class="custom-metric">'
                    f'<span class="cm-help">?</span>'
                    f'<span class="cm-tooltip">포즈(Pose)란 억제 단백질이 타겟 단백질 표면에 도킹했을 때 형성되는 3차원적인 결합 위치와 각도(자세)를 의미합니다. LightDock 점수가 양성이면 해당 결합 모델이 물리화학적으로 유효함을 나타냅니다.</span>'
                    f'<div class="cm-label">LightDock 양성 점수 포즈 수</div>'
                    f'<div class="cm-value">{n_passed}개</div></div>',
                    unsafe_allow_html=True,
                )

            st.divider()  # 두 지표 그룹 사이를 시각적으로 분리

            # 2. fnat 점수 파트 (에러 처리 및 파일 존재 여부 확인 추가)
            if fnat_path.exists():
                try:
                    _df_fnat_raw = parse_rank_filtered(fnat_path)
                    if not _df_fnat_raw.empty and "Scoring" in _df_fnat_raw.columns:
                        df_fnat = _df_fnat_raw.sort_values("Scoring", ascending=False).reset_index(drop=True)
                    else:
                        df_fnat = _df_fnat_raw
                    if not df_fnat.empty:
                        st.markdown("<b style='font-size:1.3125rem;'>💡 핫스팟 결합 정확도(fnat) 요약</b>", unsafe_allow_html=True)

                        _fnat_good = int((df_fnat['fnat'] >= 0.5).sum())
                        _fnat_great = int((df_fnat['fnat'] >= 0.75).sum())

                        c_a, c_b, c_c = st.columns(3)
                        c_a.markdown(
                            f'<div class="custom-metric">'
                            f'<span class="cm-help">?</span>'
                            f'<span class="cm-tooltip">LightDock 양성 필터를 통과한 전체 포즈 수입니다. 지정한 결합 잔기를 조금이라도 만족하는 포즈의 총 개수입니다.</span>'
                            f'<div class="cm-label">LightDock 양성<br>필터 통과</div>'
                            f'<div class="cm-value">{len(df_fnat)}개</div></div>',
                            unsafe_allow_html=True,
                        )
                        c_b.markdown(
                            f'<div class="custom-metric">'
                            f'<span class="cm-help">?</span>'
                            f'<span class="cm-tooltip">fnat ≥ 0.5인 포즈 수입니다. 지정한 결합 잔기의 절반 이상이 실제 인터페이스에 존재하는 합리적 수준의 도킹 결과입니다.</span>'
                            f'<div class="cm-label">핫스팟 결합<br>정확도 양호 <br> (≥ 0.5)</div>'
                            f'<div class="cm-value">{_fnat_good}개</div></div>',
                            unsafe_allow_html=True,
                        )
                        c_c.markdown(
                            f'<div class="custom-metric">'
                            f'<span class="cm-help">?</span>'
                            f'<span class="cm-tooltip">fnat ≥ 0.75인 포즈 수입니다. 지정한 결합 잔기의 75% 이상이 인터페이스에 존재하는 고품질 도킹 결과입니다.</span>'
                            f'<div class="cm-label">핫스팟 결합<br>정확도 우수 <br> (≥ 0.75)</div>'
                            f'<div class="cm-value">{_fnat_great}개</div></div>',
                            unsafe_allow_html=True,
                        )
                    else:
                        st.info("필터를 통과한 fnat 데이터가 없습니다.")
                except Exception as e:
                    st.error(f"fnat 데이터 처리 중 오류 발생: {e}")
            else:
                st.info("fnat 필터링 결과 파일이 아직 생성되지 않았습니다.")
                
        else:
            st.info("파이프라인 실행 후 점수 요약이 표시됩니다.")

        # PDB 번들 다운로드
        top5_entries: list[tuple[int, Path, float, int, int]] = []  # (rank, path, score, swarm, glowworm)
        if has_data and not df_rank.empty:
            for rank_idx, (_, row) in enumerate(df_rank.head(5).iterrows(), start=1):
                swarm_id = int(row["Swarm"])
                glowworm_id = int(row["Glowworm"])
                score_val = float(row["Scoring"])
                p = ld_dir / f"swarm_{swarm_id}" / row["PDB"]
                if p.exists():
                    top5_entries.append((rank_idx, p, score_val, swarm_id, glowworm_id))
        if top5_entries:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for rank_idx, p, score_val, swarm_id, glowworm_id in top5_entries:
                    # ZIP 내부 파일명에 순위/점수/swarm/glowworm 정보 부여
                    arcname = (
                        f"rank{rank_idx:02d}_score{score_val:.3f}_"
                        f"swarm{swarm_id}_glowworm{glowworm_id}.pdb"
                    )
                    zf.write(p, arcname)
            buf.seek(0)
            st.download_button(
                f"📥 Download PDB Bundle (Top 5) · {len(top5_entries)}개 포함",
                data=buf.getvalue(),
                file_name="lightdock_top5.zip",
                mime="application/zip",
                use_container_width=True,
                key="dl_ld_bundle",
            )
        else:
            st.download_button(
                "📥 Download PDB Bundle (Top 5)",
                data=b"",
                file_name="lightdock_top5.zip",
                mime="application/zip",
                use_container_width=True,
                key="dl_ld_bundle",
                disabled=True,
            )

st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# Restraints 만족도 (fnat) 및 도킹 결과 요약
# ─────────────────────────────────────────────────────────────────────────────

# 섹션 타이틀을 데이터 통합에 맞게 살짝 보강했습니다
st.markdown('<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">straighten</span>Top 5 우수 포즈 분석 및 핫스팟 적중률</h4>', unsafe_allow_html=True)
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)

show_all_poses = st.checkbox(
    "🔍 모든 포즈 보기 (fnat 필터 무시)",
    value=False,
    help="체크 시 fnat 계산/필터 통과 여부와 무관하게 rank_by_scoring.list 기반 상위 5개 포즈를 표시합니다.",
    key="ld_show_all_poses",
)

if not has_data:
    st.info("파이프라인 실행 후 fnat 결과가 표시됩니다.")
elif show_all_poses:
    if not df_rank.empty:
        df_all = df_rank.copy()
        df_all.insert(0, "순위", range(1, len(df_all) + 1))
        if "PDB" not in df_all.columns:
            df_all["PDB"] = df_all.apply(
                lambda r: f"swarm_{int(r['Swarm'])}/lightdock_{int(r['Glowworm'])}.pdb",
                axis=1,
            )
        rank_top5 = df_all.head(5)[["순위", "Swarm", "Glowworm", "Scoring", "PDB"]]
        st.dataframe(
            rank_top5,
            use_container_width=True, hide_index=True,
            column_config={
                "Swarm":    st.column_config.NumberColumn("탐색구역 (Swarm)"),
                "Glowworm": st.column_config.NumberColumn("포즈ID (Glowworm)"),
                "Scoring":  st.column_config.NumberColumn("LightDock Score", format="%.3f"),
                "PDB":      st.column_config.TextColumn("PDB 파일 경로", width="medium"),
            },
        )
        st.download_button(
            "📥 Download Top 5 (CSV, fnat 무시)",
            data=rank_top5.to_csv(index=False).encode("utf-8-sig"),
            file_name="lightdock_rank_top5.csv",
            mime="text/csv",
            use_container_width=True,
            key="dl_rank_top5",
        )
    else:
        st.info("rank_by_scoring.list 데이터가 비어 있습니다.")
elif fnat_path.exists():
    try:
        _df_fnat_raw = parse_rank_filtered(fnat_path)
        if not _df_fnat_raw.empty and "Scoring" in _df_fnat_raw.columns:
            df_fnat = _df_fnat_raw.sort_values("Scoring", ascending=False).reset_index(drop=True)
        else:
            df_fnat = _df_fnat_raw
        if not df_fnat.empty:
            
            # 메트릭(요약 수치) 출력 부분 삭제 완료
            
            df_fnat.insert(0, "순위", range(1, len(df_fnat) + 1))
            
            # --- [수정 포인트 1] PDB 파일 경로 동적 생성 및 추가 ---
            if "PDB" not in df_fnat.columns:
                # LightDock 기본 규칙 (swarm_X/lightdock_Y.pdb) 적용
                df_fnat["PDB"] = df_fnat.apply(lambda r: f"swarm_{int(r['Swarm'])}/lightdock_{int(r['Glowworm'])}.pdb", axis=1)

            # --- [수정 포인트 2] 표시할 컬럼에 PDB 추가 ---
            fnat_top5 = df_fnat.head(5)[["순위", "Swarm", "Glowworm", "Scoring", "fnat", "PDB"]]
            
            st.dataframe(
                fnat_top5,
                use_container_width=True, hide_index=True,
                column_config={
                    "Swarm": st.column_config.NumberColumn("탐색구역 (Swarm)"),
                    "Glowworm": st.column_config.NumberColumn("포즈ID (Glowworm)"),
                    "Scoring": st.column_config.NumberColumn("LightDock Score", format="%.3f"),
                    "fnat": st.column_config.ProgressColumn("fnat (정확도)", min_value=0.0, max_value=1.0, format="%.3f"),
                    "PDB": st.column_config.TextColumn("PDB 파일 경로", width="medium"), # 추가된 UI 설정
                },
            )
            
            st.download_button(
                "📥 Download fnat Top 5 (CSV)",
                data=fnat_top5.to_csv(index=False).encode("utf-8-sig"),
                file_name="lightdock_fnat_top5.csv",
                mime="text/csv",
                use_container_width=True,
                key="dl_fnat_top5",
            )
    except Exception as e:
        st.warning(f"fnat 파싱 오류: {e}")

elif restraints_path.exists():
    rank_path_f = ld_dir / "rank_by_scoring.list"
    if rank_path_f.exists():
        if st.button("▶ Restraints 필터 실행 (fnat 계산)", key="btn_run_filter"):
            import shutil, subprocess, sys
            from pathlib import Path as _P

            filtered_dir = ld_dir / "filtered"
            if filtered_dir.exists():
                shutil.rmtree(filtered_dir)

            _filter_script = _P(sys.executable).parent / "lgd_filter_restraints.py"
            with st.spinner("fnat 계산 중…"):
                proc = subprocess.run(
                    [sys.executable, str(_filter_script),
                     "rank_by_scoring.list", "restraints.list", "A", "B", "--fnat", "0.0"],
                    cwd=str(ld_dir),
                    capture_output=True,
                )
            if proc.returncode == 0 and (ld_dir / "filtered" / "rank_filtered.list").exists():
                st.success("✅ fnat 필터링 완료!")
                st.rerun()
            else:
                st.error(f"❌ 필터 실행 실패 (rc={proc.returncode})")
                if proc.stderr:
                    st.code(proc.stderr.decode(errors="replace")[-1000:], language="text")
    else:
        st.info("rank_by_scoring.list가 없습니다. 파이프라인을 먼저 실행하세요.")

st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)

