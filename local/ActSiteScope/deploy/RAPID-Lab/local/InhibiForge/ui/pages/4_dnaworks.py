"""
pages/5_dnaworks.py
DNAWorks 올리고 설계 페이지 — dnaworks_code_light.html 대응.
"""
import datetime
import sys
from pathlib import Path
import streamlit as st

_app_dir = Path(__file__).resolve().parents[1]
for _p in (str(_app_dir), str(_app_dir.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from lib.session import init_session, get_rf_output_dir, get_ld_output_dir, LD_OUTPUT_BASE
from lib.html_loader import load_and_inject_css
from lib.page_header import render_page_header

init_session()
load_and_inject_css(_app_dir / "styles" / "theme.css")

from amr_dock_backend import (
    build_dnaworks_input, parse_dnaworks_logfile,
    get_tool_status, DNAWORKS_EXECUTABLE,
    save_timing, run_subprocess,
)


def _extract_binder_seq(pdb_path) -> str:
    """best.pdb에서 바인더 체인(chain B, 없으면 마지막 체인) 서열만 추출."""
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))
        chains = list(structure[0].get_chains())
        if not chains:
            return ""
        # chain B 우선, 없으면 마지막 체인
        target = next((c for c in chains if c.id == "B"), chains[-1])
        seq = ""
        for residue in target:
            if residue.get_id()[0] == " ":
                try:
                    aa = seq1(residue.get_resname())
                    if aa and aa != "X" and aa != "*":
                        seq += aa
                except Exception:
                    continue
        return seq
    except Exception as e:
        return f"Error: {e}"

import pandas as pd

render_page_header("내성 억제 단백질 합성을 위한 숙주 (DNA 합성 공장) 발현 환경 설정 및 DNA 서열 설계")

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

/* 슬라이더 라벨 / 값 */
[data-testid="stMainBlockContainer"] [data-testid="stSlider"] label {
  font-size: 1.275rem !important;
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

/* 체크박스 / 라디오 */
[data-testid="stMainBlockContainer"] .stCheckbox label,
[data-testid="stMainBlockContainer"] .stRadio label {
  font-size: 1.35rem !important;
}

/* expander 헤더 */
[data-testid="stMainBlockContainer"] [data-testid="stExpander"] summary {
  font-size: 1.3125rem !important;
}

/* number input */
[data-testid="stMainBlockContainer"] .stNumberInput [data-baseweb="input"] {
  font-size: 1.3125rem !important;
}
</style>
""", unsafe_allow_html=True)

st.markdown("<div style='margin-top:1.5rem;'></div>", unsafe_allow_html=True)

@st.cache_resource
def _cached_tool_status():
    return get_tool_status()

tool_status = _cached_tool_status()
dnaworks_ok = tool_status.get("dnaworks", False)
if not dnaworks_ok:
    st.error(f"❌ DNAWorks 실행 파일을 찾을 수 없습니다: `{DNAWORKS_EXECUTABLE}`")

rf_dir = get_rf_output_dir()
ld_dir = get_ld_output_dir()
job    = st.session_state.get("ifg_job_name", "design_01")
work_dir = ld_dir  # DNAWorks 결과는 LightDock job 폴더에 저장

# ─────────────────────────────────────────────────────────────────────────────
# 섹션 1: 입력 서열 준비
# ─────────────────────────────────────────────────────────────────────────────
st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">input</span>내성 억제 단백질 아미노산 합성 후보 입력</h3>', unsafe_allow_html=True)
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)

with st.container(border=True):
    seq_mode = st.radio(
        "억제 단백질 아미노산 합성 후보",
        ["best.pdb에서 추출", "직접 입력"],
        horizontal=True,
        key="dna_seq_mode",
    )

    if seq_mode == "best.pdb에서 추출":
        best_pdb = rf_dir / "best.pdb"

        # 페이지 진입 시 서열이 비어있으면 자동 로드
        if not st.session_state.get("ifg_dnaworks_seq") and best_pdb.exists():
            st.session_state.ifg_dnaworks_seq = _extract_binder_seq(best_pdb)

        if st.button("🔄 best.pdb에서 바인더 서열 다시 로드", use_container_width=False):
            if best_pdb.exists():
                st.session_state.ifg_dnaworks_seq = _extract_binder_seq(best_pdb)
                st.rerun()
            else:
                st.error("best.pdb를 찾을 수 없습니다.")

    protein_seq = st.text_area(
        "아미노산 서열 (Protein Sequence)",
        value=st.session_state.get("ifg_dnaworks_seq", ""),
        height=130,
        help="추출된 서열을 확인하거나 직접 편집할 수 있습니다.",
        key="dna_seq_input",
    )
st.session_state.ifg_dnaworks_seq = protein_seq
st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# 섹션 2: 숙주 맞춤형 발현 환경 설정
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)
st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">tune</span>맞춤형 숙주 발현 환경 설정</h3>', unsafe_allow_html=True)
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)

with st.container(border=True):
    col1, col2, col3 = st.columns(3)
    with col1:
        _CODON_OPTIONS = {
            "E. coli (대장균)":              "E. coli",
            "P. pastoris (피키아 효모)":      "P. pastoris",
            "S. cerevesiae (사카로미세스 효모)": "S. cerevesiae",
            "H. sapiens (인간)":             "H. sapiens",
            "M. musculus (생쥐)":            "M. musculus",
        }
        _codon_label = st.selectbox(
            "코돈 역번역 생물종 (Codon)",
            options=list(_CODON_OPTIONS.keys()),
            index=0,
            key="dna_codon",
        )
        codon_org = _CODON_OPTIONS[_codon_label]
    with col2:
        melting_temp = st.slider("목표 융해 온도 (Melting Temp, °C)",
                                 min_value=58, max_value=70, value=60, step=1, key="dna_tm")
    with col3:
        oligo_length = st.slider("올리고 길이 (Length, nt)",
                                 min_value=40, max_value=60, value=50, step=1, key="dna_ol")

    with st.expander("고급 설정 파라미터"):
        col_ex1, col_ex2 = st.columns(2)
        with col_ex1:
            hairpin_check  = st.checkbox("헤어핀 발생 검사 적용", value=True, key="dna_hairpin")
            sodium_conc    = st.number_input("나트륨 이온 농도 (M)", value=0.050, format="%.3f", key="dna_na")
            magnesium_conc = st.number_input("마그네슘 이온 농도 (M)", value=0.002, format="%.3f", key="dna_mg")
            n_solutions    = st.slider("설계 후보 수 (Solutions)", min_value=3, max_value=10, value=5, step=1, key="dna_sol")
            st.caption(
                        "💡 DNAWorks는 시뮬레이티드 어닐링 최적화를 이 횟수만큼 독립적으로 실행합니다. "
                        "많이 만들어서 거르는 게 아니라, 처음부터 이 수만큼만 생성하고 끝납니다."
                        )
        with col_ex2:
            misprime_check = st.checkbox("오결합/밀림(Misprime) 방어 적용", value=True, key="dna_misprime")
            repeat_limit   = st.number_input("반복 패턴 서열 제한 (nt)", min_value=4, max_value=20, value=8, step=1, key="dna_repeat")

st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# 실행 버튼
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
btn_r1, btn_r2, _ = st.columns([1, 1, 3])
with btn_r1:
    run_dna = st.button(
        "▶ 설계 실행",
        key="dna_run",
        type="primary",
        use_container_width=True,
        disabled=not dnaworks_ok or not protein_seq.strip(),
    )
with btn_r2:
    if st.button("🔄 초기화", key="dna_reset", use_container_width=True):
        st.session_state.ifg_dnaworks_seq = ""
        st.rerun()

# ─────────────────────────────────────────────────────────────────────────────
# 실행 버튼 바로 아래: 인라인 상태 / 소요 시간 표시
# ─────────────────────────────────────────────────────────────────────────────
# 1) 디스크 기준 상태/타이밍을 매 렌더마다 새로 계산한다.
#    세션 상태는 "running"/"error" 같은 일시 상태 전달용으로만 사용.
_timing_file = work_dir / "step_timing.json"
_timing_data: dict = {}
if _timing_file.exists():
    try:
        import json as _json
        _timing_data = _json.loads(_timing_file.read_text()).get("dnaworks", {})
    except Exception:
        pass

_dnaworks_done = (work_dir / "LOGFILE.txt").exists()


def _fmt_elapsed(dur_str) -> str:
    """'117.59s' → '0:01:57' 형식. 비어있으면 '--:--:--'."""
    if not dur_str or dur_str in ("—", None):
        return "--:--:--"
    try:
        secs = float(str(dur_str).rstrip("s"))
        h, rem = divmod(int(secs), 3600)
        m, sec = divmod(rem, 60)
        return f"{h}:{m:02d}:{sec:02d}"
    except Exception:
        return "--:--:--"


def _resolve_status() -> str:
    """현재 상태를 결정. running/error 는 세션 플래그 우선, 그 외엔 디스크 기준."""
    sess = st.session_state.get("ifg_dnaworks_status")
    if sess in ("running", "error"):
        return sess
    return "done" if _dnaworks_done else "waiting"


# 세션 상태를 디스크 기준으로 동기화 (stale 한 "waiting" 이 "done" 을 가리지 않도록).
st.session_state.ifg_dnaworks_status = _resolve_status()


st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
_inline_status_slot = st.empty()


def _render_inline_status(status_key: str | None = None, elapsed: str | None = None) -> None:
    sts = status_key if status_key is not None else _resolve_status()
    lbl = {"waiting": "실행 전", "running": "Running", "done": "Success", "error": "Error"}.get(sts, "실행 전")
    color = {"waiting": "#888888", "running": "#005ac1", "done": "#006c4c", "error": "#b02638"}.get(sts, "#888888")

    if elapsed is None:
        if sts == "done":
            # 세션 값 우선, 없으면 디스크에서 직접 로드
            td = st.session_state.get("ifg_dnaworks_timing") or _timing_data or {}
            elapsed = _fmt_elapsed(td.get("duration"))
        else:
            elapsed = "--:--:--"
    elif not elapsed:
        elapsed = "--:--:--"

    _inline_status_slot.markdown(
        f"""
        <div style="font-size:1.35rem;color:#333;line-height:1.8;">
          상태 : <b style="color:{color};">{lbl}</b><br>
          소요 시간 : <b>{elapsed}</b>
        </div>
        """,
        unsafe_allow_html=True,
    )


_render_inline_status()
st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)

if run_dna and protein_seq.strip():
    work_dir.mkdir(parents=True, exist_ok=True)
    inp_path = work_dir / "DNAWORKS.inp"
    log_path = work_dir / "dnaworks_run.log"

    inp_path.write_text(
        build_dnaworks_input({
            "title":           f"InhibiForge DNAWorks — {job}",
            "logfile":         "LOGFILE.txt",
            "melting_temp":    melting_temp,
            "oligo_length":    oligo_length,
            "codon_org":       codon_org,
            "sodium_conc":     sodium_conc,
            "magnesium_conc":  magnesium_conc,
            "repeat_limit":    int(repeat_limit),
            "n_solutions":     n_solutions,
            "misprime_check":  misprime_check,
            "hairpin_check":   hairpin_check,
            "protein_seq":     protein_seq.strip(),
        }),
        encoding="utf-8",
    )

    st.session_state.ifg_dnaworks_status = "running"
    st.session_state.ifg_dnaworks_timing = {}
    _render_inline_status("running", "--:--:--")

    start_t = datetime.datetime.now()
    rc = run_subprocess(
        [str(DNAWORKS_EXECUTABLE), "DNAWORKS.inp"],
        cwd=work_dir,
        log_path=log_path,
    )
    end_t = datetime.datetime.now()
    save_timing(work_dir, "dnaworks", start_t, end_t)

    # 타이밍 데이터 갱신
    if _timing_file.exists():
        try:
            import json as _json2
            st.session_state.ifg_dnaworks_timing = _json2.loads(_timing_file.read_text()).get("dnaworks", {})
        except Exception:
            pass

    if rc == 0:
        st.session_state.ifg_dnaworks_status = "done"
        _render_inline_status()
        st.success("✅ DNAWorks 실행 완료!")
        st.rerun()
    else:
        st.session_state.ifg_dnaworks_status = "error"
        _render_inline_status()
        st.error(f"❌ DNAWorks 실행 실패 (return code {rc})")
        if log_path.exists():
            with st.expander("실행 로그"):
                st.code(log_path.read_text(errors="replace")[-2000:], language="text")

# ─────────────────────────────────────────────────────────────────────────────
# 섹션 3: 결과 리포트
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)
st.markdown('<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">assessment</span>최종 DNA 합성 설계 리포트</h3>', unsafe_allow_html=True)
st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
st.info(
    "**결과 읽는 법** — "
    "각 솔루션은 동일한 단백질에 대한 독립적인 설계안입니다. "
    "하나의 솔루션 안에 있는 올리고 N개는 그 단백질을 PCR로 조립하기 위한 DNA 조각들이에요. "
    "레고 비유로 말하면, 솔루션 = 완성된 레고 작품 / 올리고 = 그걸 만드는 블록 조각들입니다. "
    "**Score가 낮을수록 최적**이며, 🏆 Best 솔루션의 올리고 세트를 우선 사용하세요."
)

dnaworks_log = work_dir / "LOGFILE.txt"
if not dnaworks_log.exists():
    st.info("💡 '설계 실행' 버튼을 눌러 DNAWorks를 시작하세요.")
else:
    try:
        trial_oligos, summary_rows = parse_dnaworks_logfile(dnaworks_log)
    except Exception as e:
        st.error(f"로그 파싱 오류: {e}")
        trial_oligos, summary_rows = {}, []

    if not summary_rows:
        st.warning(
            "⚠️ LOGFILE.txt는 존재하지만 점수 데이터(FINAL SUMMARY)를 찾을 수 없습니다. "
            "아래 '📄 LOGFILE 원본 다운로드'를 눌러 내용을 확인하세요. "
            "단백질 서열 오류 또는 DNAWorks 실행 중 문제가 원인일 수 있습니다."
        )
    else:
        # Score 오름차순 → TmRange 오름차순 (동점 시 더 균일한 Tm 범위가 Best)
        summary_df = (
            pd.DataFrame(summary_rows)
            .sort_values(["Score", "TmRange"])
            .reset_index(drop=True)
        )
        summary_df.insert(0, "순위", ["🏆" if i == 0 else str(i + 1) for i in range(len(summary_df))])

        # Score 0.000이 2개 이상일 때 TmRange 2차 정렬 안내
        zero_score_count = (summary_df["Score"] == 0.0).sum()
        if zero_score_count > 1:
            st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
            st.info(
                f"ℹ️ Score 0.000인 솔루션이 {zero_score_count}개 감지되었습니다. "
                "TmRange(융해 온도 범위) 오름차순으로 2차 정렬을 적용합니다. "
                "TmRange가 작을수록 올리고 간 Tm이 균일하여 PCR 조립 효율이 높습니다."
            )

        st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)
        st.markdown('<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">format_list_bulleted</span> 내성 억제 단백질 합성용 DNA 서열 목록</h4>', unsafe_allow_html=True)
        st.markdown("<div style='margin-top:2rem;'></div>", unsafe_allow_html=True)
        for rank_i, row in summary_df.iterrows():
            trial_num = int(row["Trial"])
            score_val = row["Score"]
            tmrange_val = row["TmRange"]
            rank_label = "🏆 Best" if rank_i == 0 else f"#{rank_i + 1}"
            oligos = trial_oligos.get(trial_num, [])

            # Score 0.000일 때 TmRange 정보도 표시
            if score_val == 0.0 and zero_score_count > 1:
                label_detail = (
                    f"설계 적합도 점수: {score_val:.3f} | "
                    f"융해 온도 허용 범위: {tmrange_val:.1f}°C"
                )
            else:
                label_detail = f"설계 적합도 점수: {score_val:.3f}"

            if rank_i > 0:
                st.markdown("<div style='margin-top:3rem;'></div>", unsafe_allow_html=True)
            st.markdown(
                f"<div style='font-weight:600;font-size:2.25rem;margin-bottom:1.5rem;'>"
                f"{rank_label}&nbsp;&nbsp;Trial {trial_num}&nbsp;&nbsp;—&nbsp;&nbsp;{label_detail}"
                f"&nbsp;&nbsp;|&nbsp;&nbsp;DNA 조각 {len(oligos)}개</div>",
                unsafe_allow_html=True,
            )
            with st.container(border=True):
                if oligos:
                    try:
                        from Bio.SeqUtils import MeltingTemp as _mt
                        _tm_available = True
                    except Exception:
                        _tm_available = False

                    def _calc_gc(seq: str) -> float:
                        if not seq:
                            return 0.0
                        s = seq.upper()
                        return (s.count("G") + s.count("C")) / len(s) * 100.0

                    def _calc_tm(seq: str) -> float:
                        if not seq:
                            return 0.0
                        if _tm_available:
                            try:
                                return float(_mt.Tm_NN(
                                    seq.upper(),
                                    Na=sodium_conc * 1000,
                                    Mg=magnesium_conc * 1000,
                                ))
                            except Exception:
                                pass
                        # Fallback: Wallace rule (부정확하지만 동작은 보장)
                        s = seq.upper()
                        at = s.count("A") + s.count("T")
                        gc = s.count("G") + s.count("C")
                        return 2 * at + 4 * gc

                    rows = []
                    for o in oligos:
                        seq = o.get("Sequence (DNA)") or o.get("Sequence") or ""
                        oid = int(o.get("ID", 0))
                        length = int(o.get("Length", len(seq)))
                        rows.append({
                            "DNA 조각 ID": oid,
                            "DNA 서열": seq,
                            "서열길이": length,
                            "융해 온도 (Tm, °C)": round(_calc_tm(seq), 1),
                            "GC 함량 (%)": round(_calc_gc(seq), 1),
                            "가닥 구분": "Top" if oid % 2 == 1 else "Bottom",
                        })
                    oligo_df = pd.DataFrame(
                        rows,
                        columns=["DNA 조각 ID", "DNA 서열", "서열길이", "융해 온도 (Tm, °C)", "GC 함량 (%)", "가닥 구분"],
                    )
                    _oligo_styled = (
                        oligo_df.style
                        .set_properties(**{
                            "background-color": "white",
                            "text-align": "left",
                        })
                        .set_table_styles([
                            {"selector": "th", "props": [("text-align", "center")]},
                            {"selector": "thead th", "props": [("text-align", "center")]},
                        ])
                    )
                    st.dataframe(
                        _oligo_styled,
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            "융해 온도 (Tm, °C)": st.column_config.NumberColumn(
                                "융해 온도 (Tm, °C)", format="%.1f",
                            ),
                            "GC 함량 (%)": st.column_config.NumberColumn(
                                "GC 함량 (%)", format="%.1f",
                            ),
                        },
                    )
                    st.download_button(
                        f"📥 Trial {trial_num} CSV 다운로드",
                        data=oligo_df.to_csv(index=False).encode("utf-8-sig"),
                        file_name=f"dna_oligos_trial{trial_num}_score{score_val:.3f}.csv",
                        mime="text/csv",
                        use_container_width=True,
                        key=f"dna_csv_{trial_num}",
                    )
                else:
                    st.info("이 Trial의 올리고 데이터를 찾을 수 없습니다.")

    # 원본 로그 다운로드
    st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)
    raw_log = dnaworks_log.read_text(encoding="utf-8", errors="replace")
    st.download_button(
        "LOGFILE 원본 다운로드",
        data=raw_log.encode("utf-8"),
        file_name="LOGFILE.txt",
        mime="text/plain",
        use_container_width=False,
        key="dna_log_dl",
    )
