"""
pages/5_file_storage.py
파일 보관함 페이지 — 프로젝트 브라우저 (Stage 2).

통합 Projects/ 루트 아래의 각 프로젝트 폴더를 버튼으로 나열한다:
  - 버튼: 프로젝트 이름 (클릭 시 ifg_job_name 업데이트 후 Pipeline 페이지로 이동)
  - 버튼 옆: 생성 시간, 파이프라인 종료 시간 (KST)
  - 종료 시간은 step_timing.json 이 있으면 거기의 가장 마지막 end, 없으면 '—'

하단에는 현재 선택된 프로젝트의 로그 파일 섹션을 유지해 다운로드 가능.
"""
import datetime
import json
import sys
from datetime import timezone, timedelta
from pathlib import Path

import streamlit as st

_app_dir = Path(__file__).resolve().parents[1]
for _p in (str(_app_dir), str(_app_dir.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from lib.session import (
    init_session, PROJECTS_BASE, get_project_dir,
    get_rf_output_dir, get_ld_output_dir,
)
from lib.html_loader import load_and_inject_css
from lib.file_browser import list_job_files, format_file_size
from lib.page_header import render_page_header

init_session()
load_and_inject_css(_app_dir / "styles" / "theme.css")

KST = timezone(timedelta(hours=9))

render_page_header(
    "파일 보관함",
    "프로젝트를 선택해 파이프라인 설정부터 DNAWorks 까지 모든 결과를 불러옵니다",
)

st.markdown("<div style='margin-top:1.5rem;'></div>", unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────────────────────
# 시간 유틸 (KST)
# ─────────────────────────────────────────────────────────────────────────────
def _fmt_kst(dt: datetime.datetime | None) -> str:
    if dt is None:
        return "—"
    return dt.astimezone(KST).strftime("%Y-%m-%d %H:%M")


def _folder_created_kst(p: Path) -> datetime.datetime | None:
    try:
        st_ = p.stat()
    except OSError:
        return None
    # Linux 에서는 st_ctime 이 inode 변경 시각. 최신 성공 실행 기준은 st_mtime 이 더 의미있음.
    # 생성 시각 근사치: min(st_ctime, st_mtime).
    ts = min(st_.st_ctime, st_.st_mtime)
    return datetime.datetime.fromtimestamp(ts, tz=KST)


def _parse_kst_or_naive(s: str) -> datetime.datetime | None:
    """step_timing.json 이 v14 는 '%Y-%m-%d %H:%M:%S%z', 구버전은 tz 없는 문자열."""
    s = s.strip()
    if not s:
        return None
    # 1) +0900 / +09:00 포함하는 aware 포맷 먼저 시도
    for fmt in ("%Y-%m-%d %H:%M:%S%z", "%Y-%m-%d %H:%M:%S.%f%z"):
        try:
            return datetime.datetime.strptime(s, fmt)
        except ValueError:
            pass
    # 2) naive → KST 로 간주
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%S.%f"):
        try:
            return datetime.datetime.strptime(s, fmt).replace(tzinfo=KST)
        except ValueError:
            pass
    return None


def _project_finished_kst(project_dir: Path) -> datetime.datetime | None:
    """step_timing.json 의 가장 마지막 end 값을 KST datetime 으로. 없으면 None."""
    tf = project_dir / "step_timing.json"
    if not tf.exists():
        return None
    try:
        data = json.loads(tf.read_text(encoding="utf-8"))
    except Exception:
        return None
    latest: datetime.datetime | None = None
    for step_info in data.values():
        if not isinstance(step_info, dict):
            continue
        end_s = step_info.get("end") or ""
        dt = _parse_kst_or_naive(str(end_s))
        if dt is None:
            continue
        if latest is None or dt > latest:
            latest = dt
    return latest


# ─────────────────────────────────────────────────────────────────────────────
# 프로젝트 목록 수집
# ─────────────────────────────────────────────────────────────────────────────
def _discover_projects() -> list[dict]:
    if not PROJECTS_BASE.exists():
        return []
    projects: list[dict] = []
    for p in PROJECTS_BASE.iterdir():
        if not p.is_dir():
            continue
        projects.append({
            "name":     p.name,
            "path":     p,
            "created":  _folder_created_kst(p),
            "finished": _project_finished_kst(p),
        })
    # 최신 생성순 내림차순
    projects.sort(
        key=lambda d: d["created"] or datetime.datetime.min.replace(tzinfo=KST),
        reverse=True,
    )
    return projects


# ─────────────────────────────────────────────────────────────────────────────
# 로그 수집 유틸 (하단 로그 섹션용)
# ─────────────────────────────────────────────────────────────────────────────
_LOG_NAME_WHITELIST = {"pipeline.log", "LOGFILE.txt"}
_LOG_EXTENSIONS     = {".log", ".txt", ".out"}


def _collect_logs(files: dict) -> list[Path]:
    logs: list[Path] = list(files.get("logs", []))
    for f in files.get("intermediate", []):
        if f.name in _LOG_NAME_WHITELIST or f.suffix.lower() in _LOG_EXTENSIONS:
            logs.append(f)
    seen: set[Path] = set()
    unique: list[Path] = []
    for f in logs:
        if f not in seen:
            seen.add(f)
            unique.append(f)
    try:
        unique.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    except Exception:
        pass
    return unique


def _render_log_table(title: str, files: list[Path], icon: str = "description") -> None:
    if not files:
        return
    st.markdown(
        f'<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">{icon}</span> {title}</h4>',
        unsafe_allow_html=True,
    )
    for i, f in enumerate(files):
        try:
            size_str = format_file_size(f.stat().st_size)
            mtime_dt = datetime.datetime.fromtimestamp(f.stat().st_mtime, tz=KST)
            mtime_str = mtime_dt.strftime("%Y.%m.%d %H:%M")
        except Exception:
            size_str = "—"
            mtime_str = "—"

        cols = st.columns([5, 2, 2, 1])
        with cols[0]:
            st.markdown(
                f'<span style="font-size:0.85rem;font-weight:500;">{f.name}</span>',
                unsafe_allow_html=True,
            )
        with cols[1]:
            st.markdown(
                f'<span style="font-size:0.78rem;color:#888;">{mtime_str}</span>',
                unsafe_allow_html=True,
            )
        with cols[2]:
            st.markdown(
                f'<span style="font-size:0.78rem;color:#888;">{size_str}</span>',
                unsafe_allow_html=True,
            )
        with cols[3]:
            try:
                st.download_button(
                    "⬇",
                    data=f.read_bytes(),
                    file_name=f.name,
                    key=f"dl_storage_{title}_{i}_{f.name}",
                    use_container_width=True,
                )
            except Exception:
                pass


# ─────────────────────────────────────────────────────────────────────────────
# 상단: 프로젝트 브라우저
# ─────────────────────────────────────────────────────────────────────────────
projects      = _discover_projects()
current_job   = st.session_state.get("ifg_job_name", "")

with st.container(border=True):
    st.markdown(
        f'<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">folder_open</span> '
        f'프로젝트 선택 — 총 {len(projects)} 건</h3>',
        unsafe_allow_html=True,
    )
    if current_job:
        st.caption(f"현재 선택된 프로젝트: **{current_job}**")

    if not projects:
        st.info(
            f"아직 프로젝트가 없습니다. `{PROJECTS_BASE}` 아래에 프로젝트 폴더가 "
            "생성되면 이 목록에 표시됩니다."
        )
    else:
        # 헤더
        header = st.columns([5, 3, 3])
        header[0].markdown("**프로젝트**")
        header[1].markdown("**생성 시간 (KST)**")
        header[2].markdown("**파이프라인 종료 (KST)**")
        st.divider()

        for proj in projects:
            cols = st.columns([5, 3, 3])
            is_current = (proj["name"] == current_job)
            btn_type = "primary" if is_current else "secondary"
            with cols[0]:
                if st.button(
                    proj["name"],
                    key=f"proj_btn_{proj['name']}",
                    use_container_width=True,
                    type=btn_type,
                ):
                    st.session_state["ifg_job_name"] = proj["name"]
                    # Pipeline 페이지로 이동 — 이 페이지부터 DNAWorks 까지 해당 프로젝트 데이터가 표시됨
                    try:
                        st.switch_page("pages/1_pipeline.py")
                    except Exception:
                        # switch_page 미지원 환경: rerun 으로 현재 페이지만 갱신
                        st.rerun()
            with cols[1]:
                st.markdown(
                    f'<div style="font-size:0.95rem;color:#555;padding-top:0.5rem;">'
                    f'{_fmt_kst(proj["created"])}</div>',
                    unsafe_allow_html=True,
                )
            with cols[2]:
                st.markdown(
                    f'<div style="font-size:0.95rem;color:#555;padding-top:0.5rem;">'
                    f'{_fmt_kst(proj["finished"])}</div>',
                    unsafe_allow_html=True,
                )


# ─────────────────────────────────────────────────────────────────────────────
# 하단: 현재 선택된 프로젝트의 실행 로그 섹션 (RF / LD)
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("<div style='margin-top:3rem;'></div>", unsafe_allow_html=True)

rf_dir = get_rf_output_dir()
ld_dir = get_ld_output_dir()
job    = current_job or "—"

with st.container(border=True):
    st.markdown(
        f'<h3><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">description</span> '
        f'현재 선택된 프로젝트 로그 — <code>{job}</code></h3>',
        unsafe_allow_html=True,
    )
    col_main, col_side = st.columns([8, 4], gap="large")

    with col_main:
        st.markdown(
            '<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">biotech</span> RFdiffusion / MPNN / AF2 실행 로그</h4>',
            unsafe_allow_html=True,
        )
        rf_files = list_job_files(rf_dir)
        rf_logs  = _collect_logs(rf_files)
        _render_log_table("실행 로그", rf_logs, "description")
        if not rf_logs:
            st.info("파이프라인을 실행하면 여기에 RFdiffusion 관련 로그가 표시됩니다.")

    with col_side:
        st.markdown(
            '<h4><span class="material-symbols-outlined" style="color:#005ac1;font-size:inherit;vertical-align:middle">anchor</span> LightDock 실행 로그</h4>',
            unsafe_allow_html=True,
        )
        ld_files = list_job_files(ld_dir)
        ld_logs  = _collect_logs(ld_files)
        _render_log_table("실행 로그", ld_logs, "description")
        if not ld_logs:
            st.info("LightDock 완료 후 로그가 표시됩니다.")

st.markdown("<div style='margin-top:4rem;'></div>", unsafe_allow_html=True)
