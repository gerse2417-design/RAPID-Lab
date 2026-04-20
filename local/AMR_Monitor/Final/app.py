import streamlit as st
import streamlit.components.v1 as components
import os
import time
import json
import shutil
import base64
import io
from datetime import datetime
import pandas as pd
from pipeline import FinalWorkflowPipeline
from database import init_db
import py3Dmol
from stmol import showmol
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

# --- Set Global Font for Korean ---
font_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "NanumGothic.ttf")
if os.path.exists(font_path):
    fm.fontManager.addfont(font_path)
    font_name = fm.FontProperties(fname=font_path).get_name()
    plt.rc("font", family=font_name)
    plt.rcParams["axes.unicode_minus"] = False

# --- Helper for Execution Timing ---
def render_execution_log(log_data):
    """실행 로그를 커스텀 스타일의 HTML 테이블로 출력합니다."""
    import textwrap
    
    total_seconds = 0
    if log_data:
        try:
            for entry in log_data:
                h, m, s = map(int, entry['elapsed'].split(':'))
                total_seconds += h * 3600 + m * 60 + s
        except: pass
    
    hrs, rem = divmod(total_seconds, 3600)
    mins, secs = divmod(rem, 60)
    time_str = f"{hrs:02d}:{mins:02d}:{secs:02d}"

    html_str = f"""
<div style="background-color: #ffffff; border-radius: 0.75rem; border: 1px solid rgba(0,0,0,0.1); box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05); padding: 1.5rem; margin-bottom: 2rem;">
    <div style="display: flex; align-items: center; gap: 0.5rem; margin-bottom: 1.25rem;">
        <span class="material-symbols-outlined" style="color: #2563eb; font-size: 1.25rem; font-variation-settings: 'FILL' 0;">timer</span>
        <h2 style="margin: 0; font-size: 1.125rem; font-weight: 700; color: #111827; font-family: 'Space Grotesk', sans-serif;">프로그램 실행 정보</h2>
    </div>
    
    <div style="background-color: #f9fafb; border: 1px solid #e5e7eb; border-radius: 0.375rem; overflow: hidden;">
        <table style="width: 100%; text-align: left; border-collapse: collapse;">
            <thead style="background-color: #f3f4f6; border-bottom: 1px solid #e5e7eb;">
                <tr>
                    <th style="padding: 0.75rem 1rem; font-size: 0.6875rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: #4b5563;">작업명(TASK)</th>
                    <th style="padding: 0.75rem 1rem; font-size: 0.6875rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: #4b5563;">시작 시간(START)</th>
                    <th style="padding: 0.75rem 1rem; font-size: 0.6875rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: #4b5563;">종료 시간(END)</th>
                    <th style="padding: 0.75rem 1rem; font-size: 0.6875rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: #4b5563;">경과 시간(ELAPSED)</th>
                    <th style="padding: 0.75rem 1rem; font-size: 0.6875rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: #4b5563;">상태(STATUS)</th>
                </tr>
            </thead>
            <tbody style="font-size: 0.875rem; color: #4b5563;">
"""
    
    if not log_data:
        html_str += """
                <tr>
                    <td colspan="5" style="padding: 2rem 1rem; text-align: center; color: #6b7280; font-weight: 500;">서열 데이터가 입력되지 않았습니다.</td>
                </tr>
"""
    else:
        for entry in log_data:
            status = entry.get('status', '')
            status_color = "#059669" if "Success" in status else ("#e11d48" if "Fail" in status or "Error" in status else "#4b5563")
            
            html_str += f"""
                <tr style="border-bottom: 1px solid #e5e7eb; background-color: #ffffff;">
                    <td style="padding: 0.75rem 1rem; font-weight: 500; color: #111827;">{entry.get('task', '')}</td>
                    <td style="padding: 0.75rem 1rem;">{entry.get('start', '')}</td>
                    <td style="padding: 0.75rem 1rem;">{entry.get('end', '')}</td>
                    <td style="padding: 0.75rem 1rem; font-family: monospace;">{entry.get('elapsed', '')}</td>
                    <td style="padding: 0.75rem 1rem; font-weight: 600; color: {status_color};">{status}</td>
                </tr>
"""
            
    html_str += f"""
            </tbody>
        </table>
    </div>
    
    <div style="margin-top: 1rem; text-align: right; font-size: 0.875rem; color: #4b5563; font-family: monospace;">
        ※ 전체 경과 시간: <span style="color: #2563eb; font-weight: 700;">{time_str}</span>
    </div>
</div>
"""
    
    # 스트림릿 마크다운 파서가 들여쓰기를 코드 블록으로 잘못 해석하지 않도록 줄바꿈을 완전히 제거합니다.
    clean_html = ''.join(html_str.splitlines())
    st.markdown(clean_html, unsafe_allow_html=True)

    # Prodigal이 실행된 경우에만 설명 카드 표시
    has_prodigal = any("Prodigal" in str(entry.get("task", "")) for entry in (log_data or []))
    if has_prodigal:
        st.markdown("""
<div style="background:#f0f7ff;border:1px solid #bfdbfe;border-radius:0.5rem;padding:1rem 1.25rem;margin-top:-1rem;margin-bottom:1.5rem;">
    <div style="display:flex;align-items:flex-start;gap:0.75rem;">
        <span class="material-symbols-outlined" style="color:#2563eb;font-size:1.2rem;margin-top:0.1rem;">biotech</span>
        <div>
            <span style="font-size:0.9rem;font-weight:700;color:#1e3a5f;">Prodigal</span>
            <span style="font-size:0.85rem;color:#374151;margin-left:0.5rem;">세균의 유전체 서열을 주석화하여 유전적 요소를 식별하고 단백질 코딩 유전자를 예측하여 아미노산 서열을 제공하는 소프트웨어 도구</span>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

    # BLASTp가 실행된 경우에만 설명 카드 표시
    has_blastp = any("BLASTp" in str(entry.get("task", "")) for entry in (log_data or []))
    if has_blastp:
        st.markdown("""
<div style="background:#f0f7ff;border:1px solid #bfdbfe;border-radius:0.5rem;padding:1rem 1.25rem;margin-top:-1rem;margin-bottom:1.5rem;">
    <div style="display:flex;align-items:flex-start;gap:0.75rem;">
        <span class="material-symbols-outlined" style="color:#2563eb;font-size:1.2rem;margin-top:0.1rem;">search</span>
        <div style="display:flex;flex-wrap:wrap;align-items:center;gap:0.5rem;">
            <span style="font-size:0.9rem;font-weight:700;color:#1e3a5f;">BLASTp</span>
            <span style="font-size:0.85rem;color:#374151;">Usearch 소프트웨어를 통해 로컬 정렬 검색 도구인 BLAST를 수행 시, 매칭된 서열의 수(hit 수)와 기준 유사도 이상의 단백질을 찾아내는 서열의 유사성 검색 도구</span>
            <a href="https://bio-kcs.tistory.com/entry/BLAST-BLAST-%EC%95%8C%EA%B3%A0%EB%A6%AC%EC%A6%98%EC%97%90-%EB%8C%80%ED%95%B4-%EC%95%8C%EC%95%84%EB%B3%B4%EC%9E%90#4.%20BLAST-1" target="_blank" style="font-size:0.75rem;color:#2563eb;text-decoration:none;display:inline-flex;align-items:center;gap:0.2rem;background:white;padding:0.15rem 0.4rem;border:1px solid #bfdbfe;border-radius:0.25rem;margin-left:0.2rem;">
                <span class="material-symbols-outlined" style="font-size:0.85rem;">open_in_new</span>BLAST 알고리즘 이해하기
            </a>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

    # BOARDS DB 연동이 실행된 경우에만 설명 카드 표시
    has_boards = any("BOARDS DB 연동" in str(entry.get("task", "")) for entry in (log_data or []))
    if has_boards:
        st.markdown("""
<div style="background:#f0f7ff;border:1px solid #bfdbfe;border-radius:0.5rem;padding:1rem 1.25rem;margin-top:-1rem;margin-bottom:1.5rem;">
    <div style="display:flex;align-items:flex-start;gap:0.75rem;">
        <span class="material-symbols-outlined" style="color:#2563eb;font-size:1.2rem;margin-top:0.1rem;">database</span>
        <div>
            <span style="font-size:0.9rem;font-weight:700;color:#1e3a5f;">BOARDS DB</span>
            <span style="font-size:0.85rem;color:#374151;margin-left:0.5rem;">UNIST의 Systems Biology and Machine Learning Lab에서 구축한 항생제 내성 유전자 정보 및 항생제 내성 유발 단백질의 구조 예측 정보 데이터베이스</span>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

def get_plddt_color(val):
    if val >= 90: return "#0053D6"
    elif val >= 70: return "#65CBF3"
    elif val >= 50: return "#FFE066"
    else: return "#F08041"

# --- Helper for Structural Visualization ---
def show_structural_dashboard(data):
    """
    Displays the structural analytics dashboard matching BOARDS publication standards.
    Handles both BOARDS data (mapped from pipeline) and AF2 data.
    """
    if isinstance(data, dict):
        raw_id = data.get('target_id') or data.get('id')
        if raw_id:
            if "_BOARDS" in str(data.get('target_name', '')):
                radar_id = f"BOARDS DB Index #{raw_id}"
            else:
                radar_id = f"Model ID: {raw_id}"
        else:
            radar_id = "AlphaFold2"
        models = data.get('models', [])
        pae_image = data.get('pae_image')
    else:
        # Fallback for direct model list access
        radar_id = "AlphaFold2"
        models = data
        pae_image = None
    
    # Determine if this is a BOARDS result to handle conditional header
    is_boards = "BOARDS" in radar_id
    
    if is_boards:
        st.markdown(f"""
<div style="margin-bottom:1rem;">
    <p style="font-size:0.85rem;font-weight:600;color:#2563eb;letter-spacing:0.02em;margin:0 0 0.25rem 0;display:flex;align-items:center;gap:0.4rem;">
        <span style="font-size:1.1rem; transform: translateY(1px);">🔍</span>
        항생제 내성 유발 단백질의 예측 구조 및 신뢰도 지표
    </p>
    <h2 style="font-size:1.5rem;font-weight:700;color:#1e293b;margin:0;display:flex;align-items:center;gap:0.5rem;">
        <span style="transform: translateY(2px);">📊</span> {radar_id}의 단백질 3차원 구조 및 신뢰도 결과
    </h2>
</div>""", unsafe_allow_html=True)
    else:
        st.markdown(f"""
<div style="margin-bottom:1rem;">
    <h2 style="font-size:1.5rem;font-weight:700;color:#1e293b;margin:0;display:flex;align-items:center;gap:0.5rem;">
        <span style="transform: translateY(2px);">📊</span> {radar_id}의 단백질 3차원 구조 및 신뢰도 결과
    </h2>
</div>""", unsafe_allow_html=True)

    # [Section 1] PAE Plot

    # Check if we have ANY form of PAE results
    has_pae_data = False
    
    # Check for raw JSON data in current models
    active_models = models if models else st.session_state.get("af2_models", [])
    if isinstance(active_models, dict): active_models = active_models.get("models", [])
    
    potential_pae_data = any(np.array(m.get("pae_data", [])).size > 0 for m in active_models) if active_models else False
    has_pae_image = pae_image and os.path.exists(pae_image)
    
    if potential_pae_data or has_pae_image:
        html_content = """
<div class="sci-table-card" style="margin-bottom: 1rem;">
    <div class="sci-table-header" style="justify-content: space-between;">
        <div style="display:flex;align-items:center;gap:0.5rem;">
            <span style="font-size:1.2rem; display: flex; align-items: center;">📈</span>
            <h4 style="margin:0;font-size:0.95rem;font-weight:700;color:#1e293b; position: relative; top: 5px; line-height: 1;">전체 단백질 구조의 신뢰도: 예측 정렬 오차 (Predicted Aligned Error, PAE) 도표</h4>
        </div>
        <a href="https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/pae-a-measure-of-global-confidence-in-alphafold-predictions/" target="_blank" style="font-size:0.75rem;color:#2563eb;text-decoration:none;display:inline-flex;align-items:center;gap:0.25rem;background:white;padding:0.25rem 0.5rem;border:1px solid #bfdbfe;border-radius:0.25rem;transition:all 0.2s;box-shadow:0 1px 2px rgba(0,0,0,0.05);">
            <span class="material-symbols-outlined" style="font-size:0.875rem;">open_in_new</span>PAE 이해하기
        </a>
    </div>
    <div class="sci-table-body" style="text-align:center;">
"""
        
        # 1. Professional Grid View (from JSON)
        if potential_pae_data:
            n_models = min(5, len(active_models))
            fig_grid, axes = plt.subplots(1, n_models, figsize=(18, 5))
            if n_models == 1: axes = [axes]
            
            for i in range(n_models):
                ax = axes[i]
                m = active_models[i]
                pae = np.array(m.get("pae_data", []))
                
                if pae.size > 0 and len(pae.shape) == 2:
                    im = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=30)
                    ax.set_title(f"Model Rank {m.get('rank', i)} (pLDDT: {m.get('plddt', 0):.1f})", fontsize=10)
                    ax.set_xlabel("Scored residue", fontsize=9)
                    if i == 0:
                        ax.set_ylabel("Aligned residue", fontsize=9)
                    else:
                        ax.set_yticklabels([])
                else:
                    ax.text(0.5, 0.5, "PAE 데이터 없음", ha='center', va='center')

            plt.tight_layout(rect=[0, 0.1, 1, 0.95])
            cbar = fig_grid.colorbar(im, ax=axes.ravel().tolist(), orientation='horizontal', pad=0.15, fraction=0.04, aspect=40)
            cbar.set_label("Expected position error (Ångströms, Å)", fontsize=10, fontfamily='sans-serif')
            
            # Save fig to base64
            buf = io.BytesIO()
            fig_grid.savefig(buf, format="png", bbox_inches='tight', dpi=100)
            plt.close(fig_grid)
            buf.seek(0)
            img_b64 = base64.b64encode(buf.read()).decode('utf-8')
            
            html_content += f'<img src="data:image/png;base64,{img_b64}" style="max-width:100%; border-radius:0.25rem;" />'
            html_content += "</div></div>"
            st.markdown(html_content, unsafe_allow_html=True)

        
        # 2. Official Grid View (ONLY IF JSON is missing)
        elif has_pae_image:
            with open(pae_image, "rb") as f:
                img_b64 = base64.b64encode(f.read()).decode('utf-8')
            html_content += f'<img src="data:image/png;base64,{img_b64}" style="max-width:100%; border-radius:0.25rem;" />'
            html_content += "</div></div>"
            st.markdown(html_content, unsafe_allow_html=True)
            
        # --- PAE Explanation Guide ---
        st.markdown("""
<div style="background-color: #f0fdf4; border: 1px solid #bbf7d0; border-radius: 0.5rem; padding: 1.25rem; margin-top: 1rem; margin-bottom: 2rem;">
    <div style="display: flex; align-items: flex-start; gap: 0.75rem;">
        <span class="material-symbols-outlined" style="color: #166534; font-size: 1.4rem; font-variation-settings: 'FILL' 1;">info</span>
        <div style="flex: 1;">
            <p style="margin: 0 0 0.6rem 0; font-size: 0.9rem; font-weight: 700; color: #166534;">PAE에 대한 설명 및 그래프 해석 방법</p>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1.5rem;">
                <div>
                    <p style="margin: 0 0 0.5rem 0; font-size: 0.825rem; font-weight: 700; color: #166534;">PAE:</p>
                    <ul style="margin: 0; padding-left: 1.2rem; font-size: 0.8rem; color: #14532d; line-height: 1.6; list-style-type: disc;">
                        <li>잔기 쌍별 신뢰도 지표</li>
                        <li>예측된 구조에서 두 잔기 사이의 상대적인 위치와 방향에 대해 얼마나 확신하는지를 나타내는 척도</li>
                        <li>예측된 구조와 실제 구조가 잔기 Y를 기준으로 정렬되었을 때 잔기 X에서의 예측 위치 오차</li>
                    </ul>
                </div>
                <div>
                    <p style="margin: 0 0 0.5rem 0; font-size: 0.825rem; font-weight: 700; color: #166534;">해석:</p>
                    <ul style="margin: 0; padding-left: 1.2rem; font-size: 0.8rem; color: #14532d; line-height: 1.6; list-style-type: disc;">
                        <li><b>녹색 음영:</b> 두 잔기 사이의 예상 거리 오차 (단위: 옹스트롬, Å)</li>
                        <li>PAE 값이 낮을수록 잘 정의된 상대적인 위치와 방향을 예측</li>
                        <li>PAE 값이 높을수록 상대적인 위치 또는 방향이 불확실</li>
                        <li><b>진한 녹색:</b> 낮은 오차, 높은 신뢰도</li>
                        <li><b>연한 녹색:</b> 높은 오차, 낮은 신뢰도</li>
                    </ul>
                </div>
            </div>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)
            
    else:
        # No PAE data available
        st.warning("⚠️ BOARDS DB에 해당 단백질의 Predicted Aligned Error(PAE) 데이터가 포함되어 있지 않습니다.")

    
    # [Section 2] pLDDT Quality Metrics (Profiles & Averages)
    st.markdown("""
<div style="margin:2.5rem 0 1.25rem 0;display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:0.6rem;">
    <div style="display:flex;align-items:center;gap:0.6rem;">
        <span style="font-size:1.4rem; display: flex; align-items: center;">🧪</span>
        <h3 style="margin:0;font-size:1.15rem;font-weight:700;color:#1e293b; position: relative; top: 2px;">단백질 구조의 영역별 신뢰도: pLDDT (Predicted Local Distance Difference Test) 점수</h3>
    </div>
    <a href="https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/" target="_blank" style="font-size:0.75rem;color:#2563eb;text-decoration:none;display:inline-flex;align-items:center;gap:0.25rem;background:white;padding:0.25rem 0.5rem;border:1px solid #bfdbfe;border-radius:0.25rem;transition:all 0.2s;box-shadow:0 1px 2px rgba(0,0,0,0.05);">
        <span class="material-symbols-outlined" style="font-size:0.875rem;">open_in_new</span>pLDDT 점수 이해하기
    </a>
</div>
""", unsafe_allow_html=True)
    
    col1, col2 = st.columns([3, 2])
    
    # Extract data for plotting
    all_plddts = [m['plddt_list'] for m in models if 'plddt_list' in m]
    avg_plddts = [m.get('plddt', 0) for m in models]
    titles = [f"Rank {m['rank']}" for m in models]
    
    if all_plddts:
        with col1:
            # pLDDT Profile Plot
            fig1, ax1 = plt.subplots(figsize=(8, 4.5))
            colors = plt.cm.tab10(np.linspace(0, 1, 10))
            for i, (plddt, title) in enumerate(zip(all_plddts, titles)):
                ax1.plot(plddt, label=title, alpha=0.8, linewidth=1.2, color=colors[i%10])
            
            ax1.axhline(y=70, color='red', linestyle='--', alpha=0.5, label='Threshold (70)')
            ax1.set_xlabel("Positions", fontsize=9)
            ax1.set_ylabel("Predicted LDDT(pLDDT)", fontsize=9)
            ax1.set_ylim(0, 105)
            ax1.legend(loc='lower right', fontsize='x-small', ncol=2)
            ax1.grid(True, linestyle=':', alpha=0.6)
            
            buf1 = io.BytesIO()
            fig1.savefig(buf1, format="png", bbox_inches='tight', dpi=100)
            plt.close(fig1)
            buf1.seek(0)
            img1_b64 = base64.b64encode(buf1.read()).decode('utf-8')
            
            st.markdown(f"""
<div class="sci-table-card" style="height:100%; margin-bottom:0;">
    <div class="sci-table-header">
        <span style="font-size:1.2rem; display: flex; align-items: center;">📈</span>
        <h4 style="margin:0;font-size:0.95rem;font-weight:700;color:#1e293b; position: relative; top: 5px; line-height: 1;">아미노산 잔기별 신뢰도 도표 (Per-residue Confidence Plot)</h4>
    </div>
    <div class="sci-table-body" style="text-align:center;">
        <img src="data:image/png;base64,{img1_b64}" style="max-width:100%; border-radius:0.25rem;" />
    </div>
</div>
""", unsafe_allow_html=True)
            
        with col2:
            # Average pLDDT Histogram
            fig2, ax2 = plt.subplots(figsize=(6, 4.5))
            # Use Viridis/Spectral color map for a professional look
            bar_colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(avg_plddts)))
            bars = ax2.bar(titles, avg_plddts, color=bar_colors, edgecolor='black', alpha=0.7)
            
            ax2.set_ylim(0, 105)
            for bar in bars:
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height + 1, f"{height:.1f}", 
                        ha='center', va='bottom', fontsize=9, fontweight='bold')
            
            ax2.set_ylabel("Average pLDDT", fontsize=9)
            plt.xticks(rotation=45, fontsize=9)
            ax2.grid(axis='y', linestyle='--', alpha=0.3)
            
            buf2 = io.BytesIO()
            fig2.savefig(buf2, format="png", bbox_inches='tight', dpi=100)
            plt.close(fig2)
            buf2.seek(0)
            img2_b64 = base64.b64encode(buf2.read()).decode('utf-8')
            
            st.markdown(f"""
<div class="sci-table-card" style="height:100%; margin-bottom:0;">
    <div class="sci-table-header">
        <span style="font-size:1.2rem; display: flex; align-items: center;">📈</span>
        <h4 style="margin:0;font-size:0.95rem;font-weight:700;color:#1e293b; position: relative; top: 5px; line-height: 1;">예측 모델별 평균 pLDDT 점수 (Average pLDDT Score)</h4>
    </div>
    <div class="sci-table-body" style="text-align:center;">
        <img src="data:image/png;base64,{img2_b64}" style="max-width:100%; border-radius:0.25rem;" />
    </div>
</div>
""", unsafe_allow_html=True)

    # --- pLDDT Explanation Guide ---
    st.markdown("""
<div style="background-color: #eff6ff; border: 1px solid #bfdbfe; border-radius: 0.5rem; padding: 1.25rem; margin-top: 1rem; margin-bottom: 2rem;">
    <div style="display: flex; align-items: flex-start; gap: 0.75rem;">
        <span class="material-symbols-outlined" style="color: #2563eb; font-size: 1.4rem; font-variation-settings: 'FILL' 1;">info</span>
        <div style="flex: 1;">
            <p style="margin: 0 0 0.6rem 0; font-size: 0.9rem; font-weight: 700; color: #1e40af;">pLDDT 점수에 대한 설명 및 해석 방법</p>
            <p style="margin: 0 0 1rem 0; font-size: 0.85rem; color: #1e3a8a; line-height: 1.6; word-break: keep-all;">
                <b>pLDDT 점수:</b> 잔기별 신뢰도 점수, 예측된 구조의 각 아미노산 잔기가 얼마나 정확하게 예측되었는지를 0에서 100 사이의 점수로 나타내는 신뢰도 지표
            </p>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 0.75rem; border-top: 1px dashed #bfdbfe; padding-top: 1rem;">
                <div style="display: flex; align-items: center; gap: 0.6rem;">
                    <div style="width: 12px; height: 12px; background-color: #0053D6; border-radius: 2px; flex-shrink: 0;"></div>
                    <span style="font-size: 0.85rem; color: #1e3a8a;"><b>90 이상:</b> 높은 정확도</span>
                </div>
                <div style="display: flex; align-items: center; gap: 0.6rem;">
                    <div style="width: 12px; height: 12px; background-color: #65CBF3; border-radius: 2px; flex-shrink: 0;"></div>
                    <span style="font-size: 0.85rem; color: #1e3a8a;"><b>70 이상 90 미만:</b> 전반적으로 양호한 정확도</span>
                </div>
                <div style="display: flex; align-items: center; gap: 0.6rem;">
                    <div style="width: 12px; height: 12px; background-color: #FFE915; border-radius: 2px; flex-shrink: 0;"></div>
                    <span style="font-size: 0.85rem; color: #1e3a8a;"><b>50 이상 70 미만:</b> 낮은 신뢰도</span>
                </div>
                <div style="display: flex; align-items: center; gap: 0.6rem;">
                    <div style="width: 12px; height: 12px; background-color: #FF7D45; border-radius: 2px; flex-shrink: 0;"></div>
                    <span style="font-size: 0.85rem; color: #1e3a8a;"><b>50 미만:</b> 무질서 또는 비구조적 상태</span>
                </div>
            </div>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

    # [Section 3] Structure Downloads
    dl_title = "항생제 내성 유발 단백질의 3차원 구조 예측 모델 다운로드 (.pdb)" if is_boards else "해당 단백질의 3차원 구조 예측 모델 다운로드 (.pdb)"
    st.markdown(f"""
<div style="margin:2.5rem 0 1.25rem 0;display:flex;align-items:center;gap:0.6rem;">
    <span style="font-size:1.4rem; display: flex; align-items: center;">📁</span>
    <h3 style="margin:0;font-size:1.15rem;font-weight:700;color:#1e293b; position: relative; top: 2px;">{dl_title}</h3>
</div>
""", unsafe_allow_html=True)
    
    dl_cols = st.columns(max(1, len(models)))
    for i, m in enumerate(models):
        with dl_cols[i]:
            actual_rank = m.get('rank', i)
            if os.path.exists(m.get('pdb_path', '')):
                st.markdown(f'<div style="font-size:0.75rem;font-weight:700;color:#334155;margin-bottom:0.3rem;">RANK {actual_rank}</div>', unsafe_allow_html=True)
                # Using the original filename from BOARDS
                orig_name = f"index_1980_rank{actual_rank}.pdb" if "1980" in m.get('file_name','') else m.get('file_name', f"ranked_{actual_rank}.pdb")
                # `use_container_width` is being deprecated in favor of `width`
                st.download_button("📄 PDB", open(m['pdb_path']).read(), file_name=orig_name, key=f"dl_dash_{actual_rank}_{i}", use_container_width=True)
            else:
                st.markdown(f'<div style="font-size:0.75rem;font-weight:700;color:#334155;margin-bottom:0.3rem;">RANK {actual_rank}</div>', unsafe_allow_html=True)
                st.write("❌ 없음")

def render_model_selection_and_analysis_options(models, is_boards=True):
    """Encapsulates the model selection, 3D viz, and analysis settings."""
    import shutil
    
    st.markdown("""
    <style>
        /* NO DANGEROUS :has() selectors here anymore to prevent bubbling */
        
        /* Make Expanders white */
        div[data-testid="stExpander"], 
        div[data-testid="stExpander"] > details {
            background-color: white !important;
            background: white !important;
            border-color: #e2e8f0 !important;
        }
        /* Ensure selectboxes have a solid grey background for visibility */
        div[data-testid="stSelectbox"] div[data-baseweb="select"] > div {
            background-color: #e2e8f0 !important;
            border: 1px solid #cbd5e1 !important;
            border-radius: 0.35rem !important;
            cursor: pointer !important;
            color: #1e293b !important;
        }
        /* Standard secondary buttons styling (PDB buttons) */
        button[kind="secondary"] {
            background-color: white !important;
            border: 1px solid #e2e8f0 !important;
            color: #334155 !important;
        }
    </style>
    """, unsafe_allow_html=True)
    
    # Bulletproof JS Injection using Marker pattern
    components.html("""
    <script>
    (function() {
        function applyWhite() {
            // Find our specific inserted flag
            const flags = window.parent.document.querySelectorAll('.viz-bg-flag');
            flags.forEach(function(flag) {
                // Find nearest Streamlit container block
                let block = flag.closest('[data-testid="stVerticalBlock"]');
                if (block) {
                    block.style.setProperty('background-color', 'white', 'important');
                    block.style.setProperty('background', 'white', 'important');
                    if (block.parentElement) {
                        block.parentElement.style.setProperty('background-color', 'white', 'important');
                        block.parentElement.style.setProperty('background', 'white', 'important');
                    }
                }
                let wrapper = flag.closest('[data-testid="stVerticalBlockBorderWrapper"]');
                if (wrapper) {
                    wrapper.style.setProperty('background-color', 'white', 'important');
                    wrapper.style.setProperty('background', 'white', 'important');
                    const cacheKids = wrapper.querySelectorAll('div[class^="st-emotion-cache"]');
                    cacheKids.forEach(k => {
                        k.style.setProperty('background-color', 'white', 'important');
                    });
                }
            });
        }
        applyWhite();
        setTimeout(applyWhite, 100);
        setTimeout(applyWhite, 500);
        setTimeout(applyWhite, 1000);
        const observer = new MutationObserver(applyWhite);
        observer.observe(window.parent.document.body, { childList: true, subtree: true });
    })();
    </script>
    """, height=0)
    
    with st.container(border=True):
        st.markdown("<span class='viz-bg-flag'></span>", unsafe_allow_html=True)
        # 1. Header Row - Standardized height and font sizes
        # Using [3.5, 1, 1.2] to give enough room and keep elements close
        hdr_col1, hdr_col2, hdr_col3 = st.columns([3.5, 1, 1.2], vertical_alignment="center")
        
        with hdr_col1:
            sel_title = "항생제 내성 유발 단백질의 3차원 구조 예측 모델 선택 및 시각화" if is_boards else "해당 단백질의 3차원 구조 예측 모델 선택 및 시각화"
            st.markdown(f"""
            <div style="display:flex; align-items:center; gap:0.6rem;">
                <span class="material-symbols-outlined" style="font-size:1.5rem; color:#ef4444; position: relative; top: 1px;">visibility</span>
                <h3 style="margin:0; font-size:1.15rem; font-weight:700; color:#1e293b; white-space:nowrap; position: relative; top: 1px;">{sel_title}</h3>
            </div>
            """, unsafe_allow_html=True)
            
        with hdr_col2:
            st.markdown("<div style='font-size:1.1rem; font-weight:700; color:#1e293b; text-align:right; white-space:nowrap; position: relative; top: 1px;'>예측 모델 선택:</div>", unsafe_allow_html=True)
            
        with hdr_col3:
            # Pushing the selectbox down slightly to match the text height
            st.markdown("<div style='margin-top: 5px;'></div>", unsafe_allow_html=True)
            selected_rank = st.selectbox(
                "예측 모델 선택", 
                [m['rank'] for m in models], 
                format_func=lambda r: f"Rank {r} (pLDDT: {next(m.get('plddt', 0) for m in models if m['rank'] == r):.1f})",
                index=0,
                key=f"select_rank_{st.session_state['pipeline_step']}",
                label_visibility="collapsed"
            )
        
        selected_model = next(m for m in models if m['rank'] == selected_rank)
        
        # 2. 3D View (pLDDT coloring)
        view = py3Dmol.view(width=860, height=500)
        pdb_render_str = selected_model.get('pdb_data') or (open(selected_model['pdb_path'],'r').read() if selected_model.get('pdb_path') else "")
        if pdb_render_str:
            view.addModel(pdb_render_str, 'pdb')
            
            # Base scheme
            view.setStyle({'cartoon': {'color': 'white'}})
            
            # Exactly mapping pLDDT thresholds
            p_list = selected_model.get('plddt_list', [])
            blue_resi = [i+1 for i, p in enumerate(p_list) if p >= 90]
            cyan_resi = [i+1 for i, p in enumerate(p_list) if 70 <= p < 90]
            yellow_resi = [i+1 for i, p in enumerate(p_list) if 50 <= p < 70]
            orange_resi = [i+1 for i, p in enumerate(p_list) if p < 50]
    
            if blue_resi: view.addStyle({'resi': blue_resi}, {'cartoon': {'color': '#0053D6'}})
            if cyan_resi: view.addStyle({'resi': cyan_resi}, {'cartoon': {'color': '#65CBF3'}})
            if yellow_resi: view.addStyle({'resi': yellow_resi}, {'cartoon': {'color': '#FFE915'}})
            if orange_resi: view.addStyle({'resi': orange_resi}, {'cartoon': {'color': '#FF7D45'}})
            
            hover_js = """
            function(atom,viewer) {
                if(!atom.label) {
                    atom.label = viewer.addLabel(atom.resn + atom.resi + " (pLDDT: " + atom.b.toFixed(1) + ")",
                    {position: atom, backgroundColor: '#202020', fontColor: 'white', fontSize: 14});
                }
            }"""
            unhover_js = """
            function(atom,viewer) {
                if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                }
            }"""
            view.setHoverable({}, True, hover_js, unhover_js)
            view.setBackgroundColor('white')
            view.zoomTo()
        showmol(view, height=500, width=860)
        
        st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
        
        # 3. Legend Row (Two columns)
        col_leg1, col_leg2 = st.columns([1, 1])
        with col_leg1:
            st.markdown("""
            <div style="border: 1px solid #e2e8f0; border-radius: 0.25rem; padding: 1rem 1.25rem; background-color: #f8fafc; height: 100%;">
                <div style="font-weight: 700; font-size:0.75rem; color: #475569; letter-spacing: 0.05em; margin-bottom: 0.75rem;">MODEL CONFIDENCE</div>
                <div style="display: flex; flex-wrap: wrap; gap: 0.8rem 1rem; font-size: 0.75rem; color: #334155;">
                    <div style="display:flex; align-items:center; gap:0.4rem; width: 45%;">
                        <div style="width: 14px; height: 14px; background-color: #0053D6; border-radius: 2px;"></div> Very high (pLDDT > 90)
                    </div>
                    <div style="display:flex; align-items:center; gap:0.4rem; width: 45%;">
                        <div style="width: 14px; height: 14px; background-color: #65CBF3; border-radius: 2px;"></div> Confident (90 > pLDDT > 70)
                    </div>
                    <div style="display:flex; align-items:center; gap:0.4rem; width: 45%;">
                        <div style="width: 14px; height: 14px; background-color: #FFE915; border-radius: 2px;"></div> Low (70 > pLDDT > 50)
                    </div>
                    <div style="display:flex; align-items:center; gap:0.4rem; width: 45%;">
                        <div style="width: 14px; height: 14px; background-color: #FF7D45; border-radius: 2px;"></div> Very low (pLDDT < 50)
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
        with col_leg2:
            st.markdown("""
            <div style="border: 1px solid #bfdbfe; border-radius: 0.25rem; padding: 1.5rem 1.5rem; background-color: #eff6ff; display: flex; align-items: center; height: 100%;">
                <div style="font-size: 0.825rem; color: #1e40af; line-height: 1.5;">AlphaFold2는 각 아미노산 잔기마다 0에서 100 사이의 값을 갖는 신뢰도 점수(pLDDT)를 산출한다.</div>
            </div>
            """, unsafe_allow_html=True)
            
        st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)
        
        # 4. Expandable section inside the border
        exp_title = "📁 선택된 예측 모델 다운로드 (.pdb)"
        with st.expander(exp_title, expanded=True):
            exp_col1, exp_col2 = st.columns([2.5, 1], vertical_alignment="top")
            with exp_col1:
                plddt_val = selected_model.get('plddt', 0)
                st.markdown(f"""
                <div style="display:flex; align-items:flex-start; gap:1rem; padding-top: 4px;">
                    <div style="width:2.5rem; height:2.5rem; background-color:#eff6ff; border-radius:50%; display:flex; align-items:center; justify-content:center; color:#2563eb;">
                        <span class="material-symbols-outlined" style="font-size:1.4rem;">description</span>
                    </div>
                    <div style="line-height:1.2;">
                        <div style="font-weight:700; font-size:0.9rem; color:#1e293b;">현재 선택된 모델: Rank {selected_rank}</div>
                        <div style="font-size:0.75rem; color:#64748b; margin-top:0.25rem;">pLDDT 점수: {plddt_val:.1f}%</div>
                    </div>
                </div>
                """, unsafe_allow_html=True)
            
            with exp_col2:
                if os.path.exists(selected_model.get('pdb_path', '')):
                    with open(selected_model['pdb_path'], 'rb') as f:
                        st.download_button(
                            label="📕 PDB 다운로드",
                            data=f.read(),
                            file_name=os.path.basename(selected_model['pdb_path']),
                            key=f"dl_btn_{selected_rank}_pdb_exp",
                            use_container_width=True
                        )

# --- [Result History & Recovery Helpers] ---
def parse_existing_af2_results(af_out_abs):
    """지정한 결과 폴더에서 PDB/JSON 데이터를 파싱하여 모델 리스트를 반환합니다."""
    import glob
    import re
    
    pdb_files = sorted(glob.glob(os.path.join(af_out_abs, "*.pdb")))
    if not pdb_files: return []

    def extract_rank_from_pdb(f):
        m = re.search(r'rank_?(\d+)', f)
        if m: return int(m.group(1))
        m = re.search(r'model_(\d+)', f)
        if m: return int(m.group(1))
        return 999
    pdb_files.sort(key=extract_rank_from_pdb)
    
    models_data = []
    all_json_files = glob.glob(os.path.join(af_out_abs, "*.json"))
    pipe = FinalWorkflowPipeline(init_db())
    
    for i, pdb_f in enumerate(pdb_files):
        rank_id = str(extract_rank_from_pdb(pdb_f))
        plddt_list = pipe.get_plddt_stats_fixed(pdb_f)
        avg_plddt = np.mean(plddt_list) if plddt_list else 0.0
        with open(pdb_f, "r") as f: pdb_content = f.read()
        
        pae_data = []
        potential_jsons = [j for j in all_json_files if f"model_{rank_id}" in j or f"rank_{rank_id}" in j]
        if not potential_jsons:
             potential_jsons = [j for j in all_json_files if f"_{i}." in j or f"_{i}_" in j]
        if potential_jsons:
            try:
                with open(potential_jsons[0], 'r') as jf:
                    js = json.load(jf)
                    if isinstance(js, list): js = js[0]
                    pae_data = js.get('pae') or js.get('predicted_aligned_error') or []
            except: pass
        
        models_data.append({
            "rank": i, "display_rank": int(rank_id) if rank_id != "999" else i,
            "file_name": os.path.basename(pdb_f), "pdb_path": pdb_f,
            "plddt": avg_plddt, "plddt_list": plddt_list, "pae_data": pae_data, "pdb_data": pdb_content
        })
    # [ADD] 공식 PAE 그리드 이미지 탐색 (*_pae.png)
    pae_image_path = None
    png_files = glob.glob(os.path.join(af_out_abs, "*_pae.png"))
    if png_files:
        pae_image_path = png_files[0]

    return {
        "models": models_data,
        "pae_image": pae_image_path
    }


# --- Main App ---
st.set_page_config(page_title="Final Workflow Pipeline Dashboard", layout="wide")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@500&family=Inter:wght@400;600&display=swap');
@import url('https://fonts.googleapis.com/css2?family=Material+Symbols+Outlined:opsz,wght,FILL,GRAD@20..48,100..700,0..1,-50..200&display=swap');


/* Typography Global */
h1, h2, h3, h4, h5 {
    font-family: 'Space Grotesk', sans-serif !important;
}
div, p, span {
    font-family: 'Inter', sans-serif;
}

/* Background & Containers */
div[data-testid="stAppViewContainer"] {
    background-color: #f3f4f6;
    color: #111827;
}
/* In Streamlit 1.55.0, bordered containers use class stVerticalBlock.
   Target via the secondaryBackgroundColor (set in config.toml to #ffffff) */
.stVerticalBlock {
    background-color: transparent;
}

/* Base Input Elements Override (Text box, textarea, number input) */
div[data-baseweb="input"] > div, 
div[data-baseweb="textarea"] > div, 
div[data-baseweb="base-input"] {
    background-color: #eef0f3 !important;
    border-color: #d1d5db !important;
    border-radius: 0.5rem !important;
}
div[data-baseweb="input"] > div:focus-within, 
div[data-baseweb="textarea"] > div:focus-within {
    border-color: #2563eb !important;
    box-shadow: 0 0 0 1px #2563eb !important;
}

/* Compact Selectbox for Dashboard alignment */
div[data-testid="stSelectbox"] {
    transform: translateY(1px); /* Minor nudge for baseline alignment in center mode */
}
div[data-testid="stSelectbox"] div[data-baseweb="select"] > div {
    min-height: 38px !important;
    height: 38px !important;
    background-color: #f8fafc !important;
    border-radius: 0.5rem !important;
    padding: 0 12px !important;
    display: flex !important;
    align-items: center !important;
}
div[data-testid="stSelectbox"] [data-baseweb="select"] {
    font-size: 0.95rem !important;
    font-weight: 500 !important;
}

/* Download Grid Styling - Matching Mockup */
.rank-grid-label {
    font-size: 0.7rem !important;
    font-weight: 800 !important;
    color: #475569 !important;
    text-transform: uppercase !important;
    letter-spacing: 0.05em !important;
    margin-bottom: 0.4rem !important;
}

div[data-testid="stExpander"] div[data-testid="stDownloadButton"] {
    transform: translateY(0) !important; 
}

div[data-testid="stExpander"] div[data-testid="stDownloadButton"] button {
    background-color: #ffffff !important;
    color: #1e293b !important;
    border: 1px solid #e2e8f0 !important;
    border-radius: 4px !important;
    font-weight: 600 !important;
    font-size: 0.85rem !important;
    height: 42px !important;
    box-shadow: none !important;
    transition: all 0.2s ease !important;
}

div[data-testid="stExpander"] div[data-testid="stDownloadButton"] button:hover {
    border-color: #cbd5e1 !important;
    background-color: #f8fafc !important;
}

/* File uploader drop zone - more visible grey */
[data-testid="stFileUploaderDropzone"],
section[data-testid="stFileUploaderDropzone"] {
    background-color: #eef0f3 !important;
    border: 1.5px dashed #9ca3af !important;
    border-radius: 0.5rem !important;
}


/* Streamlit Button (Primary) */
button[kind="primary"] {
    background: #005AC2 !important;
    color: white !important;
    border-radius: 0.5rem !important;
    padding: 0.65rem 0rem !important;
    font-size: 1.0rem !important;
    font-weight: 600 !important;
    border: none !important;
    box-shadow: 0 1px 2px rgba(0, 0, 0, 0.05) !important;
    transition: all 0.2s ease-in-out !important;
}
button[kind="primary"]:hover {
    background: #004494 !important;
    transform: translateY(-1px);
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1) !important;
}
button[kind="primary"]:active {
    transform: translateY(0);
}

/* Sidebar history buttons: make it look like a nav item */
[data-testid="stSidebarContent"] button:not([kind="primary"]),
[data-testid="stSidebar"] button[kind="secondary"] {
    background-color: transparent !important;
    color: #64748b !important;
    border: none !important;
    border-radius: 0.5rem !important;
    box-shadow: none !important;
    font-size: 0.9rem !important;
    font-weight: 600 !important;
    padding: 0.65rem 1.25rem !important;
    display: flex !important;
    justify-content: flex-start !important;
}
[data-testid="stSidebarContent"] button:not([kind="primary"]) div,
[data-testid="stSidebarContent"] button:not([kind="primary"]) p,
[data-testid="stSidebar"] button[kind="secondary"] div,
[data-testid="stSidebar"] button[kind="secondary"] p {
    text-align: left !important;
    justify-content: flex-start !important;
    display: flex !important;
    width: 100% !important;
    margin: 0 !important;
}
[data-testid="stSidebarContent"] button:not([kind="primary"]):hover,
[data-testid="stSidebar"] button[kind="secondary"]:hover {
    background-color: #f1f5f9 !important;
    color: #1d4ed8 !important;
}

/* Sidebar Background & Border Enhancement */
[data-testid="stSidebar"] {
    background-color: #e5e7eb !important;
    border-right: 1px solid #d1d5db !important;
}

/* Scientific card & table styling */
.sci-table-card { 
    background: white; 
    border: 1px solid #e2e8f0; 
    border-radius: 0.25rem; 
    overflow: hidden; 
    box-shadow: 0 1px 2px rgba(0,0,0,0.05); 
    margin-bottom: 1.5rem; 
}
.sci-table-header { 
    background: #f8fafc; 
    padding: 0.85rem 1.25rem; 
    font-size: 0.9rem; 
    font-weight: 700; 
    color: #1e293b; 
    border-bottom: 1px solid #e2e8f0; 
    display: flex;
    align-items: center;
    gap: 0.5rem;
}
.sci-table-header h4 {
    line-height: 1;
    position: relative;
    top: 5px; /* Offset to visually align with icons */
}
.sci-table-body {
    padding: 1.25rem;
}

.exec-table { width: 100%; border-collapse: collapse; font-size: 0.8125rem; }
.exec-table th { background: #f1f5f9; font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.05em; color: #475569; padding: 0.65rem 1rem; border-bottom: 2px solid #e2e8f0; text-align: left; }
.exec-table td { padding: 0.65rem 1rem; border-bottom: 1px solid #f1f5f9; color: #334155; font-family: monospace; }
.exec-table tr:last-child td { border-bottom: none; }
.status-ok { color: #22c55e; font-weight: 700; font-family: inherit; }

.back-btn-wrap a, .back-btn-wrap button { 
    display: inline-flex; 
    align-items: center; 
    gap: 0.5rem; 
    padding: 0.5rem 1rem; 
    background: white; 
    border: 1px solid #e2e8f0; 
    border-radius: 9999px; 
    font-size: 0.875rem; 
    font-weight: 700; 
    color: #374151; 
    box-shadow: 0 1px 2px rgba(0,0,0,0.05); 
    cursor: pointer; 
    text-decoration: none; 
}

.eval-table { width: 100%; border-collapse: collapse; font-size: 0.8125rem; }
.eval-table th { background: rgba(239,246,255,0.3); font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.05em; color: #1d4ed8; padding: 0.65rem 1rem; border-bottom: 2px solid #bfdbfe; text-align: left; }
.eval-table td { padding: 0.65rem 1rem; border-bottom: 1px solid #eff6ff; color: #475569; }
.eval-table tr:last-child td { border-bottom: none; }
.eval-mono { color: #2563eb; font-family: monospace; font-weight: 700; }

/* Force white background for bordered containers */
div[data-testid="stVerticalBlockBorderWrapper"],
div[data-testid="stVerticalBlockBorderWrapper"] > div {
    background-color: white !important;
}

/* Red home button in sidebar - Natural bottom placement */
.st-key-back_to_main_pselector {
    margin-top: 0.4rem !important;
    padding-bottom: 1.5rem !important;
}

.st-key-back_to_main_pselector button {
    background-color: #ef4444 !important;
    border-color: #ef4444 !important;
    color: white !important;
    font-weight: 700 !important;
    border-radius: 0.35rem !important;
    height: 3rem !important;
    font-size: 1.05rem !important;
    transition: all 0.2s !important;
}
.st-key-back_to_main_pselector button:hover {
    background-color: #dc2626 !important;
    border-color: #dc2626 !important;
}

</style>
</style>
""", unsafe_allow_html=True)

# 사이드바 렌더링
with st.sidebar:
    st.markdown("""
    <div style='padding-top: 0.5rem; padding-bottom: 1.5rem; padding-left: 0.5rem;'>
        <div style='font-size: 1.5rem; font-weight: 800; color: #1d4ed8; margin-bottom: 0.2rem; letter-spacing: -0.02em;'>AMR Monitor</div>
        <div style='font-size: 0.8rem; color: #64748b; font-weight: 500;'>by ResistBreakers</div>
    </div>
    """, unsafe_allow_html=True)
    
    step_now = st.session_state.get("pipeline_step", "INTAKE")
    
    def render_nav_item(label, icon, is_active):
        bg = "#eff6ff" if is_active else "transparent"
        color = "#1d4ed8" if is_active else "#64748b"
        border = "4px solid #1d4ed8" if is_active else "4px solid transparent"
        font_weight = "700" if is_active else "500"
        icon_style = "font-variation-settings: 'FILL' 1;" if is_active else ""
        
        st.markdown(f"""
        <div style="
            display: flex; align-items: center; gap: 0.8rem;
            padding: 0.75rem 1rem;
            margin-bottom: 0.25rem;
            background-color: {bg};
            border-left: {border};
            color: {color};
            border-radius: 0 0.5rem 0.5rem 0;
            font-size: 0.9rem;
            font-weight: {font_weight};
        ">
            <span class="material-symbols-outlined" style="font-size: 1.25rem; {icon_style}">{icon}</span>
            <span>{label}</span>
        </div>
        """, unsafe_allow_html=True)
        
    render_nav_item("DNA 또는 아미노산 서열 입력", "upload_file", step_now == "INTAKE")
    render_nav_item("분석 도구 선택", "route", step_now == "ROUTE_SELECT")
    render_nav_item("항생제 내성 유전자의 단백질 구조 예측 정보 (BOARDS DB)", "database", step_now == "BOARDS_HIT")
    render_nav_item("AlphaFold2 결과", "view_in_ar", step_now in ["AF2_RUNNING", "AF2_RESULTS"])
    
    st.markdown("<hr style='margin:1.5rem 0 0.5rem 0; border:none; border-top:1px solid #e2e8f0; opacity:0.5;'>", unsafe_allow_html=True)
    
    # 파일 보관함: 단일 버튼으로 전체 페이지 이동
    sb_hist_icon = "📂" if step_now == "HISTORY" else "📁"
    if st.button(f"{sb_hist_icon} 파일 보관함", use_container_width=True, key="btn_open_history"):
        st.session_state["pipeline_step"] = "HISTORY"
        st.rerun()
        
    st.markdown("<div style='margin: 0.4rem 0;'></div>", unsafe_allow_html=True)
    if st.button("새 프로젝트", type="primary", use_container_width=True, key="new_project_sidebar"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()

    # --- Home Return Button (Placeholder for multi-app integration) ---
    if st.button("🏠 서비스 홈", type="primary", use_container_width=True, key="back_to_main_pselector"):
        # This will be integrated with the main selection screen later
        pass


# Custom Header matching 1_start.png
st.markdown("""
<div style="display: flex; align-items: center; margin-bottom: 2.5rem; margin-top: 1rem;">
    <h1 style="margin: 0; padding: 0; font-size: 1.5rem; font-weight: 800; color: #111827; letter-spacing: -0.02em; white-space: nowrap;">AMR Monitor</h1>
    <span style="color: #d1d5db; margin: 0 1rem; font-size: 1.4rem;">|</span>
    <span style="color: #2563eb; font-size: 1.05rem; font-weight: 600; letter-spacing: -0.01em;">항생제 내성 유전자 탐지 및 항생제 내성 유발 단백질의 3차원 구조 예측 모델 생성</span>
</div>
""", unsafe_allow_html=True)

# Top Reset Button
if st.session_state.get("pipeline_step", "INTAKE") != "INTAKE":
    if st.button("🔄 처음으로 돌아가기", key="reset_top"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()

# 프로그램을 실행 정보 조건부 노출 (파일 보관함 및 히스토리 뷰에서는 숨김)
_hide_exec_log = st.session_state.get("is_from_history", False) or st.session_state.get("pipeline_step") == "HISTORY"
if not _hide_exec_log:
    render_execution_log(st.session_state.get("execution_log", []))

# ====== MAIN PIPELINE FLOW ======

if "pipeline_step" not in st.session_state:
    st.session_state["pipeline_step"] = "INTAKE"

# --- [파일 보관함] 전체 화면 목록 페이지 ---
if st.session_state["pipeline_step"] == "HISTORY":
    st.markdown("""
<div style="display:flex; align-items:center; gap:0.75rem; margin-bottom:2rem;">
    <span class="material-symbols-outlined" style="color:#2563eb; font-size:2rem;">folder_open</span>
    <h2 style="margin:0; font-size:1.4rem; font-weight:700; color:#111827;">&#xD30C;&#xC77C; &#xBCF4;&#xAD00;&#xD568;</h2>
</div>
""", unsafe_allow_html=True)
    
    results_dir = "results"
    all_dirs = []
    if os.path.exists(results_dir):
        all_dirs = sorted(
            [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d)) and "RUN_RADAR" in d],
            reverse=True
        )
    
    if not all_dirs:
        st.info("저장된 결과 파일이 없습니다. 새 프로젝트를 시작해 주세요.")
    else:
        # 헤더 행
        # 헤더 행
        hdrc1, hdrc2, hdrc3, hdrc4 = st.columns([3.4, 2, 0.6, 0.8])
        with hdrc1:
            st.markdown("<div style='font-size:0.75rem; font-weight:600; color:#6b7280; text-transform:uppercase; letter-spacing:0.05em; padding:0.5rem 0;'>&#xD30C;&#xC77C;&#xBA85;</div>", unsafe_allow_html=True)
        with hdrc2:
            st.markdown("<div style='font-size:0.75rem; font-weight:600; color:#6b7280; text-transform:uppercase; letter-spacing:0.05em; padding:0.5rem 0;'>&#xC0DD;&#xC131;&#xC77C;</div>", unsafe_allow_html=True)
        with hdrc3:
            st.markdown("<div style='font-size:0.75rem; font-weight:600; color:#6b7280; text-transform:uppercase; letter-spacing:0.05em; padding:0.5rem 0; text-align:center;'>&#xACB0;&#xACFC;</div>", unsafe_allow_html=True)
        with hdrc4:
            st.markdown("<div style='font-size:0.75rem; font-weight:600; color:#6b7280; text-transform:uppercase; letter-spacing:0.05em; padding:0.5rem 0; text-align:center;'>&#xAD00;&#xB9AC;</div>", unsafe_allow_html=True)
        
        st.markdown("<hr style='margin:0 0 0.5rem 0; border-color:#e5e7eb;'>", unsafe_allow_html=True)
        
        for d in all_dirs:
            display_time = d.replace("RUN_RADAR_", "")
            try:
                dt_obj = datetime.strptime(display_time, "%Y%m%d%H%M")
                created_str = dt_obj.strftime("%Y.%m.%d %H:%M")
                folder_name = dt_obj.strftime("%y년%m월%d일 %H:%M")
            except:
                created_str = display_time
                folder_name = d
            
            jdir = os.path.join(results_dir, d)
            # 폴더 크기 계산
            total_size = 0
            try:
                for dirpath, dirnames, filenames in os.walk(jdir):
                    for fname in filenames:
                        fp = os.path.join(dirpath, fname)
                        total_size += os.path.getsize(fp)
                size_mb = total_size / (1024 * 1024)
                size_str = f"{size_mb:.1f} MB" if size_mb >= 1 else f"{total_size / 1024:.0f} KB"
            except:
                size_str = "-"
            
            # PDB 파일 유무로 상태 판단
            import glob
            pdb_list = glob.glob(os.path.join(jdir, "**/*.pdb"), recursive=True)
            status_html = "<span style='display:inline-flex;align-items:center;gap:0.35rem;font-size:0.8rem;color:#059669;'><span style='width:8px;height:8px;border-radius:50%;background:#059669;display:inline-block;'></span>완료</span>" if pdb_list else "<span style='display:inline-flex;align-items:center;gap:0.35rem;font-size:0.8rem;color:#d97706;'><span style='width:8px;height:8px;border-radius:50%;background:#d97706;display:inline-block;'></span>보관</span>"
            
            row_c1, row_c2, row_c3, row_c4 = st.columns([3.4, 2, 0.6, 0.8], vertical_alignment="center")
            
            with row_c1:
                # 파일 정보 (파일명 제목 아래)
                st.markdown(f"""
<div style="display:flex;align-items:center;gap:0.75rem;padding:0.6rem 0;">
    <div style="width:2.2rem;height:2.2rem;background-color:#eff6ff;border-radius:0.4rem;display:flex;align-items:center;justify-content:center;flex-shrink:0;">
        <span class="material-symbols-outlined" style="font-size:1.2rem;color:#2563eb;">description</span>
    </div>
    <div>
        <div style="font-weight:600;font-size:0.9rem;color:#111827;">{folder_name}</div>
        <div style="font-size:0.75rem;color:#6b7280;">AlphaFold2 분석 결과 &bull; {size_str}</div>
    </div>
</div>""", unsafe_allow_html=True)
                
            with row_c2:
                # 생성일
                st.markdown(f"<div style='padding:0.75rem 0; font-size:0.875rem; color:#4b5563;'>{created_str}</div>", unsafe_allow_html=True)
            
            with row_c3:
                # 로드 버튼 (결과 제목 아래)
                if st.button("📂", key=f"hist_load_{d}", help="이 결과 열기", use_container_width=True):
                    afdirs = [ad for ad in os.listdir(jdir) if os.path.isdir(os.path.join(jdir, ad))]
                    af_out = os.path.join(jdir, afdirs[0]) if afdirs else jdir
                    recovered = parse_existing_af2_results(af_out)
                    if recovered and recovered.get("models"):
                        st.session_state["af2_models"] = recovered
                        st.session_state["job_dir"] = jdir
                        st.session_state["pipeline_step"] = "AF2_RESULTS"
                        st.session_state["is_from_history"] = True
                        st.rerun()
                    else:
                        st.error("결과를 찾을 수 없습니다.")
                        
            with row_c4:
                # 삭제 버튼 (관리 제목 아래)
                if st.button("🗑️", key=f"hist_del_{d}", help="삭제", use_container_width=True):
                    shutil.rmtree(jdir, ignore_errors=True)
                    st.rerun()
            
            st.markdown("<hr style='margin:0; border-color:#f3f4f6;'>", unsafe_allow_html=True)
    
    st.stop()

# --- [입력] WGS FASTA (.fna) 타겟 서열 업로드 ---
if st.session_state["pipeline_step"] == "INTAKE":
    # JS injection: find bordered containers via window.parent.document and force white background
    components.html("""
<script>
(function() {
    function applyWhite() {
        try {
            var doc = window.parent.document;
            // Find all divs and check if they have a visible border (the st.container(border=True) elements)
            var divs = doc.querySelectorAll('div[data-testid="stVerticalBlock"]');
            divs.forEach(function(el) {
                var style = window.parent.getComputedStyle(el);
                // Bordered containers have a non-transparent, non-none border
                if (style.borderTopWidth !== '0px' && style.borderTopWidth !== '') {
                    el.style.setProperty('background-color', '#ffffff', 'important');
                }
            });
        } catch(e) {}
    }
    // Apply immediately and after delays to handle React render cycles
    applyWhite();
    setTimeout(applyWhite, 300);
    setTimeout(applyWhite, 800);
    setTimeout(applyWhite, 2000);
    // Watch for DOM changes
    try {
        var observer = new MutationObserver(function() { applyWhite(); });
        observer.observe(window.parent.document.body, {childList: true, subtree: true});
    } catch(e) {}
})();
</script>
""", height=0)
    st.markdown("<br>", unsafe_allow_html=True)
    col_a, col_b = st.columns(2, gap="large")
    
    with col_a:
        with st.container(border=True):
            st.markdown("""
<div style="display: flex; align-items: center; gap: 0.5rem; margin-bottom: 0;">
    <span class="material-symbols-outlined" style="color: #2563eb; font-size: 1.25rem; font-variation-settings: 'FILL' 0;">terminal</span>
    <h4 style="margin: 0; font-size: 1.1rem; font-weight: 700; color: #111827; font-family: 'Space Grotesk', sans-serif;">DNA 또는 아미노산 서열 입력</h4>
</div>
""", unsafe_allow_html=True)
            st.markdown("<div style='margin-bottom: 0.8rem;'></div>", unsafe_allow_html=True)
            radar_input_method_label = st.segmented_control("입력 방식 선택", ["파일(File)", "텍스트(Text)"], default="파일(File)", label_visibility="collapsed")
            
            radar_input_method = "Text_Input" if radar_input_method_label == "텍스트(Text)" else "File_Upload"
            
            st.markdown("<div style='margin-top: 1rem;'></div>", unsafe_allow_html=True)
            
            uploaded_file = None
            pasted_seq = ""
            
            # 선택한 라디오 값에 따라 해당 UI만 표시
            if radar_input_method == "File_Upload":
                uploaded_file = st.file_uploader("여기로 파일을 끌고 오거나 클릭하여 업로드하세요.\n\n(지원 파일 형식: .fna, .faa, .fasta)", type=["fna", "faa", "fasta"])
            else:
                pasted_seq = st.text_area("DNA 또는 아미노산 서열을 입력하세요.", height=192, placeholder="MNIVEN...")

    with col_b:
        with st.container(border=True):
            st.markdown("""
<div style="display: flex; align-items: center; gap: 0.5rem; margin-bottom: 1.5rem;">
    <span class="material-symbols-outlined" style="color: #2563eb; font-size: 1.3rem;">podcasts</span>
    <h4 style="margin: 0; font-size: 1.1rem; font-weight: 700; color: #111827; font-family: 'Space Grotesk', sans-serif;">검색 설정</h4>
</div>
""", unsafe_allow_html=True)

            # --- Identity Cutoff ---
            st.markdown("""
<div style="margin-bottom: 0.4rem;">
    <div style="font-size: 0.95rem; font-weight: 500; color: #1f2937;">기준 유사도 (Identity Cutoff)</div>
    <div style="font-size: 0.8rem; color: #6b7280; margin-top: 0.1rem;">• 필수 사항 (Mandatory)</div>
</div>
""", unsafe_allow_html=True)
            radar_cutoff = st.number_input("Identity Cutoff", value=0.95, min_value=0.0, max_value=1.0, step=0.01, format="%0.2f", label_visibility="collapsed")
            st.markdown("""
<div style="display: flex; align-items: flex-start; gap: 0.3rem; margin-top: -0.6rem; margin-bottom: 1.5rem;">
    <span class="material-symbols-outlined" style="color: #3b82f6; font-size: 1.1rem;">info</span>
    <span style="font-size: 0.85rem; color: #6b7280; padding-top: 0.1rem;">입력 서열과 비교하여 아미노산 서열의 유사도가 기준값 이상인 데이터를 BOARDS DB에서 추출합니다.\n\n- 유사도 범위: 0.0 ~ 1.0 (1.0 = 일치율 100%)</span>
</div>
""", unsafe_allow_html=True)

            # --- Sample Name ---
            st.markdown("""
<div style="margin-bottom: 0.4rem;">
    <div style="font-size: 0.95rem; font-weight: 500; color: #1f2937;">세균의 종 (Sample Species)</div>
    <div style="font-size: 0.8rem; color: #6b7280; margin-top: 0.1rem;">• 선택 사항 (Optional)</div>
</div>
""", unsafe_allow_html=True)
            radar_sample = st.text_input("Sample Name", value="", label_visibility="collapsed")
            st.markdown("""
<div style="display: flex; align-items: flex-start; gap: 0.3rem; margin-top: -0.6rem; margin-bottom: 1.5rem;">
    <span class="material-symbols-outlined" style="color: #3b82f6; font-size: 1.1rem;">info</span>
    <span style="font-size: 0.85rem; color: #6b7280; padding-top: 0.1rem;">분석할 세균의 학명을 이명법으로 입력하세요.\n\n- 예: Enterococcus faecalis</span>
</div>
""", unsafe_allow_html=True)

            # --- Strain Tag ---
            st.markdown("""
<div style="margin-bottom: 0.4rem;">
    <div style="font-size: 0.95rem; font-weight: 500; color: #1f2937;">균주명 (Strain)</div>
    <div style="font-size: 0.8rem; color: #6b7280; margin-top: 0.1rem;">• 선택 사항 (Optional)</div>
</div>
""", unsafe_allow_html=True)
            radar_strain = st.text_input("Strain Tag", value="", label_visibility="collapsed")
            st.markdown("""
<div style="display: flex; align-items: flex-start; gap: 0.3rem; margin-top: -0.6rem; margin-bottom: 0.5rem;">
    <span class="material-symbols-outlined" style="color: #3b82f6; font-size: 1.1rem;">info</span>
    <span style="font-size: 0.85rem; color: #6b7280; padding-top: 0.1rem;">서버 내부 데이터를 구분하는 식별자로 사용되는 균주명을 소문자로 입력하세요.\n\n- 예: Enterococcus faecalis V583에서 균주명인 V583을 v583으로 기입</span>
</div>
""", unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    
    # Bottom centered start button
    start_col1, start_col2, start_col3 = st.columns([1, 1, 1])
    with start_col2:
        start_btn = st.button("🚀 검색 시작 (Start)", type="primary", use_container_width=True)

    if start_btn:
        has_file = bool(uploaded_file)
        has_text = bool(pasted_seq.strip())

        if not has_file and not has_text:
            st.error("서열 파일을 업로드하거나, 서열 텍스트를 입력하세요!")
        else:
            # Automatically use whichever input is actually present
            effective_method = "File_Upload" if has_file else "Text_Input"
            st.session_state["run_id"] = f"RUN_RADAR_{datetime.now().strftime('%Y%m%d%H%M')}"
            st.session_state["job_dir"] = f"results/{st.session_state['run_id']}"
            st.session_state["radar_config"] = {
                "db_name": "BOARDS", "cutoff": radar_cutoff, "input_method": effective_method,
                "sample_name": radar_sample.strip(), "strain_tag": radar_strain.strip()
            }
            
            fasta_path = os.path.join("results", f"{st.session_state['run_id']}.fasta")
            os.makedirs("results", exist_ok=True)
            with open(fasta_path, "w") as f:
                raw_data = uploaded_file.read().decode() if has_file else pasted_seq.strip()
                # Ensure FASTA header exists to prevent IndexError in downstream tools (AlphaFold2, etc.)
                if not raw_data.strip().startswith(">"):
                    header_name = radar_sample.strip() if radar_sample.strip() else "query_sequence"
                    sanitized_header = header_name.replace(" ", "_").replace("\n", "").replace("\r", "")
                    raw_data = f">{sanitized_header}\n{raw_data}"
                f.write(raw_data)
            st.session_state["fasta_path"] = fasta_path
            
            with st.spinner("RADAR 엔진 실행 중..."):
                session = init_db()
                pipe = FinalWorkflowPipeline(session, custom_config={"radar": st.session_state["radar_config"]})
                pipe.save_run(st.session_state["run_id"], "Full_Pipeline", fasta_path)
                try:
                    is_hit, data = pipe.run_actual_radar_engine(fasta_path)
                    if is_hit:
                        st.session_state["pipeline_step"] = "ROUTE_SELECT"
                        st.session_state["radar_data"] = data
                        st.session_state["execution_log"] = pipe.execution_log
                        st.rerun()
                    else:
                        st.error("📡 RADAR 분석 실패")
                except Exception as e:
                    st.error(f"🚨 에러 보사: {e}")

    


# --- [단계 선택] ---
elif st.session_state["pipeline_step"] == "ROUTE_SELECT":
    # JS: force white background on bordered containers in ROUTE_SELECT too
    components.html("""
<script>
(function() {
    function applyWhite() {
        try {
            var doc = window.parent.document;
            var divs = doc.querySelectorAll('div[data-testid="stVerticalBlock"]');
            divs.forEach(function(el) {
                var style = window.parent.getComputedStyle(el);
                if (style.borderTopWidth !== '0px' && style.borderTopWidth !== '') {
                    el.style.setProperty('background-color', '#ffffff', 'important');
                }
            });
        } catch(e) {}
    }
    applyWhite();
    setTimeout(applyWhite, 300);
    setTimeout(applyWhite, 800);
    setTimeout(applyWhite, 2000);
    try {
        var observer = new MutationObserver(function() { applyWhite(); });
        observer.observe(window.parent.document.body, {childList: true, subtree: true});
    } catch(e) {}
})();
</script>
""", height=0)
    radar_res = st.session_state.get("radar_data", {})
    is_protein = radar_res.get("is_protein", False)
    status_msg = "✅ 항생제 내성 유발 단백질 탐지 완료" if is_protein else "✅ 항생제 내성 유전자 탐지 완료"
    
    st.markdown(f"""
<div style="background-color: #ecfdf5; border: 1px solid #a7f3d0; border-radius: 0.5rem; padding: 0.8rem 1rem; display: flex; align-items: center; gap: 0.5rem; color: #047857; font-size: 0.85rem; font-weight: 700; box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05); margin-bottom: 2rem;">
    {status_msg}
</div>
""", unsafe_allow_html=True)

    # --- 서열의 유사성 검색 결과 요약 ---
    st.markdown("""
<div style="display: flex; align-items: center; gap: 0.5rem; margin-bottom: 0.8rem;">
    <span class="material-symbols-outlined" style="color: #2563eb; font-size: 1.5rem;">satellite_alt</span>
    <h2 style="margin: 0; font-size: 1.25rem; font-weight: 700; color: #111827; letter-spacing: -0.02em;">서열의 유사성 검색 결과 요약</h2>
</div>
""", unsafe_allow_html=True)
    
    cnt_table1 = len(radar_res['df_table1']) if radar_res.get("df_table1") is not None else 0
    cnt_table2 = len(radar_res['df_table2']) if radar_res.get("df_table2") is not None else 0
    
    st.markdown(f"""
<div style="background: white; border-radius: 0.75rem; border: 1px solid #e5e7eb; padding: 1.5rem; display: flex; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05); margin-bottom: 2rem;">
    <div style="flex: 1; padding: 0 1.5rem; border-right: 1px solid #e5e7eb;">
        <p style="font-size: 0.75rem; color: #6b7280; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 0.5rem;">항생제 내성 유전자의 단백질 구조 예측 정보 (BOARDS DB) 기록</p>
        <div style="display: flex; align-items: baseline; gap: 0.3rem;">
            <span style="font-size: 2.5rem; font-weight: 700; color: #111827; line-height: 1;">{cnt_table1}</span>
            <span style="font-size: 0.9rem; font-weight: 700; color: #4b5563;">건</span>
        </div>
    </div>
    <div style="flex: 1; padding: 0 1.5rem;">
        <p style="font-size: 0.75rem; color: #6b7280; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 0.5rem;">임상 기록</p>
        <div style="display: flex; align-items: baseline; gap: 0.3rem;">
            <span style="font-size: 2.5rem; font-weight: 700; color: #111827; line-height: 1;">{cnt_table2}</span>
            <span style="font-size: 0.9rem; font-weight: 700; color: #4b5563;">건</span>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)
            
    # --- 다음 단계 분석 경로 선택 ---
    st.markdown("""
<div style="display: flex; align-items: center; gap: 0.5rem; margin-bottom: 1.25rem;">
    <span class="material-symbols-outlined" style="color: #2563eb; font-size: 1.5rem;">construction</span>
    <h2 style="margin: 0; font-size: 1.25rem; font-weight: 700; color: #111827; letter-spacing: -0.02em;">다음 단계 분석 경로 선택</h2>
</div>
""", unsafe_allow_html=True)
    
    col_b, col_af = st.columns([1.05, 2], gap="large")
    with col_b:
        with st.container(border=True):
            st.markdown("""
<style>
/* BOARDS button: targeted using key class (not global) */
.st-key-boards_use_btn button {
    background-color: #2563eb !important;
    border-color: #2563eb !important;
    color: white !important;
    font-weight: 700 !important;
    border-radius: 0.35rem !important;
    height: 3rem !important;
    font-size: 1.05rem !important;
    transition: all 0.2s !important;
    width: 100% !important;
}
.st-key-boards_use_btn button:hover {
    background-color: #1d4ed8 !important;
    border-color: #1d4ed8 !important;
}
.custom-header-blue {
    margin: -1rem -1rem 1rem -1rem;
    background-color: #dbeafe;
    padding: 1.25rem 1.5rem;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    color: #2563eb;
    border-bottom: 1px solid #bfdbfe;
    border-radius: 0.5rem 0.5rem 0 0;
}
</style>
<div class="custom-header-blue">
   <span class="material-symbols-outlined" style="font-size: 1.3rem; font-variation-settings: 'FILL' 1;">database</span>
   <h3 style="margin: 0; font-weight: 700; font-size: 1rem; text-transform: uppercase; letter-spacing: 0.02em;">항생제 내성 유전자의 단백질 구조 예측 정보 (BOARDS DB) 사용</h3>
</div>
<div style="min-height: 290px; padding: 0 0.5rem;">
   <p style="color: #2563eb; font-size: 0.9rem; font-weight: 700; line-height: 1.5; margin-bottom: 1.5rem; word-break: keep-all;">총 3,943개의 항생제 내성 유전자 정보와 27,395개의 예측된 단백질 구조를 포함하는 데이터베이스에 접근하여 빠르게 분석으로 넘어갑니다.</p>
   <div>
      <p style="color: #111827; font-size: 0.85rem; font-weight: 700; margin-bottom: 1rem;">구성 요소:</p>
      <ul style="color: #6b7280; font-size: 0.85rem; font-weight: 500; line-height: 2; padding-left: 1.2rem; margin: 0; list-style-type: disc;">
         <li>[표 1] BOARDS DB 검색 결과 (BOARDS DB Information)</li>
         <li>[표 2] 임상 정보 (Clinical Information)</li>
         <li>[전역적 단백질 구조 신뢰도 분석] 예측 정렬 오차 (PAE) 도표</li>
         <li>[국소적 단백질 구조 신뢰도 분석] pLDDT 점수: 평균 pLDDT 점수, 아미노산 잔기별 신뢰도 도표</li>
         <li>[항생제 내성 유발 단백질의 3차원 구조] 예측 모델 시각화, 파일(.pdb) 다운로드</li>
      </ul>
   </div>
</div>
<hr style="border: none; border-top: 1px solid #e5e7eb; margin: 1.5rem -1rem 1rem -1rem;" />
""", unsafe_allow_html=True)
            if st.button("☁️ BOARDS 데이터 사용", width="stretch", key="boards_use_btn"):
                st.session_state["pipeline_step"] = "BOARDS_HIT"
                st.rerun()

    with col_af:
        with st.container(border=True):
            st.markdown("""
<style>
/* AlphaFold2 start button: targeted using key class (not global) */
.st-key-af2_start_btn button {
    background-color: #ef4444 !important;
    border-color: #ef4444 !important;
    color: white !important;
    font-weight: 700 !important;
    border-radius: 0.35rem !important;
    height: 3.5rem !important;
    font-size: 1.1rem !important;
    transition: all 0.2s !important;
    width: 100% !important;
    margin-top: 0.5rem !important;
}
.st-key-af2_start_btn button:hover {
    background-color: #dc2626 !important;
    border-color: #dc2626 !important;
}

/* AlphaFold Header */
.custom-header-amber {
    margin: -1rem -1rem 1.5rem -1rem;
    background-color: #fef3c7;
    padding: 1.25rem 1.5rem;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    color: #d97706;
    border-bottom: 1px solid #fde68a;
    border-radius: 0.5rem 0.5rem 0 0;
}

/* New inner header pulling to the edges of the inner container */
.custom-inner-settings-header {
    margin: -1rem -1rem 1.5rem -1rem;
    background-color: #f8fafc;
    border-bottom: 1px solid #e5e7eb;
    border-radius: 0.5rem 0.5rem 0 0;
    padding: 1rem 1.25rem;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    color: #111827;
}

/* SLIDERS RED */
[data-baseweb="slider"] div[role="slider"] {
    background-color: #ef4444 !important;
}
[data-baseweb="slider"] div[data-testid="stTickBar"] > div {
    background-color: #ef4444 !important;
}
[data-baseweb="slider"] div[data-testid="stThumbValue"] {
    color: #ef4444 !important;
    font-weight: 700 !important;
}
</style>

<div class="custom-header-amber">
   <span style="font-size: 1.3rem;">🧬</span>
   <h3 style="margin: 0; font-weight: 700; font-size: 1rem; text-transform: uppercase; letter-spacing: 0.02em;">ALPHAFOLD2 실행</h3>
</div>
<div style="padding: 0 0.5rem;">
<p style="color: #d97706; font-size: 1.15rem; font-weight: 700; margin-bottom: 0.5rem;">단백질의 3차원 구조를 생성합니다.</p>
<p style="color: #4b5563; font-size: 0.9rem; font-weight: 500; margin-bottom: 1.5rem;">AlphaFold2에 문제가 있을 시 <a href="https://colab.research.google.com/drive/16V5nFp-LuXAAIA7Hatiez2zzvgpm0xmr?usp=drive_link" style="color: #d97706; font-weight: 700; text-decoration: underline; text-underline-offset: 4px;" target="_blank">이 링크</a>를 이용하세요.</p>
</div>
<style>
/* Alphafold red button forced override */
div[data-testid="column"]:nth-of-type(2) button {
    background-color: #ef4444 !important;
    border-color: #ef4444 !important;
    color: white !important;
    font-weight: 700 !important;
    border-radius: 0.35rem !important;
    height: 3.5rem !important;
    font-size: 1.1rem !important;
    transition: all 0.2s !important;
    width: 100% !important;
    margin-top: 0.5rem !important;
}
div[data-testid="column"]:nth-of-type(2) button:hover {
    background-color: #dc2626 !important;
    border-color: #dc2626 !important;
    opacity: 0.9 !important;
}
/* Selectbox background color modification */
div[data-testid="stSelectbox"] div[data-baseweb="select"] > div {
    background-color: #f1f5f9 !important;
    border-radius: 0.35rem !important;
}
</style>
""", unsafe_allow_html=True)
            
            with st.container(border=True):
                st.markdown("""
<div class="custom-inner-settings-header">
<span class="material-symbols-outlined" style="font-size: 1.15rem;">settings</span>
<span style="font-size: 0.9rem; font-weight: 700; letter-spacing: -0.01em;">AlphaFold2 설정</span>
</div>
""", unsafe_allow_html=True)
                col_s1, col_s2 = st.columns(2)
                with col_s1:
                    rt_m = st.selectbox("Runtime Mode", ["Local (Conda: af_pipeline)", "Docker (Container)"], index=0)
                    msa_m = st.selectbox("MSA Mode", ["mmseqs2_uniref_env", "mmseqs2_uniref", "single_sequence"], index=0)
                    gpu_mem = st.slider("GPU Memory Fraction", 0.1, 1.0, 0.9, help="GPU 메모리 점유율을 조절합니다.")
                
                with col_s2:
                    num_m = st.slider("단백질 구조 모델 수 (Number of Models)", 1, 5, 5)
                    num_r = st.slider("재구축 횟수 (Number of Recycles)", 1, 20, 3)
                    out_name = st.text_input("결과 폴더 이름 (Output Directory Name)", value="af2_results")
                
                docker_img = "ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"
                if "Docker" in rt_m:
                    docker_img = st.text_input("Docker Image URL", value=docker_img)

            st.markdown('<hr style="margin: 1.5rem -1rem 1rem -1rem; border: none; border-top: 1px solid #e5e7eb;" />', unsafe_allow_html=True)
            if st.button("🔥 AlphaFold2 분석 과정 가동", width="stretch", type="primary", key="af2_start_btn"):

                st.session_state["af2_config"] = {
                    "mode": "local" if "Local" in rt_m else "docker",
                    "msa_mode": msa_m, 
                    "num_models": num_m, 
                    "num_recycles": num_r, 
                    "gpu_mem_fraction": gpu_mem,
                    "docker_image": docker_img,
                    "output_dir": out_name
                }
                st.session_state["pipeline_step"] = "AF2_RUNNING"
                st.rerun()
    
# --- [Stage 1] RADAR / BOARDS ---
elif st.session_state["pipeline_step"] == "BOARDS_HIT":
    radar_res = st.session_state["radar_data"]
    execution_log = st.session_state.get("execution_log", [])

    # ── CSS for this page ──────────────────────────────────────────────
    # (Moved to global)

    # ── 3. 서열의 유사성 검색 결과 요약 ────────────────────────────────────────
    st.markdown("""
<div style="display:flex;align-items:center;gap:0.75rem;margin-bottom:1.25rem;">
    <div style="width:2.25rem;height:2.25rem;background:#005AC2;border-radius:0.25rem;display:flex;align-items:center;justify-content:center;box-shadow:0 1px 3px rgba(0,0,0,0.1);">
        <span class="material-symbols-outlined" style="color:white;font-size:1.2rem;">assignment</span>
    </div>
    <h2 style="margin:0;font-size:1.1rem;font-weight:700;color:#1e293b;">서열의 유사성 검색 결과 요약</h2>
</div>""", unsafe_allow_html=True)

    # E-value formatting helper (e.g., 9.97e-237 -> 9.97 × 10⁻²³⁷)
    def format_evalue(e):
        try:
            if pd.isna(e): return ""
            e_flt = float(e)
            if e_flt == 0: return "0"
            s = f"{e_flt:.2e}"
            if "e" not in s: return s
            base, exp = s.split("e")
            exp_int = int(exp)
            superscript_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
            exp_str = str(exp_int).translate(superscript_map)
            return f"{base} × 10{exp_str}"
        except:
            return e

    # 표 1 - BOARDS DB 검색 결과
    df1 = radar_res.get("df_table1")
    if df1 is not None:
        st.markdown('<div class="sci-table-card"><div class="sci-table-header">🧾 [표 1] BOARDS DB 검색 결과 (BOARDS DB Information)</div></div>', unsafe_allow_html=True)
        rename_map1 = {
            "sample_species": "세균의 종", "strain_tag": "균주명", "seq_type": "서열의 종류",
            "q_id": "질의 ID(Query ID)", "temp_id": "BOARDS DB 식별번호",
            "q_start": "시작 위치", "q_end": "종료 위치",
            "snp_marker": "SNP 마커(SNP marker)",
            "p.identity": "단백질 서열 유사도", "bitscore": "Bitscore", "evalue": "E-value"
        }
        df1_view = df1.rename(columns=rename_map1)
        if "E-value" in df1_view.columns:
            df1_view["E-value"] = df1_view["E-value"].apply(format_evalue)
        
        st.dataframe(df1_view, use_container_width=True, hide_index=True)

    st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)

    # 표 2 - 임상 정보
    df2 = radar_res.get("df_table2")
    if df2 is not None:
        st.markdown('<div class="sci-table-card"><div class="sci-table-header">🧾 [표 2] 임상 정보 (Clinical Information)</div></div>', unsafe_allow_html=True)
        rename_map2 = {
            "sample_species": "세균의 종", "strain_tag": "균주명", "seq_type": "서열의 종류",
            "q_id": "질의 ID(Query ID)", "temp_id": "BOARDS DB 식별번호",
            "snp_marker": "SNP 마커(SNP marker)", "amr_gene_class": "내성유전자 분류",
            "amr_gene_detail": "상세 내성유전자 정보",
            "target_antibiotics_class": "표적 항생제 계열",
            "target_antibiotics_major": "주요 표적 항생제",
            "p.identity": "단백질 서열 유사도", "bitscore": "Bitscore", "evalue": "E-value"
        }
        df2_view = df2.rename(columns=rename_map2)
        if "E-value" in df2_view.columns:
            df2_view["E-value"] = df2_view["E-value"].apply(format_evalue)
            
        st.dataframe(df2_view, use_container_width=True, hide_index=True)

    # ── 4. 주요 검색 지표 가이드 ──────────────────────────────────────
    st.markdown("""
<div style="background:rgba(239,246,255,0.5);border:1px solid #bfdbfe;border-radius:0.25rem;padding:1.5rem;margin:2rem 0;box-shadow:0 1px 2px rgba(0,0,0,0.05);">
    <div style="display:flex;align-items:center;gap:0.5rem;margin-bottom:1rem;justify-content:space-between;">
        <div style="display:flex;align-items:center;gap:0.5rem;">
            <span style="font-size:1.2rem;">🧬</span>
            <h3 style="margin:0;font-size:1rem;font-weight:700;color:#1e293b;">주요 검색 지표 가이드</h3>
        </div>
        <a href="https://bio-kcs.tistory.com/entry/BLAST-BLAST-%EC%95%8C%EA%B3%A0%EB%A6%AC%EC%A6%98%EC%97%90-%EB%8C%80%ED%95%B4-%EC%95%8C%EC%95%84%EB%B3%B4%EC%9E%90#5.%20BLAST%20statistical%20significance-1" target="_blank" style="font-size:0.75rem;color:#2563eb;text-decoration:none;display:inline-flex;align-items:center;gap:0.25rem;background:white;padding:0.25rem 0.5rem;border:1px solid #bfdbfe;border-radius:0.25rem;transition:all 0.2s;">
            <span class="material-symbols-outlined" style="font-size:0.875rem;">open_in_new</span>BLAST 평가 지표 이해하기
        </a>
    </div>
    <ul style="margin:0 0 1.25rem 0;padding-left:1.25rem;font-size:0.875rem;color:#334155;line-height:2;">
        <li><span style="font-weight:700;color:#1d4ed8;">단백질 서열 유사도:</span> 질의 서열과 참조 서열 간의 아미노산 서열 일치도, 범위: 0 ~ 1(= 100%)</li>
        <li><span style="font-weight:700;color:#1d4ed8;">Bitscore(비트 점수):</span> 서열 정렬(Sequence alignment) 알고리즘에서 아미노산 서열이 서로 얼마나 유사한지 나타내는 정규화된 점수, Bitscore가 높을수록 두 아미노산 서열이 서로 더 유사하다는 것을 의미하고 E-value는 감소</li>
        <li><span style="font-weight:700;color:#1d4ed8;">E-value:</span> 두 아미노산 서열의 유사성 점수가 특정 값 이상일 확률을 기반으로 <span style="font-weight:700;color:#1d4ed8;">우연히 나타날 수 있는 정렬의 수 (기대값, expected value)</span>, E-value가 작을수록 정렬이 우연이 아닐 가능성이 높아지고, 생물학적 의미가 있을 가능성 증가</li>
    </ul>
    <div style="background:white;border:1px solid #bfdbfe;border-radius:0.25rem;overflow:hidden;">
        <table class="eval-table">
            <thead><tr><th>E-value (기대값)</th><th>Description (설명)</th><th>해석</th></tr></thead>
            <tbody>
                <tr><td class="eval-mono">&gt;&nbsp;10<sup>-1</sup></td><td style="font-style:italic;color:#64748b;">insignificant</td><td><strong style="color:#1e293b;">무의미함:</strong> 우연히 일치했을 확률이 높아 분석 결과로 신뢰하기 어려움.</td></tr>
                <tr><td class="eval-mono">10<sup>-1</sup>&nbsp;~&nbsp;10<sup>-5</sup></td><td style="font-style:italic;color:#64748b;">distantly homologous</td><td><strong style="color:#1e293b;">먼 상동성:</strong> 아주 먼 친척 관계 정도로, 진화적으로 먼 연관성이 있을 수 있음.</td></tr>
                <tr><td class="eval-mono">10<sup>-5</sup>&nbsp;~&nbsp;10<sup>-50</sup></td><td style="font-style:italic;color:#64748b;">homologous proteins</td><td><strong style="color:#1e293b;">상동성 단백질:</strong> 기능이나 구조가 유사한 단백질일 확률이 높음, 유의미한 매칭.</td></tr>
                <tr><td class="eval-mono">10<sup>-50</sup>&nbsp;~&nbsp;10<sup>-100</sup></td><td style="font-style:italic;color:#64748b;">nearly identical</td><td><strong style="color:#1e293b;">거의 일치함:</strong> 서열이 아주 강력하게 일치하며, 같은 단백질로 간주.</td></tr>
                <tr><td class="eval-mono">&lt;&nbsp;10<sup>-100</sup></td><td style="font-style:italic;color:#64748b;">identical proteins</td><td><strong style="color:#1e293b;">완벽히 일치함:</strong> 동일한 단백질로 판단.</td></tr>
            </tbody>
        </table>
    </div>
</div>""", unsafe_allow_html=True)


    # ── 5. 구조 분석 & 시각화 (기존 함수 유지) ────────────────────────
    show_structural_dashboard(radar_res)

    st.markdown("<div style='height:0.2rem'></div>", unsafe_allow_html=True)
    render_model_selection_and_analysis_options(radar_res["models"], is_boards=True)


# --- [AlphaFold2 가동] ---
elif st.session_state["pipeline_step"] == "AF2_RUNNING":
    with st.spinner("AlphaFold2 엔진 구동 중..."):
        from af2_module.af2_runner import run_alphafold_docker, run_alphafold_local
        af_conf = st.session_state["af2_config"]
        af_out_abs = os.path.join(os.getcwd(), st.session_state["job_dir"], af_conf["output_dir"])
        
        try:
            af_out_abs = os.path.join(os.getcwd(), st.session_state["job_dir"], af_conf["output_dir"])
            
            # --- [ Instant Recovery Check ] ---
            # 이미 동일 디렉토리에 결과가 완성되어 있다면 실시간 스킵 제안
            import glob
            existing_pdb = glob.glob(os.path.join(af_out_abs, "*.pdb"))
            if existing_pdb:
                st.warning("⚠️ 해당 폴더에 연산 결과가 이미 존재합니다.")
                if st.button("🚀 연산 스킵 및 결과 바로 보기", type="primary"):
                    st.session_state["af2_models"] = parse_existing_af2_results(af_out_abs)
                    st.session_state["pipeline_step"] = "AF2_RESULTS"
                    st.rerun()
                st.markdown("---")
            
            af_start = datetime.now()
            
            # 결정: RADAR 단계에서 Prodigal에 의해 번역된 번역본(query_faa)이 존재한다면 사용하고, 없으면 사용자가 입력한 원본을 사용
            af2_input_seq = st.session_state.get("radar_data", {}).get("query_faa")
            if not af2_input_seq or not os.path.exists(af2_input_seq):
                af2_input_seq = st.session_state["fasta_path"]
                
            if af_conf["mode"] == "docker":
                gen = run_alphafold_docker(
                    af2_input_seq, 
                    af_out_abs, 
                    msa_mode=af_conf["msa_mode"],
                    docker_image=af_conf["docker_image"],
                    gpu_mem_fraction=str(af_conf["gpu_mem_fraction"]),
                    num_models=af_conf["num_models"],
                    num_recycles=af_conf["num_recycles"]
                )
            else:
                gen = run_alphafold_local(
                    af2_input_seq, 
                    af_out_abs, 
                    msa_mode=af_conf["msa_mode"], 
                    num_models=af_conf["num_models"], 
                    num_recycles=af_conf["num_recycles"],
                    gpu_mem_fraction=str(af_conf["gpu_mem_fraction"])
                )
            
            log_placeholder = st.empty()
            full_log = ""
            for line in gen:
                full_log += line
                display_log = "\n".join(full_log.splitlines()[-50:])
                log_placeholder.code(display_log, language="bash")
            
            af_end = datetime.now()
            
            # [점검] 결과 파싱 로직 강화 (ColabFold 1.6.1+ 및 다양한 파일명 대응)
            import glob
            import re
            
            # PDB 파일 검색 및 정렬
            pdb_files = sorted(glob.glob(os.path.join(af_out_abs, "*.pdb")))
            if not pdb_files: 
                raise RuntimeError(f"PDB 파일을 찾을 수 없습니다. (경로: {af_out_abs})")
            
            # rank_(\d+).pdb 또는 ranked_(\d+).pdb 또는 unrelaxed_..._model_(\d+).pdb 패턴 대응
            def extract_rank_from_pdb(f):
                # 1. ranked_(\d+) or rank_(\d+)
                m = re.search(r'rank_?(\d+)', f)
                if m: return int(m.group(1))
                # 2. _model_(\d+)
                m = re.search(r'model_(\d+)', f)
                if m: return int(m.group(1))
                return 999
            
            pdb_files.sort(key=extract_rank_from_pdb)
            
            pipe = FinalWorkflowPipeline(init_db())
            # 모든 JSON 파일 (scores, predicted_aligned_error 등 모두 포함)
            all_json_files = glob.glob(os.path.join(af_out_abs, "*.json"))
            models_data = []
            
            for i, pdb_f in enumerate(pdb_files):
                # 파일명에서 랭크/모델 번호 추출
                rank_id = str(extract_rank_from_pdb(pdb_f))
                
                plddt_list = pipe.get_plddt_stats_fixed(pdb_f)
                avg_plddt = np.mean(plddt_list) if plddt_list else 0.0
                with open(pdb_f, "r") as f: pdb_content = f.read()
                
                # PAE 데이터 추출 시도
                pae_data = []
                # 1. 랭크/모델 번호가 일치하는 JSON 파일 탐색 (scores 또는 predicted_aligned_error 파일)
                potential_jsons = [j for j in all_json_files if f"model_{rank_id}" in j or f"rank_{rank_id}" in j]
                
                # 2. 만약 해당 패턴이 없으면 i(0-based) 순으로 매칭 시도
                if not potential_jsons:
                    potential_jsons = [j for j in all_json_files if f"_{i}." in j or f"_{i}_" in j]

                if potential_jsons:
                    try:
                        # 가장 관련성 높은 JSON (scores 또는 pae 가공 완료된 파일) 우선 시도
                        with open(potential_jsons[0], 'r') as jf:
                            js = json.load(jf)
                            # ColabFold/AlphaFold 다양한 키 대응 ('pae', 'predicted_aligned_error')
                            if isinstance(js, list): js = js[0]
                            pae_data = js.get('pae') or js.get('predicted_aligned_error') or []
                    except: pass

                models_data.append({
                    "rank": i,
                    "display_rank": int(rank_id) if rank_id != "999" else i,
                    "file_name": os.path.basename(pdb_f), 
                    "pdb_path": pdb_f,
                    "plddt": avg_plddt,
                    "plddt_list": plddt_list, 
                    "pae_data": pae_data, 
                    "pdb_data": pdb_content
                })
            
            # [ADD] 공식 PAE 그리드 이미지 탐색
            pae_image_path = None
            png_files = glob.glob(os.path.join(af_out_abs, "*_pae.png"))
            if png_files:
                pae_image_path = png_files[0]

            st.session_state["af2_models"] = {
                "models": models_data,
                "pae_image": pae_image_path
            }
            st.session_state["execution_log"] = st.session_state.get("execution_log", []) + [{
                "task": "AlphaFold2", "start": af_start.strftime("%H:%M:%S"), "end": af_end.strftime("%H:%M:%S"),
                "elapsed": str(af_end - af_start).split(".")[0], "status": "Success"
            }]
            st.session_state["pipeline_step"] = "AF2_RESULTS"
            st.rerun()
            
        except Exception as e:
            st.error(f"❌ AlphaFold2 에러: {e}")

# --- [AlphaFold2 결과] ---
elif st.session_state["pipeline_step"] == "AF2_RESULTS":
    # --- [데이터 자동 복구 로직] ---
    # 만약 결과창인데 모델 데이터가 없거나 PAE가 비어있다면 디스크에서 재파싱 시도
    af2_data = st.session_state.get("af2_models")
    needs_recovery = False
    if not af2_data:
        needs_recovery = True
    else:
        # Check for models list in both possible formats (list or dict)
        models_list = af2_data.get("models", []) if isinstance(af2_data, dict) else af2_data
        if not models_list:
            needs_recovery = True
        elif not models_list[0].get('pae_data') and not (isinstance(af2_data, dict) and af2_data.get('pae_image')):
            # Only recover if BOTH raw PAE data and official PNG are missing
            needs_recovery = True

    if needs_recovery:
        job_dir = st.session_state.get("job_dir")
        af_out = st.session_state.get("af2_config", {}).get("output_dir", "af2_results")
        if job_dir:
            af_out_abs = os.path.join(os.getcwd(), job_dir, af_out)
            if os.path.exists(af_out_abs):
                st.info("📡 결과 데이터를 디스크에서 복구하는 중입니다...")
                # 재파싱 로직 수행 (AF2_RUNNING의 파싱 로직과 동일)
                # ... (아래에서 구현)
                import glob
                import re
                pdb_files = sorted(glob.glob(os.path.join(af_out_abs, "*.pdb")))
                if pdb_files:
                    def extract_rank_from_pdb(f):
                        m = re.search(r'rank_?(\d+)', f)
                        if m: return int(m.group(1))
                        m = re.search(r'model_(\d+)', f)
                        if m: return int(m.group(1))
                        return 999
                    pdb_files.sort(key=extract_rank_from_pdb)
                    
                    recover_models = []
                    all_json_files = glob.glob(os.path.join(af_out_abs, "*.json"))
                    pipe = FinalWorkflowPipeline(init_db())
                    
                    for i, pdb_f in enumerate(pdb_files):
                        rank_id = str(extract_rank_from_pdb(pdb_f))
                        plddt_list = pipe.get_plddt_stats_fixed(pdb_f)
                        avg_plddt = np.mean(plddt_list) if plddt_list else 0.0
                        with open(pdb_f, "r") as f: pdb_content = f.read()
                        
                        pae_data = []
                        potential_jsons = [j for j in all_json_files if f"model_{rank_id}" in j or f"rank_{rank_id}" in j]
                        if not potential_jsons:
                             potential_jsons = [j for j in all_json_files if f"_{i}." in j or f"_{i}_" in j]
                        if potential_jsons:
                            try:
                                with open(potential_jsons[0], 'r') as jf:
                                    js = json.load(jf)
                                    if isinstance(js, list): js = js[0]
                                    pae_data = js.get('pae') or js.get('predicted_aligned_error') or []
                            except: pass
                        
                        recover_models.append({
                            "rank": i, "display_rank": int(rank_id) if rank_id != "999" else i,
                            "file_name": os.path.basename(pdb_f), "pdb_path": pdb_f,
                            "plddt": avg_plddt, "plddt_list": plddt_list, "pae_data": pae_data, "pdb_data": pdb_content
                        })
                if pdb_files:
                    # ... (생략된 랭킹 추출 로직)
                    # ...
                    pae_image_path = None
                    png_files = glob.glob(os.path.join(af_out_abs, "*_pae.png"))
                    if png_files:
                        pae_image_path = png_files[0]
                    
                    st.session_state["af2_models"] = {
                        "models": recover_models,
                        "pae_image": pae_image_path
                    }
                    st.success("✅ 결과 복구 완료!")
                    st.rerun()
    
    with st.expander("📂 데이터 파일 경로 확인", expanded=False):
        af2_data = st.session_state["af2_models"]
        models = af2_data.get("models", []) if isinstance(af2_data, dict) else af2_data
        st.write(f"**결과 디렉토리:** `{os.path.dirname(models[0]['pdb_path'])}`")
        for i, m in enumerate(models[:3]):
            st.code(f"Model {i+1} PDB: {m['pdb_path']}\nAvg pLDDT: {m['plddt']:.2f}")

    show_structural_dashboard(st.session_state["af2_models"])
    st.markdown("<hr style='margin:0.5rem 0; border-top:1px solid #f1f5f9; opacity:0.3;'>", unsafe_allow_html=True)
    af2_data = st.session_state["af2_models"]
    models_to_pass = af2_data.get("models") if isinstance(af2_data, dict) else af2_data
    render_model_selection_and_analysis_options(models_to_pass, is_boards=False)
