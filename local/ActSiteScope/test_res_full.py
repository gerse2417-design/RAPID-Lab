import sys
import numpy as np

# Mock implementation to bypass streamlit issues
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

def format_residue_info_columns(residue_str, res_map):
    if not residue_str: return "", ""
    # 쉼표로 분리하여 각 잔기를 개별 처리
    items = [i.strip() for i in str(residue_str).split(',')]
    ids = []
    names = []
    
    for item in items:
        # 공백으로 나누어 마지막 부분을 번호로 간주 (예: "THR 70" -> "70", "70A" -> "70A")
        parts = item.split()
        if not parts: continue
        n = parts[-1]
        
        # 1. 체인 'A'에서 우선 조회
        r_name = res_map.get(('A', n))
        target_id = f"A{n}"
        
        # 2. 실패 시 전체 체인을 검색하여 번호가 매칭되는 아미노산명 및 체인ID 확보
        if r_name is None:
            found = False
            for (ch, seq), name in res_map.items():
                if seq == n:
                    r_name = name
                    target_id = f"{ch}{n}"
                    found = True
                    break
            if not found:
                r_name = "UNK"
                target_id = f"A{n}"
        
        ids.append(target_id)
        names.append(f"{r_name} {n}")
        
    return ", ".join(ids), ", ".join(names)

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
            html += f'<th{w_attr_res}>잔기 이름</th>'
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

with open("results/tmp_input/1m40.pdb1", "r") as f:
    pdb_str = f.read()

res_map = get_pdb_residue_map(pdb_str)
hotspots = [{"site": "Site 1", "prob": 0.9, "residues": "UNK 26, UNK 71, UNK 98"}]
html = generate_styled_table(hotspots, ["site","prob","residues"], res_map, widths=[15,20,35,30])
print(html)
