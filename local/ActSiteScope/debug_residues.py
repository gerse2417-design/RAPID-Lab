import re

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

def format_residue_info_columns_FIXED(residue_str, res_map):
    if not residue_str: return "", ""
    
    # 1. 쉼표와 공백을 모두 처리하기 위해 정규식으로 분리
    # 먼저 쉼표 주변 공백 정리 후 쉼표로 분리하거나, 그냥 공백/쉼표 복합 분리
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
                # 3. 최후의 보루: 원본 이름이 있으면 사용, 아니면 UNK
                r_name = orig_name if orig_name else "UNK"
                target_id = f"{ch_in}{n}"
        
        ids.append(target_id)
        names.append(f"{r_name} {n}")
        
    return ", ".join(ids), ", ".join(names)

res_map = {('A', '70'): 'TYR', ('A', '104'): 'THR'}

print("--- Testing MPBind format (THR 70, UNK 104) ---")
print(format_residue_info_columns_FIXED("THR 70, UNK 104", res_map))

print("--- Testing P2Rank format (A_70 A_104) ---")
print(format_residue_info_columns_FIXED("A_70 A_104", res_map))

print("--- Testing fPocket format (SER 70, LYS 73) ---")
print(format_residue_info_columns_FIXED("SER 70, LYS 73", res_map))

