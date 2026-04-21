import os
import re
import pandas as pd
import numpy as np

class ConsensusEngine:
    """
    Consensus Analysis Engine for Antibiotic Inhibitor Hotspot Selection.
    Implements advanced CASE A/B/C logic with multi-level tie-breaking.
    """
    def __init__(self, pdb_path, mp_hotspots, p2r_results, fpc_results, mcsa_res, **kwargs):
        self.pdb_path = pdb_path
        self.mp_hotspots = mp_hotspots
        self.p2r_results = p2r_results
        self.fpc_results = fpc_results
        self.mcsa_res = mcsa_res
        
        # Keyword arguments with safe defaults (handles old/new names)
        # Try both 'p2r_res_scores_dict' and 'p2r_res_scores' for maximum compatibility
        self.p2r_res_scores = kwargs.get('p2r_res_scores_dict', kwargs.get('p2r_res_scores'))
        if self.p2r_res_scores is None: self.p2r_res_scores = {}
        
        self.mp_th = kwargs.get('mp_th', 0.5)
        p2r_pct = kwargs.get('p2r_pct', 30)
        self.p2r_th = float(1.0 - (p2r_pct / 100.0))
        self.fpc_th = kwargs.get('fpc_th', 0.5)
        self.apbs_pct = float(kwargs.get('apbs_pct', 30) / 100.0)
        
        self.residue_map = {} 
        self.pocket_residues = set()
        self.residue_centers = {} # (chain, res_num) -> [x, y, z]

    def _get_nums(self, res_str):
        if not res_str: return set()
        return set(re.findall(r'\d+', str(res_str)))

    def _update_map(self, res_nums, tool_name, increment=1, chain='A', raw_val=0):
        for num in res_nums:
            try:
                key = (chain, int(num))
                if key not in self.residue_map:
                    self.residue_map[key] = {
                        'count': 0, 'tools': set(), 'res_name': 'UNK', 
                        'is_mcsa': False, 'raw_apbs': 0, 'p2r_res_score': 0
                    }
                
                if tool_name not in self.residue_map[key]['tools']:
                    self.residue_map[key]['tools'].add(tool_name)
                    self.residue_map[key]['count'] += increment
                    
                if tool_name == "APBS":
                    self.residue_map[key]['raw_apbs'] = raw_val
                if tool_name == "P2Rank":
                    self.residue_map[key]['p2r_res_score'] = self.p2r_res_scores.get(key, 0)
            except: continue

    def _get_distance(self, p1, p2):
        if p1 is None or p2 is None: return 999.0
        return np.linalg.norm(np.array(p1) - np.array(p2))

    def run(self):
        # 0. Prep: Get residue centers and names from PDB manually (No Biopython)
        if os.path.exists(self.pdb_path):
            with open(self.pdb_path, 'r') as f:
                res_coords = {} # (chain, res_num) -> list of coords
                res_names = {}  # (chain, res_num) -> name
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        try:
                            ch = line[21].strip() or 'A'
                            num = int(line[22:26].strip())
                            name = line[17:20].strip()
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            key = (ch, num)
                            if key not in res_coords: 
                                res_coords[key] = []
                                res_names[key] = name
                            res_coords[key].append([x, y, z])
                        except: continue
                
                for key, coords in res_coords.items():
                    self.residue_centers[key] = np.mean(coords, axis=0)
                    # Initialize names in residue_map
                    if key not in self.residue_map:
                        self.residue_map[key] = {
                            'count': 0, 'tools': set(), 'res_name': res_names[key], 
                            'is_mcsa': False, 'raw_apbs': 0, 'p2r_res_score': 0
                        }
                    else:
                        self.residue_map[key]['res_name'] = res_names[key]

        # 1. Primary Filtering
        self.mp_hotspots = [h for h in self.mp_hotspots if h.get("prob", 0) >= self.mp_th]
        if self.p2r_results:
            max_p2r = max(p.get("score", 0) for p in self.p2r_results)
            self.p2r_results = [p for p in self.p2r_results if p.get("score", 0) >= max_p2r * self.p2r_th]
        self.fpc_results = [f for f in self.fpc_results if f.get("score", 0) >= self.fpc_th]

        # 2. Define Pocket Residues
        for p in self.p2r_results: self.pocket_residues.update(self._get_nums(p.get("residues", "")))
        for f in self.fpc_results: self.pocket_residues.update(self._get_nums(f.get("residues", "")))

        # 3. Aggregate results
        mp_all = set()
        for h in self.mp_hotspots: mp_all.update(self._get_nums(h.get("residues", "")))
        self._update_map(mp_all, "MPBind", increment=1.0)
        
        p2r_all = set()
        for p in self.p2r_results: p2r_all.update(self._get_nums(p.get("residues", "")))
        self._update_map(p2r_all, "P2Rank", increment=1.5)
        
        fpc_all = set()
        for f in self.fpc_results: fpc_all.update(self._get_nums(f.get("residues", "")))
        self._update_map(fpc_all, "fPocket", increment=1.0)

        # 4. Integrate APBS
        res_ep_map = {} 
        for h_list in [self.mp_hotspots, self.p2r_results, self.fpc_results]:
            for item in h_list:
                raw_ep = item.get("raw_apbs_score", 0)
                norm_ep = item.get("apbs_score", 0)
                nums = self._get_nums(item.get("residues", ""))
                for n in nums:
                    if n not in res_ep_map or abs(raw_ep) > abs(res_ep_map[n][0]):
                        res_ep_map[n] = (raw_ep, norm_ep)

        if res_ep_map:
            all_abs_eps = sorted([abs(v[0]) for v in res_ep_map.values()], reverse=True)
            threshold_idx = max(0, int(len(all_abs_eps) * self.apbs_pct) - 1)
            abs_threshold = all_abs_eps[threshold_idx] if all_abs_eps else 0
            for res_num, (raw_ep, norm_ep) in res_ep_map.items():
                if abs(raw_ep) >= abs_threshold and raw_ep != 0:
                    if norm_ep >= 0.5:
                        inc = 2.0 if norm_ep >= 0.9 else 1.0
                        if res_num in self.pocket_residues:
                            self._update_map([res_num], "APBS", increment=inc, raw_val=raw_ep)

        # 5. Case Logic Implementation
        final_selected = []
        mcsa_keys = []
        for r in self.mcsa_res:
            try: 
                key = ('A', int(r['res_num']))
                mcsa_keys.append(key)
                
                # M-CSA 잔기의 EP 수치 병합 (st.session_state에서 넘어온 값 활용)
                if key in self.residue_map:
                    if 'raw_apbs_score' in r and r['raw_apbs_score'] != 0:
                        self.residue_map[key]['raw_apbs'] = r['raw_apbs_score']
                    elif 'ep_value' in r and r['ep_value'] != 0:
                        self.residue_map[key]['raw_apbs'] = r['ep_value']
            except: continue

        # Helper for sorting
        def get_tiebreak_data(key):
            info = self.residue_map.get(key, {'count': 0, 'raw_apbs': 0, 'p2r_res_score': 0})
            res_coord = self.residue_centers.get(key)
            
            # Dist to M-CSA Anchor
            min_mcsa_dist = 999.0
            for mk in mcsa_keys:
                d = self._get_distance(res_coord, self.residue_centers.get(mk))
                if d < min_mcsa_dist: min_mcsa_dist = d
            
            # Dist to Nearest Pocket Center
            min_pocket_dist = 999.0
            for p in self.p2r_results + self.fpc_results:
                d = self._get_distance(res_coord, p.get("center"))
                if d < min_pocket_dist: min_pocket_dist = d
                
            return {
                "score": info['count'],
                "norm_score": min(info['count'] / 5.5, 1.0),
                "mcsa_dist": min_mcsa_dist,
                "pocket_dist": min_pocket_dist,
                "apbs_ep": abs(info['raw_apbs']),
                "p2r_res": info['p2r_res_score'],
                "index": key[1]
            }

        if mcsa_keys:
            # --- CASE A: M-CSA Based ---
            # 1. Include all Anchors
            anchors = [self._build_res_info(k, get_tiebreak_data(k)) for k in mcsa_keys]
            if len(anchors) > 8:
                final_selected = anchors # User: "초과할 경우 앵커만 반영하고 마무리"
            else:
                # 2. Add Neighbors (6A, pocket, score >= 2.0)
                neighbors = []
                seen_keys = set(mcsa_keys)
                for key, center in self.residue_centers.items():
                    if key in seen_keys: continue
                    tb = get_tiebreak_data(key)
                    if tb['mcsa_dist'] <= 6.0 and str(key[1]) in self.pocket_residues and tb['score'] >= 2.0:
                        neighbors.append((key, tb))
                
                # Sort neighbors: Score(D), MCSA_dist(A), Pocket_dist(A), APBS(D), Index(A)
                neighbors.sort(key=lambda x: (-x[1]['score'], x[1]['mcsa_dist'], x[1]['pocket_dist'], -x[1]['apbs_ep'], x[1]['index']))
                
                final_selected = anchors + [self._build_res_info(n[0], n[1]) for n in neighbors[:8-len(anchors)]]
        else:
            # --- CASE B: De novo ---
            candidates = []
            for key, center in self.residue_centers.items():
                tb = get_tiebreak_data(key)
                if str(key[1]) in self.pocket_residues and tb['score'] >= 2.0:
                    candidates.append((key, tb))
            
            # Sort: Score(D), Pocket_dist(A), APBS(D), P2R_ResScore(D), Index(A)
            candidates.sort(key=lambda x: (-x[1]['score'], x[1]['pocket_dist'], -x[1]['apbs_ep'], -x[1]['p2r_res'], x[1]['index']))
            final_selected = [self._build_res_info(c[0], c[1]) for c in candidates[:8]]

        # --- CASE C: Fallback ---
        if len(final_selected) < 3 and self.p2r_results:
            top_p = self.p2r_results[0]
            p1_res_nums = self._get_nums(top_p.get("residues", ""))
            fb_candidates = []
            for n in p1_res_nums:
                k = ('A', int(n))
                tb = get_tiebreak_data(k)
                fb_candidates.append((k, tb))
            
            # Sort: Pocket_dist(A), P2R_ResScore(D), Score(D), APBS(D)
            fb_candidates.sort(key=lambda x: (x[1]['pocket_dist'], -x[1]['p2r_res'], -x[1]['score'], -x[1]['apbs_ep']))
            final_selected = [self._build_res_info(c[0], c[1]) for c in fb_candidates[:6]]

        return final_selected

    def _build_res_info(self, key, tb):
        info = self.residue_map.get(key, {'tools': set(), 'res_name': 'UNK', 'is_mcsa': False})
        
        is_mcsa = key in [('A', int(r['res_num'])) for r in self.mcsa_res]
        tools_list = list(info.get('tools', []))
        if is_mcsa and "M-CSA" not in tools_list:
            tools_list.append("M-CSA")
            
        res_info = {
            "Chain": key[0],
            "ResNo": key[1],
            "ResName": info.get('res_name', 'UNK'),
            "M-CSA": "Yes" if is_mcsa else "No",
            "Count_pred": tb['score'],
            "Tools": ", ".join(sorted(tools_list)),
            "Score": round(tb['norm_score'], 1) if not is_mcsa else 1.0,
            "ep_value": info.get('raw_apbs', 0)
        }
        
        center = self.residue_centers.get(key)
        if center is not None:
            res_info["center"] = center.tolist() if hasattr(center, "tolist") else center
            
        return res_info

    def generate_report(self, hotspots, output_path, **kwargs):
        rfd_str = ",".join([f"{h['Chain']}{h['ResNo']}" for h in hotspots])
        res_str = ", ".join([f"{h['ResName']} {h['ResNo']}" for h in hotspots])
        
        with open(output_path, "w") as f:
            f.write("ActSiteScope: Final Inhibitor Hotspot Report\n")
            f.write(f"RFdiffusion Format: [{rfd_str}]\n")
            f.write(f"Researcher Format: {res_str}\n\n")
            df = pd.DataFrame(hotspots)
            if not df.empty: f.write(df[["Chain", "ResNo", "ResName", "M-CSA", "Tools", "Score"]].to_string(index=False))
        return rfd_str, res_str
