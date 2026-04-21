import logging

def calculate_motif_concordance(predicted_residues, mcsa_anchor_residues):
    """
    Calculate similarity between predicted pocket residues and M-CSA 'Answer Key'.
    Returns a score from 0.0 to 1.0.
    """
    if not mcsa_anchor_residues:
        return 0.0
    
    # Standardize residues for comparison (e.g. "SER 70")
    pred_set = set([r.strip().upper() for r in predicted_residues.split(",")])
    
    # M-CSA residues list of dicts: [{'res_name': 'SER', 'res_num': 70}, ...]
    anchor_set = set([f"{r['res_name'].upper()} {r['res_num']}" for r in mcsa_anchor_residues])
    
    intersection = pred_set.intersection(anchor_set)
    if not anchor_set:
        return 0.0
        
    # Recall of anchor residues is a good measure for "concordance"
    concordance = len(intersection) / len(anchor_set)
    return concordance

def rank_pockets(pockets, mcsa_residues, electro_scores=None):
    """
    Rank pockets based on Motif Concordance, P2Rank score, and APBS potential.
    Adds 'tier', 'bio_score', and 'concordance' to each pocket.
    """
    ranked_pockets = []
    
    for i, p in enumerate(pockets):
        concordance = calculate_motif_concordance(p["residues"], mcsa_residues)
        
        # Decide Tier: 1 if it contains any catalytic residue
        tier = 1 if concordance > 0 else 2
        
        # Composite Bio-Score calculation
        # BioScore = w1 * Concordance + w2 * (Prob or Score) + w3 * APBS_Score
        w1, w2, w3 = 0.5, 0.3, 0.2
        apbs_s = electro_scores[i] if electro_scores and i < len(electro_scores) else 0.5
        
        # Determine base score from available keys
        base_score = p.get("prob", p.get("score", 0.5))
        bio_score = (w1 * concordance) + (w2 * base_score) + (w3 * apbs_s)
        
        p_updated = p.copy()
        p_updated.update({
            "tier": tier,
            "bio_score": bio_score,
            "apbs_score": apbs_s,
            "concordance": concordance,
            "status": "🏆 Tier 1 (Motif Anchored)" if tier == 1 else "🔬 Tier 2 (Pure Geometric)"
        })
        ranked_pockets.append(p_updated)
        
    # Final sorting: Tier 1 first, then by Bio-Score
    return sorted(ranked_pockets, key=lambda x: (x["tier"], -x["bio_score"]))

if __name__ == "__main__":
    test_pockets = [
        {"id": 1, "residues": "SER 70, LYS 73, SER 130", "score": 0.85},
        {"id": 2, "residues": "ASP 101, GLY 102", "score": 0.72}
    ]
    test_mcsa = [{"res_name": "SER", "res_num": 70}, {"res_name": "LYS", "res_num": 73}]
    
    ranked = rank_pockets(test_pockets, test_mcsa)
    for p in ranked:
        print(f"Site {p['id']} | Tier: {p['tier']} | Score: {p['bio_score']:.2f}")
