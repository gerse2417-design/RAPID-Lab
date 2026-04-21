import numpy as np
from Bio import pairwise2


########## Process PDB file ##########
def get_pdb_xyz(pdb_file):
    current_pos = -1000
    X = []
    current_aa = {} # N, CA, C, O, R

    def save_residue(curr_aa, X_list):
        if curr_aa != {}:
            atom_mm = ["N", "CA", "C", "O"]
            cur_mean = None
            for aa in atom_mm:
                if aa in curr_aa:
                    cur_mean = curr_aa[aa]
                    break
            if not isinstance(cur_mean, np.ndarray):
                for aa in curr_aa:
                    cur_mean = curr_aa[aa]
                    break
            for aa in atom_mm:
                if aa not in curr_aa:
                    curr_aa[aa] = cur_mean
            
            R_group = []
            for atom in curr_aa:
                if atom not in ["N", "CA", "C", "O"]:
                    R_group.append(curr_aa[atom])
            if R_group == []:
                R_group = [curr_aa["CA"]]
            R_group = np.array(R_group).mean(0)
            X_list.append([curr_aa["N"], curr_aa["CA"], curr_aa["C"], curr_aa["O"], R_group])

    for line in pdb_file:
        line_atom = line[0:4].strip()
        if line_atom == "ATOM":
            try:
                res_id = int(line[22:26].strip())
            except:
                res_id = current_pos
            
            if res_id != current_pos:
                save_residue(current_aa, X)
                current_aa = {}
                current_pos = res_id
            
            atom = line[12:16].strip()
            if atom != "H":
                try:
                    xyz = np.array([line[30:38].strip(), line[38:46].strip(), line[46:54].strip()]).astype(np.float32)
                    current_aa[atom] = xyz
                except: pass
        elif line_atom == "TER":
            save_residue(current_aa, X)
            current_aa = {}
            current_pos = -1000
    
    # Save the last residue
    save_residue(current_aa, X)
    
    if len(X) == 0:
        return np.zeros((1, 5, 3), dtype=np.float32)
        
    return np.array(X)


########## Get DSSP ##########
def process_dssp(dssp_file):
    aa_type = "ACDEFGHIKLMNPQRSTVWY"
    SS_type = "HBEGITSC"
    rASA_std = [115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
                185, 160, 145, 180, 225, 115, 140, 155, 255, 230]

    try:
        with open(dssp_file, "r") as f:
            lines = f.readlines()
    except:
        return "", []

    seq = ""
    dssp_feature = []

    p = 0
    while p < len(lines) and lines[p].strip().find("#") == -1:
        p += 1
    
    if p >= len(lines):
        return "", []

    for i in range(p + 1, len(lines)):
        if len(lines[i]) < 14: continue
        aa = lines[i][13]
        if aa == "!" or aa == "*":
            continue
        seq += aa
        SS = lines[i][16] if len(lines[i]) > 16 else " "
        if SS == " ":
            SS = "C"
        SS_vec = np.zeros(8)
        try:
            idx = SS_type.find(SS)
            if idx != -1: SS_vec[idx] = 1
        except: pass

        try:
            ASA = float(lines[i][34:38].strip())
            RSA = min(1, ASA / rASA_std[aa_type.find(aa)]) # relative solvent accessibility
        except:
            RSA = 0.0
            
        dssp_feature.append(np.concatenate((np.array([RSA]), SS_vec)))

    return seq, dssp_feature


def match_dssp(seq, dssp, ref_seq):
    if not seq or not dssp:
        return [np.zeros(9) for _ in range(len(ref_seq))]
        
    alignments = pairwise2.align.globalxx(ref_seq, seq)
    if not alignments:
         return [np.zeros(9) for _ in range(len(ref_seq))]
         
    target_ref_seq = alignments[0].seqA
    target_seq = alignments[0].seqB

    padded_item = np.zeros(9)

    new_dssp = []
    dssp_ptr = 0
    for aa in target_seq:
        if aa == "-":
            new_dssp.append(padded_item)
        else:
            if dssp_ptr < len(dssp):
                new_dssp.append(dssp[dssp_ptr])
                dssp_ptr += 1
            else:
                new_dssp.append(padded_item)

    matched_dssp = []
    for i in range(len(target_ref_seq)):
        if target_ref_seq[i] == "-":
            continue
        matched_dssp.append(new_dssp[i])

    # Ensure it matches ref_seq length
    if len(matched_dssp) < len(ref_seq):
        matched_dssp += [padded_item] * (len(ref_seq) - len(matched_dssp))
    elif len(matched_dssp) > len(ref_seq):
        matched_dssp = matched_dssp[:len(ref_seq)]

    return matched_dssp
