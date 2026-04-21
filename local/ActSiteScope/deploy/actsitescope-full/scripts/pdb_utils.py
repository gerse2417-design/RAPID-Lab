import re
import os

def extract_pdb_info(pdb_content, pdb_filename=None):
    """
    Extracts PDB metadata from file content.
    Returns: {"pdb_id": str or None, "title_keyword": str or None, "uniprot_id": str or None}
    """
    pdb_id = None
    title_lines = []
    uniprot_id = None
    
    lines = pdb_content.splitlines()
    for line in lines:
        # PDB ID extraction from HEADER
        if line.startswith("HEADER"):
            id_field = line[62:66].strip()
            if len(id_field) == 4 and id_field.isalnum():
                pdb_id = id_field
            else:
                parts = line.split()
                if len(parts) >= 2:
                    potential_id = parts[-1]
                    if len(potential_id) == 4 and potential_id.isalnum():
                        pdb_id = potential_id
        
        # UniProt ID extraction from DBREF
        if line.startswith("DBREF"):
            # Format usually: DBREF  1M40 A   26   290  UNP    P62593   BLAT_ECOLX      26    290
            # Some files might have different spacing.
            parts = line.split()
            if "UNP" in parts or "unp" in (p.lower() for p in parts):
                # The UniProt ID is usually the first alnum string after UNP
                try:
                    idx = -1
                    if "UNP" in parts: idx = parts.index("UNP")
                    else: idx = [p.upper() for p in parts].index("UNP")
                    
                    if idx != -1 and idx + 1 < len(parts):
                        potential_unp = parts[idx + 1]
                        if potential_unp and len(potential_unp) >= 6:
                            uniprot_id = potential_unp
                except: pass
        
        # Protein Name extraction from TITLE
        if line.startswith("TITLE"):
            content = line[10:].strip()
            if len(content) > 0 and content[0].isdigit():
                title_lines.append(content[1:].strip())
            else:
                title_lines.append(content)

    full_title = " ".join(title_lines).strip()
    title_keyword = None
    if full_title:
        # Use the most significant word (usually last or second to last in beta-lactamase titles)
        words = [w.strip(",.;") for w in full_title.split() if len(w) > 3]
        if words:
            title_keyword = words[-1]
    
    # Fallback for PDB ID if header has 'XXXX' or is missing
    if (not pdb_id or pdb_id.upper() == "XXXX") and pdb_filename:
        # Look for 4-char PDB ID in filename (e.g. 1M40)
        matches = re.findall(r'\b([0-9][a-zA-Z0-9]{3})\b', os.path.basename(pdb_filename))
        if matches:
            pdb_id = matches[0].upper()

    # Fallback to filename if title is missing (Common in AlphaFold models)
    if not title_keyword and pdb_filename:
        # Extract meaningful keywords from filename
        # e.g., "orig_1M40_1_Chain_A_BETA-LACTAMASE_TEM_..."
        clean_name = os.path.basename(pdb_filename).replace("_", " ").replace("-", " ").upper()
        # Look for common resistance keywords
        for key in ["TEM", "NDM", "KPC", "OXA", "SHV", "VIM", "IMP", "LACTAMASE", "HYDROLASE"]:
            if key in clean_name:
                title_keyword = key
                break
        # If still no keyword, take the longest word that isn't the file extension or numeric
        if not title_keyword:
            words = [w for w in clean_name.split() if w and len(w) > 2 and not w.endswith('.PDB')]
            if words:
                title_keyword = max(words, key=len)

    return {
        "pdb_id": pdb_id, 
        "title_keyword": title_keyword, 
        "uniprot_id": uniprot_id
    }

def clean_pdb_file(input_path, output_path):
    """
    Cleans a PDB file by:
    1. Keeping only ATOM and HETATM records.
    2. Removing water (HOH) and other common solvents.
    3. Ensuring proper formatting.
    Returns: True if successful, False otherwise.
    """
    try:
        with open(input_path, 'r') as f:
            lines = f.readlines()
        
        cleaned_lines = []
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Filter out water and other common solvents
                res_name = line[17:20].strip().upper()
                if res_name in ['HOH', 'WAT', 'SOL', 'TIP3']:
                    continue
                cleaned_lines.append(line)
        
        if not cleaned_lines:
            return False
            
        with open(output_path, 'w') as f:
            f.writelines(cleaned_lines)
            
        return True
    except Exception as e:
        print(f"Error cleaning PDB file: {e}")
        return False


if __name__ == "__main__":
    # Test with sample content
    sample = """HEADER    HYDROLASE/HYDROLASE INHIBITOR           29-AUG-02   1M40
TITLE     ULTRA HIGH RESOLUTION CRYSTAL STRUCTURE OF TEM-1
TITLE    2 BETA-LACTAMASE"""
    print(extract_pdb_info(sample))
