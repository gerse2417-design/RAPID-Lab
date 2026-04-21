import streamlit as st
import py3Dmol
from stmol import showmol
import os
import glob

# Page config
st.set_page_config(page_title="AlphaFold2 Structure Viewer", layout="wide")

st.title("🧬 AlphaFold2 Structure Viewer")
st.markdown("""
This tool visualizes the 3D structure predicted by AlphaFold2.
The coloring represents the **pLDDT (confidence score)**:
- 🔵 **Blue (>90)**: Very high confidence
- 🧊 **Light Blue (70-90)**: Confident
- 🟡 **Yellow (50-70)**: Low confidence
- 🟠 **Orange (<50)**: Very low confidence
""")

def render_pdb(pdb_file):
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()
    
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, 'pdb')
    
    # Apply AlphaFold coloring (pLDDT is stored in B-factor column)
    # Color schemes: Dark Blue (>90), Light Blue (80), Yellow (70), Orange (50)
    view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 50, 'max': 90}}})
    
    view.zoomTo()
    showmol(view, height=600, width=1000)

# Sidebar for file selection
st.sidebar.header("Settings")

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RESULTS_DIR = os.path.join(SCRIPT_DIR, "af2_results")

results_dir = st.sidebar.text_input("Results Directory", DEFAULT_RESULTS_DIR)

if os.path.exists(results_dir):
    pdb_files = glob.glob(os.path.join(results_dir, "*.pdb"))
    if pdb_files:
        # Sort to show ranked_0 (best) first
        pdb_files.sort() 
        selected_file = st.sidebar.selectbox("Select PDB File", pdb_files)
        
        if selected_file:
            st.subheader(f"Visualizing: {os.path.basename(selected_file)}")
            render_pdb(selected_file)
            
            # Show some stats if JSON exists
            json_file = selected_file.replace(".pdb", ".json")
            if os.path.exists(json_file):
                st.sidebar.info("Metadata found for this structure.")
    else:
        st.warning(f"No PDB files found in {results_dir}")
else:
    st.info(f"Directory {results_dir} not found. Please run the prediction first.")

st.sidebar.markdown("---")
st.sidebar.markdown("Built for AlphaFold2 Local Pipeline")
