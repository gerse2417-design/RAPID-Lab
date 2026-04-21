import sys, os
sys.path.append('.')
from app.streamlit_app import render_3d_view
try:
    render_3d_view("ATOM")
    print("SUCCESS")
except Exception as e:
    print(f"ERROR: {e}")
