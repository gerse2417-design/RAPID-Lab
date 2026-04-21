import sys, json
sys.path.append('.')
from app.streamlit_app import render_3d_view
import unittest.mock as mock

with mock.patch('app.streamlit_app.components.html') as mock_html:
    render_3d_view("ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N")
    html_code = mock_html.call_args[0][0]
    with open('output_html_dump.html', 'w') as f:
        f.write(html_code)
print("Dumped HTML.")
