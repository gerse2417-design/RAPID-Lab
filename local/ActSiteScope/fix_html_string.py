import sys, re

file_path = "app/streamlit_app.py"
with open(file_path, "r") as f:
    content = f.read()

# Define the HTML block
start_marker = "    html_code = f\"\"\""
end_marker = "    \"\"\"\n    components.html(html_code, height=600)"

if start_marker in content and end_marker in content:
    idx_start = content.find(start_marker)
    idx_end = content.find(end_marker) + len("    \"\"\"")
    
    html_section = content[idx_start + len(start_marker):content.find(end_marker)]
    
    # Unescape all {{ and }} back to { and }
    html_section = html_section.replace("{{", "{").replace("}}", "}")
    
    # Replace variables with VAR_
    html_section = html_section.replace("{data_json}", "VAR_DATA_JSON")
    html_section = html_section.replace("{viewer_bg}", "VAR_VIEWER_BG")
    html_section = html_section.replace("{text_color}", "VAR_TEXT_COLOR")
    html_section = html_section.replace("{panel_bg}", "VAR_PANEL_BG")
    html_section = html_section.replace("{panel_bdr}", "VAR_PANEL_BDR")
    html_section = html_section.replace("{shadow}", "VAR_SHADOW")
    html_section = html_section.replace("{legend_display}", "VAR_LEGEND_DISPLAY")

    # Replace the section block
    new_start = """    html_code = \"\"\"""" + html_section + """    \"\"\"
    html_code = html_code.replace("VAR_DATA_JSON", data_json)
    html_code = html_code.replace("VAR_VIEWER_BG", viewer_bg)
    html_code = html_code.replace("VAR_TEXT_COLOR", text_color)
    html_code = html_code.replace("VAR_PANEL_BG", panel_bg)
    html_code = html_code.replace("VAR_PANEL_BDR", panel_bdr)
    html_code = html_code.replace("VAR_SHADOW", shadow)
    html_code = html_code.replace("VAR_LEGEND_DISPLAY", legend_display)
    
    components.html(html_code, height=600)"""

    new_content = content[:idx_start] + new_start + content[content.find(end_marker) + len(end_marker):]
    
    with open(file_path, "w") as f:
        f.write(new_content)
    print("SUCCESS")
else:
    print("MARKERS NOT FOUND")

