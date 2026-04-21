import re

file_path = "app/streamlit_app.py"
with open(file_path, "r") as f:
    content = f.read()

# Extract the script tag content
scripts = re.findall(r'<script>(.*?)</script>', content, re.DOTALL)
if scripts:
    script_content = scripts[0]
    # Replace VAR_... placeholders with something valid
    script_content = script_content.replace("VAR_DATA_JSON", "{}")
    script_content = script_content.replace("VAR_VIEWER_BG", "'white'")
    script_content = script_content.replace("VAR_TEXT_COLOR", "'black'")
    # ... and so on if needed
    
    with open("temp_check.js", "w") as f:
        f.write(script_content)
    print("Extracted script to temp_check.js")
else:
    print("No script found")
