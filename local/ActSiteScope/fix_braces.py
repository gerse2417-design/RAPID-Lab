import sys, os
file_path = "app/streamlit_app.py"

with open(file_path, "r") as f:
    lines = f.readlines()

def fix_line(line):
    # If the line already has double braces, don't double them again.
    # Actually, the safest way is to just replace all { and } with {{ and }}
    # EXCEP the ones already doubled.
    import re
    # Temporary replace {{ and }} to unique tokens
    line = line.replace("{{", "@@LBR@@").replace("}}", "@@RBR@@")
    # Replace single { and } with {{ and }}
    line = line.replace("{", "{{").replace("}", "}}")
    # Restore {{ and }}
    line = line.replace("@@LBR@@", "{{").replace("@@RBR@@", "}}")
    return line

for i in range(459, 565):  # 0-indexed: lines 460 to 565
    lines[i] = fix_line(lines[i])

with open(file_path, "w") as f:
    f.writelines(lines)

