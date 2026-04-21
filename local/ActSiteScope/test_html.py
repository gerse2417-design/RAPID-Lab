import urllib.request
import time
import sys

time.sleep(2)  # Wait for streamlit to run
try:
    with urllib.request.urlopen('http://localhost:8501') as response:
        html = response.read().decode('utf-8')
        with open('streamlit_dump.html', 'w') as f:
            f.write(html)
        print("Success")
except Exception as e:
    print(f"Error: {e}")
