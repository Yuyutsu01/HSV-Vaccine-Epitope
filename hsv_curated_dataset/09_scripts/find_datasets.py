"""
Find NCBI Datasets Windows download URL
"""
import urllib.request
import re

base_url = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/"
try:
    req = urllib.request.Request(base_url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        content = resp.read().decode("utf-8", errors="ignore")
        for match in re.finditer(r'href=["\']([^"\']+)["\']', content):
            href = match.group(1)
            print("Found:", href)
except Exception as e:
    print("Error:", e)
