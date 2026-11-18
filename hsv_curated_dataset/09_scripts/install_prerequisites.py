import os
import json
import urllib.request
import zipfile
import shutil

BIN_DIR = "hsv_curated_dataset/bin"
os.makedirs(f"{BIN_DIR}/hyphy", exist_ok=True)
os.makedirs(f"{BIN_DIR}/iqtree", exist_ok=True)

# 1. Check HyPhy releases
print("Checking HyPhy releases for Windows assets...")
req = urllib.request.Request("https://api.github.com/repos/veg/hyphy/releases", headers={"User-Agent": "Mozilla/5.0"})
try:
    with urllib.request.urlopen(req) as resp:
        releases = json.loads(resp.read().decode())
        found_hyphy = False
        for r in releases:
            for a in r.get("assets", []):
                name = a["name"].lower()
                if "win" in name or "exe" in name:
                    print(f"Found HyPhy Windows asset in {r['tag_name']}: {a['name']}")
                    found_hyphy = True
                    break
            if found_hyphy:
                break
        if not found_hyphy:
            print("No Windows binary found in recent HyPhy GitHub releases.")
except Exception as e:
    print(f"Error checking HyPhy releases: {e}")

# 2. Check IQ-TREE Windows download
print("\nDownloading IQ-TREE 2.4.0 for Windows...")
iqtree_url = "https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Windows.zip"
iqtree_zip = f"{BIN_DIR}/iqtree.zip"
try:
    urllib.request.urlretrieve(iqtree_url, iqtree_zip)
    print("Downloaded IQ-TREE. Extracting...")
    with zipfile.ZipFile(iqtree_zip, 'r') as zip_ref:
        zip_ref.extractall(f"{BIN_DIR}/iqtree")
    print("Extracted IQ-TREE.")
    # Find iqtree2.exe
    for root, dirs, files in os.walk(f"{BIN_DIR}/iqtree"):
        for f in files:
            if f.lower().startswith("iqtree2") and f.endswith(".exe"):
                print("Found iqtree2 binary at:", os.path.join(root, f))
except Exception as e:
    print(f"Error downloading IQ-TREE: {e}")
