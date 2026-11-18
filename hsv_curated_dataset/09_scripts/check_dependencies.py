"""
Check system dependencies and Python packages for HSV curated dataset pipeline.
"""
import sys
import shutil

print("=== Python Environment ===")
print("Python version:", sys.version)

packages = ["Bio", "pandas", "numpy", "scipy"]
for pkg in packages:
    try:
        mod = __import__(pkg)
        ver = getattr(mod, "__version__", "installed")
        print(f"  {pkg:10s}: {ver}")
    except ImportError:
        print(f"  {pkg:10s}: NOT INSTALLED")

print("\n=== External CLI Tools ===")
tools = ["datasets", "dataformat", "seqkit", "mafft", "macse", "blastn"]
for t in tools:
    path = shutil.which(t)
    if path:
        print(f"  {t:10s}: Found at {path}")
    else:
        print(f"  {t:10s}: NOT on PATH")
