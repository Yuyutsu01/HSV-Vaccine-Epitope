"""
Step 11: Evolutionary Selection Analysis Pipeline (Tier A, Tier B, Tier C).
Integrates HyPhy analyses (BUSTED, FEL, MEME, FUBAR, aBSREL) across:
  - Combined (HSV-1 + HSV-2) alignments
  - HSV-1 only alignments
  - HSV-2 only alignments

Tier A (13 clean genes):
  US6, UL27, UL44, UL22, UL1, UL53, UL10, UL49A, UL54, UL48, UL47, UL23, UL30
Tier C (5 flagged genes with low complexity / high divergence):
  US4, US7, US8, UL46, UL49
Tier B (2 repeat-rich genes):
  RL2, RS1 (Masking investigation in step11b_masking_investigation.py)
"""

import os
import sys
import json
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO

# Target directories
DIR_ALN_H1 = "hsv_curated_dataset/06_alignments/HSV1"
DIR_ALN_H2 = "hsv_curated_dataset/06_alignments/HSV2"
DIR_ALN_COMB = "hsv_curated_dataset/06_alignments/combined"
DIR_MASKED_COMB = "hsv_curated_dataset/06_alignments/masked/combined"
DIR_REPORTS = "hsv_curated_dataset/07_qc_reports"
DIR_HYPHY_RAW = "hsv_curated_dataset/10_hyphy_outputs"

os.makedirs(DIR_REPORTS, exist_ok=True)
os.makedirs(DIR_HYPHY_RAW, exist_ok=True)

TIER_A_GENES = [
    "US6", "UL27", "UL44", "UL22", "UL1", "UL53", "UL10",
    "UL49A", "UL54", "UL48", "UL47", "UL23", "UL30"
]

TIER_C_GENES = [
    "US4", "US7", "US8", "UL46", "UL49"
]

TIER_B_GENES = [
    "RL2", "RS1"
]

METHODS = ["BUSTED", "FEL", "MEME", "FUBAR", "aBSREL"]

# Pre-computed within-type identities from Step 10 diagnostic
WITHIN_TYPE_ID = {
    "US4": (69.02, 98.76),
    "US7": (73.57, 93.00),
    "US8": (79.61, 78.39),
    "UL46": (86.26, 98.96),
    "UL49": (97.98, 72.32),
}

# Pre-computed low-complexity coordinates for Tier C (SEG / entropy < 2.2 bits)
TIER_C_LC_REGIONS = {
    "US4": {"HSV1": [(77, 98)], "HSV2": [(438, 486)]},
    "US7": {"HSV1": [(207, 231)], "HSV2": [(205, 226)]},
    "US8": {"HSV1": [], "HSV2": [(183, 207), (411, 434)]},
    "UL46": {"HSV1": [], "HSV2": []},
    "UL49": {"HSV1": [(94, 113)], "HSV2": [(93, 113), (268, 297)]},
}


def check_hyphy_installed() -> bool:
    """Checks whether the HyPhy CLI binary is available on system PATH."""
    return shutil.which("hyphy") is not None


def run_hyphy_method(method: str, alignment_path: str, output_json: str, tree_path: str = None) -> dict:
    """
    Executes a HyPhy method CLI command on a codon alignment.
    Example: hyphy busted --alignment <aln> --output <json>
    """
    if not check_hyphy_installed():
        return {"status": "SKIPPED", "error": "hyphy executable not found on PATH"}

    cmd = ["hyphy", method.lower(), "--alignment", alignment_path, "--output", output_json]
    if tree_path and os.path.exists(tree_path):
        cmd.extend(["--tree", tree_path])

    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if os.path.exists(output_json):
            with open(output_json, "r") as f:
                return json.load(f)
        return {"status": "COMPLETED", "stdout": res.stdout}
    except Exception as e:
        return {"status": "ERROR", "error": str(e)}


def generate_tier_c_annotation_summary():
    """Generates the low-complexity annotation table for Tier C genes."""
    rows = []
    for g in TIER_C_GENES:
        id_h1, id_h2 = WITHIN_TYPE_ID[g]
        lc_h1 = "; ".join([f"{s}-{e}" for s, e in TIER_C_LC_REGIONS[g]["HSV1"]]) or "None"
        lc_h2 = "; ".join([f"{s}-{e}" for s, e in TIER_C_LC_REGIONS[g]["HSV2"]]) or "None"
        rows.append({
            "gene": g,
            "within_type_AA_identity_HSV1": id_h1,
            "within_type_AA_identity_HSV2": id_h2,
            "hsv1_lc_regions_aa": lc_h1,
            "hsv2_lc_regions_aa": lc_h2,
            "notes": "Low-complexity repeat tracts annotated; flagged for conservative site interpretation"
        })
    df = pd.DataFrame(rows)
    out_path = os.path.join(DIR_REPORTS, "11c_tier_c_lc_annotation.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved Tier C Low-Complexity Annotation Summary to {out_path}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    print("=" * 70)
    print("STEP 11: HYPHY SELECTION ANALYSIS HARNESS")
    print("=" * 70)
    print(f"HyPhy installed on system PATH: {check_hyphy_installed()}")
    generate_tier_c_annotation_summary()
