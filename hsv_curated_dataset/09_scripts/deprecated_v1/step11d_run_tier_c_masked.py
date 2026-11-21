"""
Step 11d: Execute Stage 1 Screening on Tier C Masked Alignments.
Runs FUBAR, aBSREL, and BUSTED on the 12 newly created masked alignments
for US4, US7, US8, and UL49 across HSV1, HSV2, and combined.
Total: 12 combinations x 3 methods = 36 runs.
"""

import os
import sys
import time
import json
import subprocess
import pandas as pd
from datetime import datetime
from Bio import SeqIO

DIR_BASE = "hsv_curated_dataset"
DIR_MASKED = f"{DIR_BASE}/06_alignments/masked"
DIR_TREES = f"{DIR_BASE}/06_alignments/trees"
DIR_STAGE1_TIER_C = f"{DIR_BASE}/07_qc_reports/hyphy/stage1/tier_c"
DIR_TMP = f"{DIR_BASE}/07_qc_reports/hyphy/_tmp"
LOG_FILE = f"{DIR_BASE}/07_qc_reports/hyphy/tier_c_masked_progress.log"
WORKSPACE_DIR = os.path.abspath(DIR_BASE)

os.makedirs(DIR_STAGE1_TIER_C, exist_ok=True)
os.makedirs(DIR_TMP, exist_ok=True)

GENES = ["US4", "US7", "US8", "UL49"]
COHORTS = ["HSV1", "HSV2", "combined"]
METHODS = ["FUBAR", "BUSTED", "aBSREL"]


def parse_iqtree_mapping(tree_log_path: str):
    """Extracts sequence ID renaming mapping from IQ-TREE log."""
    mapping = {}
    if not os.path.exists(tree_log_path):
        return mapping
    in_block = False
    with open(tree_log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "WARNING: Some sequence names are changed as follows:" in line:
                in_block = True
                continue
            if in_block:
                if " -> " in line:
                    orig, renamed = line.strip().split(" -> ")
                    mapping[orig] = renamed
                elif line.strip() == "":
                    break
    return mapping


def log_msg(msg: str):
    """Logs timestamped event to file and stdout."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted = f"[{ts}] {msg}"
    print(formatted, flush=True)
    with open(LOG_FILE, "a", encoding="utf-8") as fh:
        fh.write(formatted + "\n")


def run_tier_c_masked():
    log_msg("=== LAUNCHING STEP D: TIER C MASKED SCREENING (36 RUNS) ===")
    total_runs = len(GENES) * len(COHORTS) * len(METHODS)
    completed_count = 0
    skipped_count = 0
    error_count = 0

    for g in GENES:
        for c in COHORTS:
            aln_path = f"{DIR_MASKED}/{c}/{g}_masked.fasta" if c != "combined" else f"{DIR_MASKED}/combined/{g}_combined_masked.fasta"
            tree_path = f"{DIR_TREES}/{c}/{g}.treefile"
            tree_log = f"{DIR_TREES}/{c}/{g}.log"
            label = f"{g}_{c}_masked"

            if not os.path.exists(aln_path) or not os.path.exists(tree_path):
                log_msg(f"SKIPPED {label}: missing alignment or tree ({aln_path} / {tree_path})")
                continue

            mapping = parse_iqtree_mapping(tree_log)
            recs = list(SeqIO.parse(aln_path, "fasta"))

            tmp_fasta = os.path.join(DIR_TMP, f"{label}_tmp.fasta")
            mapping_tsv = os.path.join(DIR_TMP, f"{label}_mapping.tsv")

            if not os.path.exists(mapping_tsv):
                map_rows = []
                for r in recs:
                    s_id = mapping.get(r.id, r.id)
                    acc = r.id.split("|")[0] if "|" in r.id else r.id
                    map_rows.append({
                        "sanitized_id": s_id,
                        "original_id": r.id,
                        "accession": acc
                    })
                pd.DataFrame(map_rows).to_csv(mapping_tsv, sep="\t", index=False)

            with open(tmp_fasta, "w", encoding="utf-8") as fh:
                for r in recs:
                    s_id = mapping.get(r.id, r.id)
                    fh.write(f">{s_id}\n{r.seq}\n")

            rel_tmp_fasta = tmp_fasta.replace("hsv_curated_dataset", "/data").replace("\\", "/")
            rel_tree = tree_path.replace("hsv_curated_dataset", "/data").replace("\\", "/")

            for method in METHODS:
                out_json = os.path.join(DIR_STAGE1_TIER_C, f"{label}_{method}.json")
                rel_out_json = out_json.replace("hsv_curated_dataset", "/data").replace("\\", "/")

                if os.path.exists(out_json) and os.path.getsize(out_json) > 500:
                    skipped_count += 1
                    completed_count += 1
                    log_msg(f"CACHED [{completed_count}/{total_runs}]: {label} {method} already exists.")
                    continue

                t0 = time.time()
                cmd = [
                    "docker", "run", "--rm",
                    "--mount", f"type=bind,src={WORKSPACE_DIR},dst=/data",
                    "hyphy/hyphy",
                    "hyphy", method.lower(),
                    "--alignment", rel_tmp_fasta,
                    "--tree", rel_tree,
                    "--output", rel_out_json
                ]

                try:
                    res = subprocess.run(cmd, capture_output=True, text=True)
                    runtime = round(time.time() - t0, 1)

                    if res.returncode == 0 and os.path.exists(out_json):
                        completed_count += 1
                        log_msg(f"SUCCESS [{completed_count}/{total_runs}]: {label} {method} completed in {runtime}s.")
                    else:
                        error_count += 1
                        err_snip = (res.stderr or res.stdout or "Unknown error")[:200].replace("\n", " ")
                        log_msg(f"ERROR [{completed_count+1}/{total_runs}]: {label} {method} failed (code {res.returncode}): {err_snip}")
                except Exception as e:
                    error_count += 1
                    log_msg(f"EXCEPTION [{completed_count+1}/{total_runs}]: {label} {method} - {str(e)}")

            if os.path.exists(tmp_fasta):
                try:
                    os.remove(tmp_fasta)
                except Exception:
                    pass

    log_msg(f"=== STEP D FINISHED: {completed_count}/{total_runs} completed, {error_count} errors. ===")


if __name__ == "__main__":
    run_tier_c_masked()
