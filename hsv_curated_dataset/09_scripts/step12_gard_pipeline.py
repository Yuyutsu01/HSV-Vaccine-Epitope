"""
Step 12 GARD: Genetic Algorithm for Recombination Detection (GARD).
Screens the top-10 combined vaccine & regulatory genes for evidence of
phylogenetic discordance and recombination breakpoints.

Outputs:
  - 07_qc_reports/hyphy/gard/<gene>_combined_GARD.json
  - 07_qc_reports/hyphy/gard_progress.log
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
DIR_ALN_BASE = f"{DIR_BASE}/06_alignments"
DIR_TREES_BASE = f"{DIR_BASE}/06_alignments/trees"
DIR_GARD = f"{DIR_BASE}/07_qc_reports/hyphy/gard"
DIR_TMP = f"{DIR_BASE}/07_qc_reports/hyphy/_tmp"
LOG_PROGRESS = f"{DIR_BASE}/07_qc_reports/hyphy/gard_progress.log"
WORKSPACE_DIR = os.path.abspath(DIR_BASE)

os.makedirs(DIR_GARD, exist_ok=True)
os.makedirs(DIR_TMP, exist_ok=True)

TOP10_GENES = [
    "UL27",  # gB (primary entry fusogen)
    "US6",   # gD (receptor-binding subunit vaccine target)
    "UL22",  # gH (core fusion heterodimer partner)
    "UL1",   # gL (core fusion partner with gH)
    "UL44",  # gC (complement C3b receptor)
    "UL10",  # gM (structural envelope protein)
    "UL53",  # gK (syncytial entry regulator)
    "UL30",  # DNA polymerase catalytic subunit
    "UL47",  # VP13/14 major tegument protein
    "UL46",  # VP11/12 tegument phosphoprotein
]


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


def log_event(msg: str):
    """Appends timestamped event to progress log and prints to stdout."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted = f"[{ts}] {msg}"
    print(formatted, flush=True)
    with open(LOG_PROGRESS, "a", encoding="utf-8") as lf:
        lf.write(formatted + "\n")


def run_gard_pipeline():
    """Executes HyPhy GARD recombination screening across top-10 combined genes."""
    total_runs = len(TOP10_GENES)

    log_event("=" * 70)
    log_event(f"STARTING HYPHY GARD RECOMBINATION SCREENING ({total_runs} top genes)")
    log_event("=" * 70)

    completed_count = 0
    skipped_count = 0
    error_count = 0

    for idx, gene in enumerate(TOP10_GENES):
        aln_path = f"{DIR_ALN_BASE}/combined/{gene}_combined_aligned.fasta"
        tree_path = f"{DIR_TREES_BASE}/combined/{gene}.treefile"
        tree_log = f"{DIR_TREES_BASE}/combined/{gene}.log"
        label = f"{gene}_combined"

        if not os.path.exists(aln_path):
            log_event(f"SKIPPED {label}: Missing alignment ({aln_path})")
            continue

        mapping = parse_iqtree_mapping(tree_log)
        recs = list(SeqIO.parse(aln_path, "fasta"))

        tmp_fasta = os.path.join(DIR_TMP, f"{label}_gard_tmp.fasta")
        with open(tmp_fasta, "w", encoding="utf-8") as fh:
            for r in recs:
                s_id = mapping.get(r.id, r.id)
                fh.write(f">{s_id}\n{r.seq}\n")

        rel_tmp_fasta = tmp_fasta.replace("hsv_curated_dataset", "/data").replace("\\", "/")
        out_json = os.path.join(DIR_GARD, f"{label}_GARD.json")
        rel_out_json = out_json.replace("hsv_curated_dataset", "/data").replace("\\", "/")

        if os.path.exists(out_json) and os.path.getsize(out_json) > 500:
            skipped_count += 1
            completed_count += 1
            log_event(f"CACHED [{completed_count}/{total_runs}]: {label} GARD already exists.")
            continue

        t0 = time.time()
        cmd = [
            "docker", "run", "--rm",
            "--mount", f"type=bind,src={WORKSPACE_DIR},dst=/data",
            "hyphy/hyphy",
            "hyphy",
            "CPU=8",
            "ENV=TOLERATE_NUMERICAL_ERRORS=1;",
            "gard",
            "--alignment", rel_tmp_fasta,
            "--output", rel_out_json
        ]

        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            runtime = round(time.time() - t0, 1)

            if res.returncode == 0 and os.path.exists(out_json):
                completed_count += 1
                log_event(f"SUCCESS [{completed_count}/{total_runs}]: {label} GARD completed in {runtime}s.")
            else:
                error_count += 1
                err_snip = (res.stderr or res.stdout or "Unknown error")[:200].replace("\n", " ")
                log_event(f"ERROR [{completed_count+1}/{total_runs}]: {label} GARD failed (code {res.returncode}): {err_snip}")
        except Exception as e:
            error_count += 1
            log_event(f"EXCEPTION [{completed_count+1}/{total_runs}]: {label} GARD - {str(e)}")

        if os.path.exists(tmp_fasta):
            try:
                os.remove(tmp_fasta)
            except Exception:
                pass

    log_event("=" * 70)
    log_event(f"GARD RECOMBINATION SCREENING COMPLETE: {completed_count}/{total_runs} finished, {error_count} errors, {skipped_count} cached.")
    log_event("=" * 70)


if __name__ == "__main__":
    run_gard_pipeline()
