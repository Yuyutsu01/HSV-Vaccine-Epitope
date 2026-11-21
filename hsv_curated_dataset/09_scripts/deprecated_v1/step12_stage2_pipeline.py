"""
Step 12: Stage 2 Deep Site-Level Selection Analysis (FEL and MEME).
Executes FEL (Fixed Effects Likelihood) and MEME (Mixed Effects Model of Evolution)
on the ~25-30 triggered targets prioritized in stage2_trigger_genes_final.tsv.

Outputs:
  - 07_qc_reports/hyphy/stage2/<label>_FEL.json
  - 07_qc_reports/hyphy/stage2/<label>_MEME.json
  - 07_qc_reports/hyphy/stage2_progress.log
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
DIR_STAGE2 = f"{DIR_BASE}/07_qc_reports/hyphy/stage2"
DIR_TMP = f"{DIR_BASE}/07_qc_reports/hyphy/_tmp"
LOG_PROGRESS = f"{DIR_BASE}/07_qc_reports/hyphy/stage2_progress.log"
TRIGGER_TSV = f"{DIR_BASE}/07_qc_reports/hyphy/stage2_trigger_genes_final.tsv"
WORKSPACE_DIR = os.path.abspath(DIR_BASE)

STAGE2_METHODS = ["FEL", "MEME"]

os.makedirs(DIR_STAGE2, exist_ok=True)
os.makedirs(DIR_TMP, exist_ok=True)


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


def run_stage2_pipeline():
    """Main execution orchestrator for Stage 2 (FEL + MEME)."""
    if not os.path.exists(TRIGGER_TSV):
        log_event(f"Error: {TRIGGER_TSV} does not exist. Run Stage 1 trigger definition first.")
        return

    df_trigger = pd.read_csv(TRIGGER_TSV, sep="\t")
    total_runs = len(df_trigger) * len(STAGE2_METHODS)

    log_event("=" * 70)
    log_event(f"STARTING HYPHY STAGE 2: FEL & MEME DEEP DIVERSIFYING SELECTION ({total_runs} total runs)")
    log_event("=" * 70)

    completed_count = 0
    skipped_count = 0
    error_count = 0

    for idx, row in df_trigger.iterrows():
        gene = row["gene"]
        cohort = row["cohort"]
        tier = row["tier"]
        label = row["label"]

        # Resolve alignment path
        if "masked" in label:
            aln_path = f"{DIR_ALN_BASE}/masked/{cohort}/{gene}_masked.fasta" if cohort != "combined" else f"{DIR_ALN_BASE}/masked/combined/{gene}_combined_masked.fasta"
        else:
            aln_path = f"{DIR_ALN_BASE}/{cohort}/{gene}_aligned.fasta" if cohort != "combined" else f"{DIR_ALN_BASE}/combined/{gene}_combined_aligned.fasta"

        tree_path = f"{DIR_TREES_BASE}/{cohort}/{gene}.treefile"
        tree_log = f"{DIR_TREES_BASE}/{cohort}/{gene}.log"

        if not os.path.exists(aln_path) or not os.path.exists(tree_path):
            log_event(f"SKIPPED {label}: Missing alignment or tree ({aln_path} / {tree_path})")
            continue

        mapping = parse_iqtree_mapping(tree_log)
        recs = list(SeqIO.parse(aln_path, "fasta"))

        tmp_fasta = os.path.join(DIR_TMP, f"{label}_stage2_tmp.fasta")
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

        for method in STAGE2_METHODS:
            out_json = os.path.join(DIR_STAGE2, f"{label}_{method}.json")
            rel_out_json = out_json.replace("hsv_curated_dataset", "/data").replace("\\", "/")

            if os.path.exists(out_json) and os.path.getsize(out_json) > 500:
                skipped_count += 1
                completed_count += 1
                log_event(f"CACHED [{completed_count}/{total_runs}]: {label} {method} already exists.")
                continue

            t0 = time.time()
            cmd = [
                "docker", "run", "--rm",
                "--mount", f"type=bind,src={WORKSPACE_DIR},dst=/data",
                "hyphy/hyphy",
                "hyphy", "CPU=8", method.lower(),
                "--alignment", rel_tmp_fasta,
                "--tree", rel_tree,
                "--output", rel_out_json
            ]

            try:
                res = subprocess.run(cmd, capture_output=True, text=True)
                runtime = round(time.time() - t0, 1)

                if res.returncode == 0 and os.path.exists(out_json):
                    completed_count += 1
                    log_event(f"SUCCESS [{completed_count}/{total_runs}]: {label} {method} completed in {runtime}s.")
                else:
                    error_count += 1
                    err_snip = (res.stderr or res.stdout or "Unknown error")[:200].replace("\n", " ")
                    log_event(f"ERROR [{completed_count+1}/{total_runs}]: {label} {method} failed (code {res.returncode}): {err_snip}")
            except Exception as e:
                error_count += 1
                log_event(f"EXCEPTION [{completed_count+1}/{total_runs}]: {label} {method} - {str(e)}")

        if os.path.exists(tmp_fasta):
            try:
                os.remove(tmp_fasta)
            except Exception:
                pass

    log_event("=" * 70)
    log_event(f"STAGE 2 EXECUTION COMPLETE: {completed_count}/{total_runs} finished, {error_count} errors, {skipped_count} cached.")
    log_event("=" * 70)


if __name__ == "__main__":
    run_stage2_pipeline()

