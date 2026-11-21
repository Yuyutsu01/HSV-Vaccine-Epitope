"""
Step 11 Stage 1: Fast Evolutionary Selection Screening (FUBAR, aBSREL, BUSTED).
Automates Stage 1 selection analysis across:
  - Tier A (13 genes: unmasked)
  - Tier B (2 genes: masked + unmasked)
  - Tier C (5 genes: US4, US7, US8, UL49 masked + unmasked; UL46 unmasked)
Across 3 cohorts: HSV1_only, HSV2_only, combined.

Features:
  - Strict header sanitization matching IQ-TREE's mapping; original alignments untouched.
  - Generates temp alignment in 07_qc_reports/hyphy/_tmp/ and cleans up after run.
  - Preserves 07_qc_reports/hyphy/_tmp/<gene>_<cohort>_mapping.tsv.
  - Resumable: skips already completed JSON files.
  - Progress logging to 07_qc_reports/hyphy/stage1_progress.log.
  - Compiles 07_qc_reports/hyphy/stage1_summary.tsv and 07_qc_reports/hyphy/stage2_trigger_genes.tsv.
"""

import os
import sys
import time
import json
import subprocess
from datetime import datetime
import pandas as pd
from Bio import SeqIO

# Path configurations
WORKSPACE_DIR = os.path.abspath("hsv_curated_dataset").replace("\\", "/")
DIR_ALN_BASE = "hsv_curated_dataset/06_alignments"
DIR_TREES_BASE = "hsv_curated_dataset/06_alignments/trees"
DIR_REPORTS_BASE = "hsv_curated_dataset/07_qc_reports/hyphy"
DIR_TMP = os.path.join(DIR_REPORTS_BASE, "_tmp")
DIR_STAGE1 = os.path.join(DIR_REPORTS_BASE, "stage1")
LOG_PROGRESS = os.path.join(DIR_REPORTS_BASE, "stage1_progress.log")

os.makedirs(DIR_TMP, exist_ok=True)
os.makedirs(DIR_STAGE1, exist_ok=True)
for sub in ["tier_a", "tier_b", "tier_c"]:
    os.makedirs(os.path.join(DIR_STAGE1, sub), exist_ok=True)

# Tiers definition
TIER_A = [
    "US6", "UL27", "UL44", "UL22", "UL1", "UL53", "UL10",
    "UL49A", "UL54", "UL48", "UL47", "UL23", "UL30"
]
TIER_B = ["RL2", "RS1"]
TIER_C_MASKED = ["US4", "US7", "US8", "UL49"]
TIER_C_UNMASKED = ["UL46"]

COHORTS = ["HSV1", "HSV2", "combined"]
STAGE1_METHODS = ["FUBAR", "aBSREL", "BUSTED"]


def parse_iqtree_mapping(log_path: str) -> dict:
    """Parses IQ-TREE renaming map from the tree log file to synchronize tip names."""
    mapping = {}
    if os.path.exists(log_path):
        with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            in_block = False
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


def build_jobs_manifest():
    """Builds the comprehensive list of all Stage 1 jobs."""
    jobs = []

    # 1. Tier A: unmasked only (13 genes x 3 cohorts = 39)
    for g in TIER_A:
        for c in COHORTS:
            aln_file = f"{DIR_ALN_BASE}/{c}/{g}_aligned.fasta" if c != "combined" else f"{DIR_ALN_BASE}/combined/{g}_combined_aligned.fasta"
            tree_file = f"{DIR_TREES_BASE}/{c}/{g}.treefile"
            tree_log = f"{DIR_TREES_BASE}/{c}/{g}.log"
            jobs.append({
                "tier": "tier_a",
                "gene": g,
                "cohort": c,
                "is_masked": False,
                "aln_path": aln_file,
                "tree_path": tree_file,
                "tree_log": tree_log,
                "label": f"{g}_{c}"
            })

    # 2. Tier B: masked AND unmasked (2 genes x 3 cohorts x 2 = 12)
    for g in TIER_B:
        for c in COHORTS:
            tree_file = f"{DIR_TREES_BASE}/{c}/{g}.treefile"
            tree_log = f"{DIR_TREES_BASE}/{c}/{g}.log"
            # Unmasked
            aln_unmasked = f"{DIR_ALN_BASE}/{c}/{g}_aligned.fasta" if c != "combined" else f"{DIR_ALN_BASE}/combined/{g}_combined_aligned.fasta"
            jobs.append({
                "tier": "tier_b",
                "gene": g,
                "cohort": c,
                "is_masked": False,
                "aln_path": aln_unmasked,
                "tree_path": tree_file,
                "tree_log": tree_log,
                "label": f"{g}_{c}_unmasked"
            })
            # Masked
            aln_masked = f"{DIR_ALN_BASE}/masked/{c}/{g}_masked.fasta" if c != "combined" else f"{DIR_ALN_BASE}/masked/combined/{g}_combined_masked.fasta"
            jobs.append({
                "tier": "tier_b",
                "gene": g,
                "cohort": c,
                "is_masked": True,
                "aln_path": aln_masked,
                "tree_path": tree_file,
                "tree_log": tree_log,
                "label": f"{g}_{c}_masked"
            })

    # 3. Tier C: masked + unmasked for US4, US7, US8, UL49; unmasked only for UL46
    for g in TIER_C_MASKED:
        for c in COHORTS:
            tree_file = f"{DIR_TREES_BASE}/{c}/{g}.treefile"
            tree_log = f"{DIR_TREES_BASE}/{c}/{g}.log"
            aln_unmasked = f"{DIR_ALN_BASE}/{c}/{g}_aligned.fasta" if c != "combined" else f"{DIR_ALN_BASE}/combined/{g}_combined_aligned.fasta"
            jobs.append({
                "tier": "tier_c",
                "gene": g,
                "cohort": c,
                "is_masked": False,
                "aln_path": aln_unmasked,
                "tree_path": tree_file,
                "tree_log": tree_log,
                "label": f"{g}_{c}_unmasked"
            })
            # For Tier C masked files (if generated or unmasked baseline)
            aln_masked = f"{DIR_ALN_BASE}/masked/{c}/{g}_masked.fasta" if c != "combined" else f"{DIR_ALN_BASE}/masked/combined/{g}_combined_masked.fasta"
            if os.path.exists(aln_masked):
                jobs.append({
                    "tier": "tier_c",
                    "gene": g,
                    "cohort": c,
                    "is_masked": True,
                    "aln_path": aln_masked,
                    "tree_path": tree_file,
                    "tree_log": tree_log,
                    "label": f"{g}_{c}_masked"
                })

    for g in TIER_C_UNMASKED:
        for c in COHORTS:
            tree_file = f"{DIR_TREES_BASE}/{c}/{g}.treefile"
            tree_log = f"{DIR_TREES_BASE}/{c}/{g}.log"
            aln_unmasked = f"{DIR_ALN_BASE}/{c}/{g}_aligned.fasta" if c != "combined" else f"{DIR_ALN_BASE}/combined/{g}_combined_aligned.fasta"
            jobs.append({
                "tier": "tier_c",
                "gene": g,
                "cohort": c,
                "is_masked": False,
                "aln_path": aln_unmasked,
                "tree_path": tree_file,
                "tree_log": tree_log,
                "label": f"{g}_{c}"
            })

    return jobs


def run_stage1_pipeline():
    """Main execution orchestrator for Stage 1 screening."""
    jobs = build_jobs_manifest()
    total_runs = len(jobs) * len(STAGE1_METHODS)
    log_event(f"Initialized Stage 1 screening. Total target combinations: {len(jobs)} alignments x {len(STAGE1_METHODS)} methods = {total_runs} runs.")

    completed_count = 0
    skipped_count = 0
    error_count = 0

    for j in jobs:
        tier = j["tier"]
        gene = j["gene"]
        cohort = j["cohort"]
        aln_path = j["aln_path"]
        tree_path = j["tree_path"]
        tree_log = j["tree_log"]
        label = j["label"]

        if not os.path.exists(aln_path) or not os.path.exists(tree_path):
            log_event(f"SKIPPED {label}: missing alignment or tree file ({aln_path} / {tree_path})")
            continue

        # Refinement 1: Synchronize sanitized headers via IQ-TREE map
        mapping = parse_iqtree_mapping(tree_log)
        recs = list(SeqIO.parse(aln_path, "fasta"))

        tmp_fasta = os.path.join(DIR_TMP, f"{label}_tmp.fasta")
        mapping_tsv = os.path.join(DIR_TMP, f"{label}_mapping.tsv")

        # Write mapping table if not present
        if not os.path.exists(mapping_tsv):
            map_rows = []
            for r in recs:
                s_id = mapping.get(r.id, r.id)
                acc = r.id.split("|")[0] if "|" in r.id else r.id
                tier_info = "Tier1" if "cohort_tier=Tier1" in r.id else ("Tier2" if "cohort_tier=Tier2" in r.id else "Tier")
                map_rows.append({
                    "sanitized_id": s_id,
                    "original_id": r.id,
                    "accession": acc,
                    "cohort_tier": tier_info
                })
            pd.DataFrame(map_rows).to_csv(mapping_tsv, sep="\t", index=False)

        # Write sanitized temp alignment
        with open(tmp_fasta, "w", encoding="utf-8") as f:
            for r in recs:
                s_id = mapping.get(r.id, r.id)
                f.write(f">{s_id}\n{r.seq}\n")

        # Container paths
        rel_tmp_fasta = tmp_fasta.replace("hsv_curated_dataset", "/data").replace("\\", "/")
        rel_tree = tree_path.replace("hsv_curated_dataset", "/data").replace("\\", "/")

        for method in STAGE1_METHODS:
            out_json = os.path.join(DIR_STAGE1, tier, f"{label}_{method}.json")
            rel_out_json = out_json.replace("hsv_curated_dataset", "/data").replace("\\", "/")

            # Resumable skip check
            if os.path.exists(out_json) and os.path.getsize(out_json) > 500:
                skipped_count += 1
                completed_count += 1
                if completed_count % 25 == 0:
                    log_event(f"Progress checkpoint: {completed_count}/{total_runs} runs processed ({skipped_count} skipped/cached).")
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
                    log_event(f"SUCCESS [{completed_count+1}/{total_runs}]: {label} {method} completed in {runtime}s.")
                else:
                    error_count += 1
                    err_snippet = (res.stderr or res.stdout or "Unknown error")[:200].replace("\n", " ")
                    log_event(f"ERROR [{completed_count+1}/{total_runs}]: {label} {method} failed (code {res.returncode}): {err_snippet}")
            except Exception as e:
                error_count += 1
                log_event(f"EXCEPTION [{completed_count+1}/{total_runs}]: {label} {method} exception: {e}")

            completed_count += 1
            if completed_count % 25 == 0:
                log_event(f"Progress checkpoint: {completed_count}/{total_runs} runs processed.")

        # Clean up temp fasta per Refinement 1
        if os.path.exists(tmp_fasta):
            try:
                os.remove(tmp_fasta)
            except Exception:
                pass

    log_event(f"Stage 1 screening finished! Total runs: {completed_count}, Errors: {error_count}. Compiling summaries...")
    compile_stage1_summaries(jobs)


def compile_stage1_summaries(jobs: list):
    """Compiles Stage 1 summary table and triggers Stage 2 candidate list."""
    summary_rows = []
    trigger_rows = []

    for j in jobs:
        tier = j["tier"]
        gene = j["gene"]
        cohort = j["cohort"]
        label = j["label"]

        # Parse FUBAR
        fubar_json = os.path.join(DIR_STAGE1, tier, f"{label}_FUBAR.json")
        fubar_pp90 = 0
        fubar_pp95 = 0
        if os.path.exists(fubar_json):
            try:
                with open(fubar_json, "r") as f:
                    d = json.load(f)
                headers = [h[0] for h in d["MLE"]["headers"]]
                pp_idx = headers.index("Prob[alpha<beta]")
                for row in d["MLE"]["content"]["0"]:
                    if row[pp_idx] >= 0.90:
                        fubar_pp90 += 1
                    if row[pp_idx] >= 0.95:
                        fubar_pp95 += 1
            except Exception:
                pass

        # Parse aBSREL
        absrel_json = os.path.join(DIR_STAGE1, tier, f"{label}_aBSREL.json")
        absrel_branches = 0
        if os.path.exists(absrel_json):
            try:
                with open(absrel_json, "r") as f:
                    d = json.load(f)
                absrel_branches = d.get("test results", {}).get("positive test results", 0)
            except Exception:
                pass

        # Parse BUSTED
        busted_json = os.path.join(DIR_STAGE1, tier, f"{label}_BUSTED.json")
        busted_p = 1.0
        if os.path.exists(busted_json):
            try:
                with open(busted_json, "r") as f:
                    d = json.load(f)
                busted_p = round(float(d.get("test results", {}).get("p-value", 1.0)), 4)
            except Exception:
                pass

        any_signal = (fubar_pp90 >= 1) or (absrel_branches >= 1) or (busted_p < 0.05)
        trigger_reasons = []
        if fubar_pp90 >= 1:
            trigger_reasons.append(f"FUBAR_sites={fubar_pp90}")
        if absrel_branches >= 1:
            trigger_reasons.append(f"aBSREL_branches={absrel_branches}")
        if busted_p < 0.05:
            trigger_reasons.append(f"BUSTED_p={busted_p}")

        summary_rows.append({
            "gene": gene,
            "cohort": cohort,
            "tier": tier,
            "label": label,
            "fubar_sites_pp90": fubar_pp90,
            "fubar_sites_pp95": fubar_pp95,
            "absrel_branches": absrel_branches,
            "busted_p": busted_p,
            "any_signal": any_signal,
            "stage2_trigger": any_signal
        })

        if any_signal:
            trigger_rows.append({
                "gene": gene,
                "cohort": cohort,
                "tier": tier,
                "trigger_reason": "; ".join(trigger_reasons),
                "fubar_sites": fubar_pp90,
                "absrel_branches": absrel_branches,
                "busted_p": busted_p
            })

    # Save summary table
    df_summary = pd.DataFrame(summary_rows)
    out_summary = os.path.join(DIR_REPORTS_BASE, "stage1_summary.tsv")
    df_summary.to_csv(out_summary, sep="\t", index=False)
    log_event(f"Saved Stage 1 Summary to {out_summary}")

    # Save Stage 2 trigger table
    df_trigger = pd.DataFrame(trigger_rows)
    out_trigger = os.path.join(DIR_REPORTS_BASE, "stage2_trigger_genes.tsv")
    df_trigger.to_csv(out_trigger, sep="\t", index=False)
    log_event(f"Saved Stage 2 Trigger List to {out_trigger} (Total candidates: {len(df_trigger)})")


if __name__ == "__main__":
    run_stage1_pipeline()
