"""
Prerequisite 2: Build Maximum-Likelihood Trees with IQ-TREE 2.
Ensures equal-length multiple alignments for all 20 genes across:
  - HSV1_only (20 trees)
  - HSV2_only (20 trees)
  - combined  (20 trees)
Outputs:
  - 06_alignments/trees/HSV1/{gene}.treefile
  - 06_alignments/trees/HSV2/{gene}.treefile
  - 06_alignments/trees/combined/{gene}.treefile
  - 07_qc_reports/11_prereq2_trees_summary.tsv
"""

import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio import Phylo

IQTREE_BIN = "hsv_curated_dataset/bin/iqtree/iqtree-2.4.0-Windows/bin/iqtree2.exe"

DIR_ALN_H1 = "hsv_curated_dataset/06_alignments/HSV1"
DIR_ALN_H2 = "hsv_curated_dataset/06_alignments/HSV2"
DIR_ALN_COMB = "hsv_curated_dataset/06_alignments/combined"

DIR_TREES_H1 = "hsv_curated_dataset/06_alignments/trees/HSV1"
DIR_TREES_H2 = "hsv_curated_dataset/06_alignments/trees/HSV2"
DIR_TREES_COMB = "hsv_curated_dataset/06_alignments/trees/combined"
DIR_REPORTS = "hsv_curated_dataset/07_qc_reports"

os.makedirs(DIR_TREES_H1, exist_ok=True)
os.makedirs(DIR_TREES_H2, exist_ok=True)
os.makedirs(DIR_TREES_COMB, exist_ok=True)
os.makedirs(DIR_REPORTS, exist_ok=True)

PRIORITY_20 = [
    "US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53", "UL10", "UL49A",
    "RL2", "RS1", "UL54", "UL48", "UL46", "UL47", "UL49", "UL23", "UL30"
]

print("=" * 70)
print("EXECUTING PREREQUISITE 2: MAXIMUM-LIKELIHOOD TREE INFERENCE (IQ-TREE 2)")
print("=" * 70)

# Step A: Validate and equalize per-type alignment lengths
# If per-type alignment has varying sequence lengths due to indels,
# extract the gap-padded coordinates directly from the combined alignment.
for gname in PRIORITY_20:
    comb_path = os.path.join(DIR_ALN_COMB, f"{gname}_combined_aligned.fasta")
    if not os.path.exists(comb_path):
        continue

    comb_recs = list(SeqIO.parse(comb_path, "fasta"))
    h1_subset = []
    h2_subset = []
    for r in comb_recs:
        clean_id = r.id
        if clean_id.startswith("HSV1|"):
            clean_id = clean_id[5:]
            h1_subset.append((clean_id, str(r.seq).upper()))
        elif clean_id.startswith("HSV2|"):
            clean_id = clean_id[5:]
            h2_subset.append((clean_id, str(r.seq).upper()))

    # Check HSV-1 per-type file
    h1_path = os.path.join(DIR_ALN_H1, f"{gname}_aligned.fasta")
    h1_existing = list(SeqIO.parse(h1_path, "fasta")) if os.path.exists(h1_path) else []
    h1_lens = set(len(r.seq) for r in h1_existing)
    if len(h1_lens) > 1 or not os.path.exists(h1_path):
        # Overwrite with clean equal-length alignment from combined
        with open(h1_path, "w", encoding="utf-8") as f:
            for r_id, s in h1_subset:
                f.write(f">{r_id}\n{s}\n")
        print(f"Equalized HSV1 {gname} alignment ({len(h1_subset)} seqs, {len(h1_subset[0][1])} nt)")

    # Check HSV-2 per-type file
    h2_path = os.path.join(DIR_ALN_H2, f"{gname}_aligned.fasta")
    h2_existing = list(SeqIO.parse(h2_path, "fasta")) if os.path.exists(h2_path) else []
    h2_lens = set(len(r.seq) for r in h2_existing)
    if len(h2_lens) > 1 or not os.path.exists(h2_path):
        # Overwrite with clean equal-length alignment from combined
        with open(h2_path, "w", encoding="utf-8") as f:
            for r_id, s in h2_subset:
                f.write(f">{r_id}\n{s}\n")
        print(f"Equalized HSV2 {gname} alignment ({len(h2_subset)} seqs, {len(h2_subset[0][1])} nt)")

# Step B: Run IQ-TREE on all 60 alignments
tree_records = []

for gname in PRIORITY_20:
    cohort_jobs = [
        ("HSV1", os.path.join(DIR_ALN_H1, f"{gname}_aligned.fasta"), os.path.join(DIR_TREES_H1, gname)),
        ("HSV2", os.path.join(DIR_ALN_H2, f"{gname}_aligned.fasta"), os.path.join(DIR_TREES_H2, gname)),
        ("combined", os.path.join(DIR_ALN_COMB, f"{gname}_combined_aligned.fasta"), os.path.join(DIR_TREES_COMB, gname)),
    ]

    for cohort, aln_path, out_prefix in cohort_jobs:
        treefile_path = f"{out_prefix}.treefile"

        cmd = [
            IQTREE_BIN,
            "-s", aln_path,
            "-m", "GTR+I+G4",
            "--fast",
            "-T", "4",
            "--prefix", out_prefix,
            "--seqtype", "DNA",
            "-redo"
        ]

        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if os.path.exists(treefile_path):
                tree = Phylo.read(treefile_path, "newick")
                total_len = round(float(tree.total_branch_length()), 4)

                # Inter-type divergence check on combined trees
                inter_div = "N/A"
                if cohort == "combined":
                    h1_tips = [t for t in tree.get_terminals() if t.name.startswith("HSV1|")]
                    h2_tips = [t for t in tree.get_terminals() if t.name.startswith("HSV2|")]
                    if h1_tips and h2_tips:
                        dist = tree.distance(h1_tips[0], h2_tips[0])
                        inter_div = f"{round(float(dist), 4)} subst/site"

                tree_records.append({
                    "gene": gname,
                    "cohort": cohort,
                    "status": "SUCCESS",
                    "total_tree_length": total_len,
                    "inter_type_divergence": inter_div,
                    "tree_file": treefile_path.replace("\\", "/"),
                    "notes": f"Tree built; length={total_len}"
                })
                print(f"[{cohort:8s}] {gname:6s} -> Tree successfully generated (total_len={total_len})")
            else:
                tree_records.append({
                    "gene": gname,
                    "cohort": cohort,
                    "status": "FAILED",
                    "total_tree_length": 0.0,
                    "inter_type_divergence": "N/A",
                    "tree_file": "N/A",
                    "notes": "treefile not written"
                })
                print(f"[{cohort:8s}] {gname:6s} -> FAILED")
        except Exception as e:
            tree_records.append({
                "gene": gname,
                "cohort": cohort,
                "status": "ERROR",
                "total_tree_length": 0.0,
                "inter_type_divergence": "N/A",
                "tree_file": "N/A",
                "notes": str(e)
            })
            print(f"[{cohort:8s}] {gname:6s} -> ERROR: {e}")

df_trees = pd.DataFrame(tree_records)
out_tsv = os.path.join(DIR_REPORTS, "11_prereq2_trees_summary.tsv")
df_trees.to_csv(out_tsv, sep="\t", index=False)
print(f"\nSaved Tree Summary to {out_tsv}")

# Group summary
summary = df_trees.groupby(["cohort", "status"]).size().unstack(fill_value=0)
print("\n=== Tree Generation Summary by Cohort ===")
print(summary)
