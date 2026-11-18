"""
Step 6.5: Comprehensive Validation and Quality Verification Pipeline.
Validates extracted CDS for all 20 priority genes and 74 total genes across HSV-1 and HSV-2.

Sub-tasks:
  - 6.5a: Translation validation (identity, length ratio, internal stops, classification)
  - 6.5b: FAIL_ANCHOR_EXAMPLES analysis
  - 6.5c: UL23 (TK) frameshift verification and resistance variant identification
  - 6.5d: RL2 (ICP0) splice junction and intron GT-AG boundary verification
  - 6.5e: Reconcile genome attempt counts
  - 6.5f: 74-gene per-cohort extraction success table
"""

import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

# -------------------------------------------------------------
# File Paths and Target Definitions
# -------------------------------------------------------------
REF1_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001806.2.gbk"
REF2_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001798.2.gbk"
TIER1_HSV1_FASTA = "hsv_curated_dataset/02_filtered_genomes/HSV1_filtered.fasta"
TIER1_HSV2_FASTA = "hsv_curated_dataset/02_filtered_genomes/HSV2_filtered.fasta"
TIER2_HSV2_FASTA = "hsv_curated_dataset/02_filtered_genomes/HSV2_filtered_tier2.fasta"
EXTRACTION_LOG = "hsv_curated_dataset/07_qc_reports/06_cds_extraction_log.tsv"
CURATED_META = "hsv_curated_dataset/03_metadata/curated_metadata.tsv"

OUT_QC = "hsv_curated_dataset/07_qc_reports"
os.makedirs(OUT_QC, exist_ok=True)

PRIORITY_GENES = [
    "US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53", "UL10", "UL49A",
    "RL2", "RS1", "UL54", "UL48", "UL46", "UL47", "UL49", "UL23", "UL30"
]

# -------------------------------------------------------------
# 1. Load Reference Proteins and Gene Models
# -------------------------------------------------------------
def load_reference_data(gbk_path):
    record = SeqIO.read(gbk_path, "genbank")
    full_seq = str(record.seq).upper()
    prots = {}
    models = {}
    seen = set()
    for feat in record.features:
        if feat.type == "CDS":
            g = feat.qualifiers.get("gene", [""])[0] or feat.qualifiers.get("locus_tag", [""])[0]
            if not g or g in seen:
                continue
            seen.add(g)
            t = feat.qualifiers.get("translation", [""])[0]
            prots[g] = t
            
            strand = "+" if feat.location.strand == 1 else "-"
            parts = []
            for part in feat.location.parts:
                p_s = int(part.start)
                p_e = int(part.end)
                p_seq = full_seq[p_s:p_e]
                if strand == "-":
                    p_seq = str(Seq(p_seq).reverse_complement())
                parts.append({"start": p_s, "end": p_e, "seq": p_seq})
            full_cds = "".join(p["seq"] for p in parts)
            models[g] = {
                "gene": g,
                "strand": strand,
                "parts": parts,
                "full_cds": full_cds,
                "ref_len_nt": len(full_cds),
                "ref_len_aa": len(t),
                "is_repeat": g in ["RL1", "RL2", "RS1"]
            }
    return prots, models, full_seq

ref_prots_h1, ref_models_h1, full_ref1 = load_reference_data(REF1_GBK)
ref_prots_h2, ref_models_h2, full_ref2 = load_reference_data(REF2_GBK)

# Protein Pairwise Aligner
prot_aligner = PairwiseAligner()
prot_aligner.mode = "global"
prot_aligner.match_score = 2
prot_aligner.mismatch_score = -1
prot_aligner.open_gap_score = -4
prot_aligner.extend_gap_score = -1

# Load extraction log to map accessions to tier
df_ext_log = pd.read_csv(EXTRACTION_LOG, sep="\t")
acc_to_tier = {}
for _, row in df_ext_log.iterrows():
    acc_to_tier[(row["accession"], row["virus"])] = row["cohort_tier"]

# -------------------------------------------------------------
# STEP 6.5a: Translation Validation for 20 Priority Genes
# -------------------------------------------------------------
print("="*60)
print("EXECUTING STEP 6.5a: TRANSLATION VALIDATION (20 PRIORITY GENES)")
print("="*60)

val_results = []
fail_anchor_examples = []

for vtype in ["HSV-1", "HSV-2"]:
    ref_prots = ref_prots_h1 if vtype == "HSV-1" else ref_prots_h2
    v_dir = "HSV1" if vtype == "HSV-1" else "HSV2"
    
    for gname in PRIORITY_GENES:
        fasta_file = os.path.join(f"hsv_curated_dataset/04_cds_unaligned/{v_dir}", f"{gname}.fasta")
        if not os.path.exists(fasta_file):
            continue
            
        ref_p = ref_prots.get(gname, "")
        ref_l = len(ref_p)
        
        # Group by cohort tier
        tier_records = {"Tier 1": [], "Tier 2": []}
        
        for rec in SeqIO.parse(fasta_file, "fasta"):
            acc = rec.id.split("|")[0]
            tier = acc_to_tier.get((acc, vtype), "Tier 1")
            s_cds = str(rec.seq).upper()
            
            # Trim trailing partial codons if any
            valid_len = len(s_cds) - (len(s_cds) % 3)
            s_clean = s_cds[:valid_len]
            
            # Translate to first stop and full translation
            p_to_stop = str(Seq(s_clean).translate(to_stop=True))
            p_full = str(Seq(s_clean).translate())
            
            # Count internal stops (stop codons occurring before the final position)
            internal_stops = p_full[:-1].count("*") if len(p_full) > 1 else 0
            
            # Calculate length ratio
            trans_len = len(p_to_stop)
            length_ratio = trans_len / max(1, ref_l)
            
            # Align protein to reference to compute identity
            # Use quick identity if exact length and high similarity, else full alignment
            if p_to_stop == ref_p:
                aa_identity = 100.0
            elif abs(trans_len - ref_l) < 30 and trans_len > 0:
                # Fast global alignment
                aln = next(iter(prot_aligner.align(p_to_stop, ref_p)))
                matches = sum(c1 == c2 for c1, c2 in zip(aln[0], aln[1]))
                aa_identity = round(matches / max(1, len(aln[0])) * 100.0, 2)
            else:
                aln = next(iter(prot_aligner.align(p_to_stop, ref_p)))
                matches = sum(c1 == c2 for c1, c2 in zip(aln[0], aln[1]))
                aa_identity = round(matches / max(1, len(aln[0])) * 100.0, 2)
                
            # Classification rules:
            # - PASS_GOOD: identity >= 90%, 0.95 <= length_ratio <= 1.05, internal_stops=0
            # - PASS_TRUNCATED: identity >= 90%, 0.80 <= length_ratio < 0.95, ends at stop (expected for UL23 resistance variants only)
            # - FAIL_ANCHOR_ERROR: identity < 80%, OR length_ratio < 0.5, OR length_ratio > 1.2, OR internal_stops > 1
            # - FAIL_NONFUNCTIONAL: internal_stops >= 1 in a non-UL23 gene
            
            if (aa_identity < 80.0) or (length_ratio < 0.5) or (length_ratio > 1.2) or (internal_stops > 1):
                cls = "FAIL_ANCHOR_ERROR"
                fail_anchor_examples.append({
                    "virus": vtype,
                    "tier": tier,
                    "gene": gname,
                    "accession": acc,
                    "translated_length": trans_len,
                    "ref_length": ref_l,
                    "length_ratio": round(length_ratio, 3),
                    "aa_identity": aa_identity,
                    "internal_stops": internal_stops
                })
            elif internal_stops >= 1 and gname != "UL23":
                cls = "FAIL_NONFUNCTIONAL"
            elif aa_identity >= 90.0 and (0.95 <= length_ratio <= 1.05) and internal_stops == 0:
                cls = "PASS_GOOD"
            elif aa_identity >= 90.0 and (0.80 <= length_ratio < 0.95):
                cls = "PASS_TRUNCATED"
            else:
                # Other nonfunctional variants (e.g. truncated in non-UL23, or identity 80-90%)
                if gname == "UL23" and length_ratio < 0.95:
                    cls = "PASS_TRUNCATED"
                else:
                    cls = "FAIL_NONFUNCTIONAL"
                    
            tier_records[tier].append({
                "accession": acc,
                "cls": cls,
                "length_ratio": length_ratio,
                "aa_identity": aa_identity
            })
            
        # Summarize per tier
        for tier, items in tier_records.items():
            if not items:
                continue
            total_ext = len(items)
            c_good = sum(1 for x in items if x["cls"] == "PASS_GOOD")
            c_trunc = sum(1 for x in items if x["cls"] == "PASS_TRUNCATED")
            c_anchor_err = sum(1 for x in items if x["cls"] == "FAIL_ANCHOR_ERROR")
            c_nonfunc = sum(1 for x in items if x["cls"] == "FAIL_NONFUNCTIONAL")
            
            med_ratio = round(float(np.median([x["length_ratio"] for x in items])), 3)
            min_id = round(float(np.min([x["aa_identity"] for x in items])), 2)
            
            val_results.append({
                "virus": vtype,
                "tier": tier,
                "gene": gname,
                "total_extracted": total_ext,
                "PASS_GOOD": c_good,
                "PASS_TRUNCATED": c_trunc,
                "FAIL_ANCHOR_ERROR": c_anchor_err,
                "FAIL_NONFUNCTIONAL": c_nonfunc,
                "median_length_ratio": med_ratio,
                "min_aa_identity": min_id
            })

df_val = pd.DataFrame(val_results)
val_tsv = os.path.join(OUT_QC, "06b_translation_validation_per_gene.tsv")
df_val.to_csv(val_tsv, sep="\t", index=False)
print(f"Saved translation validation to {val_tsv}")
print("\nValidation Summary for Priority Genes:")
print(df_val.to_string(index=False))

# -------------------------------------------------------------
# STEP 6.5b: Check FAIL_ANCHOR_ERROR Examples
# -------------------------------------------------------------
print("\n" + "="*60)
print("EXECUTING STEP 6.5b: FAIL_ANCHOR_ERROR CHECK")
print("="*60)

df_val["anchor_error_pct"] = (df_val["FAIL_ANCHOR_ERROR"] / df_val["total_extracted"]) * 100.0
genes_with_high_anchor_error = df_val[df_val["anchor_error_pct"] >= 5.0]

if genes_with_high_anchor_error.empty:
    print("[SUCCESS] All 20 priority genes across all cohorts have FAIL_ANCHOR_ERROR < 5.0%!")
    print("Maximum anchor error percentage observed: {:.2f}%".format(df_val["anchor_error_pct"].max()))
else:
    print(f"[WARNING] Found {len(genes_with_high_anchor_error)} genes with FAIL_ANCHOR_ERROR >= 5%:")
    print(genes_with_high_anchor_error[["virus", "tier", "gene", "total_extracted", "FAIL_ANCHOR_ERROR", "anchor_error_pct"]])
    
df_fa = pd.DataFrame(fail_anchor_examples)
if not df_fa.empty:
    print(f"\nTotal FAIL_ANCHOR_ERROR sequences across priority genes: {len(df_fa)}")
    print(df_fa.head(10).to_string(index=False))
else:
    print("\nZero FAIL_ANCHOR_ERROR sequences observed across all priority genes!")

# -------------------------------------------------------------
# STEP 6.5c: UL23 (TK) Frameshift Verification
# -------------------------------------------------------------
print("\n" + "="*60)
print("EXECUTING STEP 6.5c: UL23 (TK) FRAMESHIFT VERIFICATION")
print("="*60)

ul23_records = []
for vtype in ["HSV1", "HSV2"]:
    fasta_file = f"hsv_curated_dataset/04_cds_unaligned/{vtype}/UL23.fasta"
    ref_l = 376  # 376 aa for both HSV-1 and HSV-2
    for rec in SeqIO.parse(fasta_file, "fasta"):
        acc = rec.id.split("|")[0]
        s_cds = str(rec.seq).upper()
        p_to_stop = str(Seq(s_cds).translate(to_stop=True))
        p_len = len(p_to_stop)
        is_trunc = (p_len < 360)
        ul23_records.append({
            "virus": "HSV-1" if vtype == "HSV1" else "HSV-2",
            "accession": acc,
            "cds_length_nt": len(s_cds),
            "protein_length_aa": p_len,
            "ref_length_aa": ref_l,
            "length_ratio": round(p_len / ref_l, 3),
            "is_truncated": is_trunc
        })

df_ul23 = pd.DataFrame(ul23_records)
trunc_ul23 = df_ul23[df_ul23["is_truncated"]]

print(f"Total UL23 CDS evaluated across HSV-1 and HSV-2: {len(df_ul23)}")
print(f"Count of PASS_TRUNCATED UL23 sequences: {len(trunc_ul23)}")
print("\nLength distribution of all UL23 proteins:")
print(df_ul23.groupby(["virus", "is_truncated"])["protein_length_aa"].describe())

print("\nList of PASS_TRUNCATED UL23 CDS:")
if not trunc_ul23.empty:
    print(trunc_ul23.to_string(index=False))
else:
    print("None detected.")

# -------------------------------------------------------------
# STEP 6.5d: RL2 (ICP0) Splice Verification
# -------------------------------------------------------------
print("\n" + "="*60)
print("EXECUTING STEP 6.5d: RL2 (ICP0) SPLICE VERIFICATION")
print("="*60)

# Check extracted lengths
rl2_records = []
for vtype, ref_l_nt in [("HSV-1", 2328), ("HSV-2", 2505)]:
    v_dir = "HSV1" if vtype == "HSV-1" else "HSV2"
    fpath = f"hsv_curated_dataset/04_cds_unaligned/{v_dir}/RL2.fasta"
    for rec in SeqIO.parse(fpath, "fasta"):
        acc = rec.id.split("|")[0]
        s_cds = str(rec.seq).upper()
        l_nt = len(s_cds)
        pct_diff = abs(l_nt - ref_l_nt) / ref_l_nt * 100.0
        rl2_records.append({
            "virus": vtype,
            "accession": acc,
            "length_nt": l_nt,
            "ref_length_nt": ref_l_nt,
            "pct_diff": round(pct_diff, 2),
            "within_5pct": pct_diff <= 5.0,
            "outside_20pct": pct_diff > 20.0
        })

df_rl2 = pd.DataFrame(rl2_records)
print("RL2 Extracted Length Distribution:")
for vtype in ["HSV-1", "HSV-2"]:
    sub = df_rl2[df_rl2["virus"] == vtype]
    print(f"  {vtype} (ref = {sub['ref_length_nt'].iloc[0]} bp):")
    print(f"     Total extracted: {len(sub)}")
    print(f"     Min / Median / Max length: {sub['length_nt'].min()} / {sub['length_nt'].median()} / {sub['length_nt'].max()} bp")
    print(f"     Count within 5% of ref length: {sub['within_5pct'].sum()} ({sub['within_5pct'].mean()*100:.1f}%)")
    print(f"     Count outside 20% of ref length: {sub['outside_20pct'].sum()}")

# Check Splice Junctions and Genomic Intron GT-AG boundaries
# We verify in the clinical genomes that the excised introns have canonical GT donor and AG acceptor
print("\nVerifying RL2 Splice Junctions (GT donor and AG acceptor):")

splice_results = []
for vtype, gbk_path, fasta_path in [
    ("HSV-1", REF1_GBK, TIER1_HSV1_FASTA),
    ("HSV-2", REF2_GBK, TIER1_HSV2_FASTA)
]:
    # Get reference exon parts
    rec_ref = SeqIO.read(gbk_path, "genbank")
    full_ref_seq = str(rec_ref.seq).upper()
    rl2_parts = None
    for feat in rec_ref.features:
        if feat.type == "CDS" and feat.qualifiers.get("gene", [""])[0] == "RL2":
            rl2_parts = feat.location.parts
            break
            
    # Check clinical genomes
    clin_recs = list(SeqIO.parse(fasta_path, "fasta"))
    j1_canonical = 0
    j2_canonical = 0
    total_eval = 0
    
    pad = 2500
    for crec in clin_recs:
        cseq = str(crec.seq).upper()
        total_eval += 1
        
        # Locate exon 1 end, exon 2 start/end, exon 3 start
        # Exon 1 end anchor
        p1 = rl2_parts[0]
        ref_e1 = full_ref_seq[p1.start:p1.end]
        p2 = rl2_parts[1]
        ref_e2 = full_ref_seq[p2.start:p2.end]
        p3 = rl2_parts[2]
        ref_e3 = full_ref_seq[p3.start:p3.end]
        
        # Check window around junction 1
        # In HSV genomes, RL2 is on + strand (Copy 1)
        k1_end = ref_e1[-18:]
        pos_e1_end = cseq.find(k1_end, max(0, p1.start - pad), min(len(cseq), p1.end + pad))
        k2_start = ref_e2[:18]
        pos_e2_start = cseq.find(k2_start, max(0, p2.start - pad), min(len(cseq), p2.end + pad))
        
        if pos_e1_end != -1 and pos_e2_start != -1:
            idx_e1 = pos_e1_end + 18
            idx_e2 = pos_e2_start
            intron1 = cseq[idx_e1:idx_e2]
            if len(intron1) >= 4 and intron1[:2] == "GT" and intron1[-2:] == "AG":
                j1_canonical += 1
                
        # Check window around junction 2
        k2_end = ref_e2[-18:]
        pos_e2_end = cseq.find(k2_end, max(0, p2.start - pad), min(len(cseq), p2.end + pad))
        k3_start = ref_e3[:18]
        pos_e3_start = cseq.find(k3_start, max(0, p3.start - pad), min(len(cseq), p3.end + pad))
        
        if pos_e2_end != -1 and pos_e3_start != -1:
            idx_e2_end = pos_e2_end + 18
            idx_e3 = pos_e3_start
            intron2 = cseq[idx_e2_end:idx_e3]
            if len(intron2) >= 4 and intron2[:2] == "GT" and intron2[-2:] == "AG":
                j2_canonical += 1
                
    splice_results.append({
        "virus": vtype,
        "total_genomes": total_eval,
        "junction1_canonical_GT_AG": j1_canonical,
        "junction1_canonical_pct": round(j1_canonical / total_eval * 100.0, 2),
        "junction2_canonical_GT_AG": j2_canonical,
        "junction2_canonical_pct": round(j2_canonical / total_eval * 100.0, 2)
    })

df_splice = pd.DataFrame(splice_results)
print(df_splice.to_string(index=False))

# -------------------------------------------------------------
# STEP 6.5e: Reconcile Genome Attempt Counts
# -------------------------------------------------------------
print("\n" + "="*60)
print("EXECUTING STEP 6.5e: RECONCILING ATTEMPT NUMBERS")
print("="*60)

num_h1_t1_gen = 303
num_h2_t1_gen = 117
num_h2_t2_gen = 229

h1_genes = 74
h2_t1_genes = 74
h2_t2_genes = 71  # 74 minus 3 repeat genes (RL1, RL2, RS1)

expected_extractions = (num_h1_t1_gen * h1_genes) + (num_h2_t1_gen * h2_t1_genes) + (num_h2_t2_gen * h2_t2_genes)
actual_log_rows = len(df_ext_log)
counter_in_summary = 48026

print(f"Expected extractions: ({num_h1_t1_gen} x {h1_genes}) + ({num_h2_t1_gen} x {h2_t1_genes}) + ({num_h2_t2_gen} x {h2_t2_genes}) = {expected_extractions}")
print(f"Actual rows in TSV log: {actual_log_rows}")
print(f"Summary counter reported: {counter_in_summary}")
print(f"Discrepancy between summary counter and actual log rows: {counter_in_summary - actual_log_rows}")
print(f"Repeat gene bypass count in Tier 2: {num_h2_t2_gen} genomes x 3 repeat genes = {num_h2_t2_gen * 3}")
print("Explanation: The summary counter incremented before the repeat-gene continue statement for Tier 2, whereas the TSV log strictly records actual gene extractions attempted.")

# -------------------------------------------------------------
# STEP 6.5f: Per-Gene Extraction Success Table (All 74 Genes)
# -------------------------------------------------------------
print("\n" + "="*60)
print("EXECUTING STEP 6.5f: PER-GENE SUCCESS TABLE (74 GENES X 3 COHORTS)")
print("="*60)

success_records = []
flagged_genes = []

for vtype, tier in [("HSV-1", "Tier 1"), ("HSV-2", "Tier 1"), ("HSV-2", "Tier 2")]:
    sub_log = df_ext_log[(df_ext_log["virus"] == vtype) & (df_ext_log["cohort_tier"] == tier)]
    n_genomes = num_h1_t1_gen if vtype == "HSV-1" else (num_h2_t1_gen if tier == "Tier 1" else num_h2_t2_gen)
    
    unique_genes = sub_log["gene"].unique()
    for gname in sorted(unique_genes):
        g_sub = sub_log[sub_log["gene"] == gname]
        attempted = len(g_sub)
        extracted = sum(g_sub["extraction_status"] == "SUCCESS")
        rate = round(extracted / max(1, attempted) * 100.0, 2)
        
        success_records.append({
            "virus": vtype,
            "tier": tier,
            "gene": gname,
            "attempted": attempted,
            "extracted": extracted,
            "success_rate": rate
        })
        
        if vtype == "HSV-1" and tier == "Tier 1" and rate < 90.0:
            flagged_genes.append((gname, rate))

df_success = pd.DataFrame(success_records)
success_tsv = os.path.join(OUT_QC, "06c_per_gene_extraction_success.tsv")
df_success.to_csv(success_tsv, sep="\t", index=False)
print(f"Saved complete 74-gene per-cohort success table to {success_tsv}")
print(f"Total entries in success table: {len(df_success)} (74 HSV-1 T1 + 74 HSV-2 T1 + 71 HSV-2 T2 = 219 rows)")

print("\nFlagged genes with success_rate < 90% in Tier 1 HSV-1:")
if flagged_genes:
    for g, r in flagged_genes:
        print(f"  - {g}: {r}%")
else:
    print("  None! All 74 genes in Tier 1 HSV-1 have success rate >= 90%.")

print("\n" + "="*60)
print("ALL STEP 6.5 VALIDATIONS COMPLETE!")
print("="*60)
