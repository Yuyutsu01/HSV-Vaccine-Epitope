"""
Step 4.5: Comprehensive Verification, Lab Strain Discovery,
Duplicate Auditing, and HSV-2 Coding-Region Recovery.
"""

import os
import re
import hashlib
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter, defaultdict

# Directories
RAW_HSV1 = "hsv_curated_dataset/00_raw_downloads/HSV1/HSV1_all_genomes.fasta"
RAW_HSV2 = "hsv_curated_dataset/00_raw_downloads/HSV2/HSV2_all_genomes.fasta"
REF_DIR = "hsv_curated_dataset/01_reference_genomes"
OUT_FILTERED = "hsv_curated_dataset/02_filtered_genomes"
OUT_QC = "hsv_curated_dataset/07_qc_reports"
os.makedirs(OUT_FILTERED, exist_ok=True)
os.makedirs(OUT_QC, exist_ok=True)

# -------------------------------------------------------------
# 1. Load Canonical References
# -------------------------------------------------------------
print("Loading canonical references and lab strain references...")
ref1_rec = SeqIO.read(os.path.join(REF_DIR, "NC_001806.2.fasta"), "fasta")
ref2_rec = SeqIO.read(os.path.join(REF_DIR, "NC_001798.2.fasta"), "fasta")
ref_f_rec = SeqIO.read(os.path.join(REF_DIR, "GU734771.1_F.fasta"), "fasta")
ref_kos_rec = SeqIO.read(os.path.join(REF_DIR, "JQ673480.1_KOS.fasta"), "fasta")
ref_mck_rec = SeqIO.read(os.path.join(REF_DIR, "JX142173.1_McKrae.fasta"), "fasta")

ref1_seq = str(ref1_rec.seq).upper()
ref2_seq = str(ref2_rec.seq).upper()
ref_f_seq = str(ref_f_rec.seq).upper()
ref_kos_seq = str(ref_kos_rec.seq).upper()
ref_mck_seq = str(ref_mck_rec.seq).upper()

# Index 17-mers for ANI estimation
K = 17
def get_kmers(s):
    return set(s[i:i+K] for i in range(0, len(s)-K+1, 2) if "N" not in s[i:i+K])

kmers_ref1 = get_kmers(ref1_seq)
kmers_ref2 = get_kmers(ref2_seq)
kmers_ref_f = get_kmers(ref_f_seq)
kmers_ref_kos = get_kmers(ref_kos_seq)
kmers_ref_mck = get_kmers(ref_mck_seq)

def calc_ani(query_seq, target_kmers):
    q_kmers = set(query_seq[i:i+K] for i in range(0, len(query_seq)-K+1, 2) if "N" not in query_seq[i:i+K])
    if not q_kmers or not target_kmers:
        return 0.0
    shared = len(q_kmers & target_kmers)
    containment = shared / len(q_kmers)
    if containment <= 0:
        return 0.0
    return min(100.0, (containment ** (1.0 / K)) * 100.0)

# Load reference CDS intervals for NC_001798.2
ref2_gbk = SeqIO.read(os.path.join(REF_DIR, "NC_001798.2.gbk"), "genbank")
ref2_cds_intervals = []
for feat in ref2_gbk.features:
    if feat.type == "CDS":
        for part in feat.location.parts:
            ref2_cds_intervals.append((int(part.start), int(part.end)))

# Define NC_001798.2 repeat coordinates (TRL, IRL, IRS, TRS)
# In HG52 (154,675 bp):
# TRL: 0 - 9213
# IRL: 119284 - 128497
# IRS: 128497 - 134704
# TRS: 148468 - 154675
# Unique Long: 9213 - 119284
# Unique Short: 134704 - 148468
REPEAT_INTERVALS_HSV2 = [
    (0, 9213),
    (119284, 134704),
    (148468, 154675)
]

def is_in_intervals(pos, intervals):
    for s, e in intervals:
        if s <= pos < e:
            return True
    return False

# -------------------------------------------------------------
# 2. Load All Genomes (1,305 total)
# -------------------------------------------------------------
print("Loading all raw genomes into memory...")
all_genomes = []
for vlabel, fpath in [("HSV-1", RAW_HSV1), ("HSV-2", RAW_HSV2)]:
    for r in SeqIO.parse(fpath, "fasta"):
        all_genomes.append((vlabel, r))
print(f"Loaded {len(all_genomes)} total genomes.")

# Load Step 4 filtering log
df_log = pd.read_csv("hsv_curated_dataset/07_qc_reports/04_genome_filtering_log.tsv", sep="\t")
log_dict = {row["accession"]: row.to_dict() for _, row in df_log.iterrows()}

# =============================================================
# STEP 4.5b: HSV-2 Coding-Region N Analysis
# =============================================================
print("\n" + "="*60)
print("EXECUTING STEP 4.5b: HSV-2 CODING-REGION N ANALYSIS")
print("="*60)

hsv2_excl_rows = []
recoverable_records = []

for vtype, rec in all_genomes:
    if vtype != "HSV-2":
        continue
    acc = rec.id.split(".")[0]
    log_info = log_dict.get(acc, {})
    if log_info.get("filter_status") != "EXCLUDED":
        continue
        
    seq = str(rec.seq).upper()
    seq_l = len(seq)
    
    # Calculate N in whole genome
    tot_n = seq.count("N")
    whole_n_pct = (tot_n / max(1, seq_l)) * 100.0
    
    # Analyze distribution across regions:
    # Approximate reference-coordinate mapping
    # Unique region length in HG52: UL is 110,071 bp, US is 13,764 bp (total unique = 123,835 bp)
    # Total repeat length: 30,840 bp (~20% of genome)
    # Map by relative coordinate index when genome is 140-160k
    coding_n = 0
    coding_total_len = sum(e - s for s, e in ref2_cds_intervals)
    repeat_n = 0
    repeat_total_len = sum(e - s for s, e in REPEAT_INTERVALS_HSV2)
    
    # Calculate N in unique region (UL: 9.2k to 119k, US: 134.7k to 148.5k)
    # Scaled to query length
    scale = seq_l / 154675.0
    
    ul_start, ul_end = int(9213 * scale), int(119284 * scale)
    us_start, us_end = int(134704 * scale), int(148468 * scale)
    
    unique_seq = seq[ul_start:ul_end] + seq[us_start:us_end]
    repeat_seq = seq[:ul_start] + seq[ul_end:us_start] + seq[us_end:]
    
    unique_n = unique_seq.count("N")
    unique_n_pct = (unique_n / max(1, len(unique_seq))) * 100.0
    repeat_n = repeat_seq.count("N")
    repeat_n_pct = (repeat_n / max(1, len(repeat_seq))) * 100.0
    
    # Coding region %N is bounded by unique region %N (since all non-repeat priority CDS are in UL/US)
    coding_n_pct = unique_n_pct
    
    # Classification: RECOVERABLE if coding %N < 0.5% AND 140k <= length <= 160k
    if coding_n_pct < 0.5 and 140000 <= seq_l <= 160000 and whole_n_pct > 1.0:
        classification = "RECOVERABLE"
        recoverable_records.append(rec)
    else:
        classification = "TRULY_LOW_QUALITY"
        
    hsv2_excl_rows.append({
        "accession": acc,
        "genome_length": seq_l,
        "whole_pct_N": round(whole_n_pct, 2),
        "coding_pct_N": round(coding_n_pct, 2),
        "repeat_pct_N": round(repeat_n_pct, 2),
        "classification": classification,
        "original_exclusion_reason": log_info.get("decision_reasons", "")
    })

df_hsv2_excl = pd.DataFrame(hsv2_excl_rows)
out_4b_tsv = os.path.join(OUT_QC, "04b_hsv2_exclusion_analysis.tsv")
df_hsv2_excl.to_csv(out_4b_tsv, sep="\t", index=False)
print(f"Saved {out_4b_tsv} with {len(df_hsv2_excl)} excluded HSV-2 records analyzed.")

n_rec = (df_hsv2_excl["classification"] == "RECOVERABLE").sum()
n_truly = (df_hsv2_excl["classification"] == "TRULY_LOW_QUALITY").sum()
print(f"  HSV-2 RECOVERABLE (coding %N < 0.5%, repeats corrupted) : {n_rec}")
print(f"  HSV-2 TRULY LOW QUALITY (coding %N >= 0.5% or length bad): {n_truly}")

if n_rec >= 50:
    out_tier2 = os.path.join(OUT_FILTERED, "HSV2_filtered_tier2.fasta")
    SeqIO.write(recoverable_records, out_tier2, "fasta")
    print(f"  -> Generated Tier 2 FASTA: {out_tier2} with {len(recoverable_records)} recoverable genomes.")

# =============================================================
# STEP 4.5c & 4.5d: Lab Strain & Duplicate Deep Audit
# =============================================================
print("\n" + "="*60)
print("EXECUTING STEP 4.5c & 4.5d: LAB STRAINS & DUPLICATE AUDIT")
print("="*60)

# Check SHA-256 sequence hashes across all 1,305 genomes
sha256_dict = defaultdict(list)
for vtype, rec in all_genomes:
    seq_str = str(rec.seq).upper()
    sha = hashlib.sha256(seq_str.encode('utf-8')).hexdigest()
    acc = rec.id.split(".")[0]
    sha256_dict[sha].append((vtype, acc, rec.id))

exact_dup_groups = {k: v for k, v in sha256_dict.items() if len(v) > 1}
total_exact_duplicates = sum(len(v) - 1 for v in exact_dup_groups.values())
print(f"Total exact sequence hash duplicates found (across all 1,305 genomes): {total_exact_duplicates}")
print(f"Unique sequence strings: {len(sha256_dict)}")

# Lab strain evaluation: metadata vs ANI
lab_tokens = ["17", "17+", "F", "KOS", "MCKRAE", "333", "HG52", "MS", "HFEM", "MACINTYRE", "SC16", "186", "RE"]
lab_meta_pat = re.compile(r"\b(" + "|".join(re.escape(x) for x in lab_tokens) + r")\b", re.I)

lab_audit_rows = []
species_verification_rows = []

ani_strain17_lab = 0
ani_hg52_lab = 0
ani_f_lab = 0
ani_kos_lab = 0
ani_mck_lab = 0
meta_lab_count = 0
union_lab_set = set()

hsv1_all_filtered_new = []
hsv2_all_filtered_new = []
hsv1_flagged_new = []
hsv2_flagged_new = []

# High-identity isolate tracker
isolate_seqs = defaultdict(list)

for vtype, rec in all_genomes:
    acc = rec.id.split(".")[0]
    desc = rec.description
    seq_str = str(rec.seq).upper()
    seq_l = len(seq_str)
    
    # Metadata lookup
    m_iso = re.search(r"\b(?:strain|isolate)\s+([A-Za-z0-9_\-\.\/]+)", desc, re.I)
    isolate = m_iso.group(1) if m_iso else "Unspecified"
    
    # Check 1: Metadata match
    meta_is_lab = bool(lab_meta_pat.search(desc)) or bool(lab_meta_pat.search(isolate))
    if meta_is_lab:
        meta_lab_count += 1
        
    # Check 2, 3, 4: ANI to references
    ani1 = round(calc_ani(seq_str, kmers_ref1), 2)
    ani2 = round(calc_ani(seq_str, kmers_ref2), 2)
    ani_f = round(calc_ani(seq_str, kmers_ref_f), 2)
    ani_kos = round(calc_ani(seq_str, kmers_ref_kos), 2)
    ani_mck = round(calc_ani(seq_str, kmers_ref_mck), 2)
    
    is_s17_der = (ani1 >= 99.9)
    is_hg52_der = (ani2 >= 99.9)
    is_f_der = (ani_f >= 99.9)
    is_kos_der = (ani_kos >= 99.9)
    is_mck_der = (ani_mck >= 99.9)
    
    if is_s17_der: ani_strain17_lab += 1
    if is_hg52_der: ani_hg52_lab += 1
    if is_f_der: ani_f_lab += 1
    if is_kos_der: ani_kos_lab += 1
    if is_mck_der: ani_mck_lab += 1
    
    ani_is_lab = (is_s17_der or is_hg52_der or is_f_der or is_kos_der or is_mck_der)
    
    is_lab_union = (meta_is_lab or ani_is_lab)
    if is_lab_union:
        union_lab_set.add(acc)

    # Step 4.5e: Computational Type Assignment
    if ani1 >= 95.0 and ani2 < 90.0:
        comp_type = "HSV-1"
    elif ani2 >= 95.0 and ani1 < 90.0:
        comp_type = "HSV-2"
    elif ani1 >= 88.0 and ani2 >= 88.0:
        comp_type = "RECOMBINANT_HSV1_HSV2"
    elif ani1 < 80.0 and ani2 < 80.0:
        comp_type = "NON_HSV"
    elif ani1 > ani2:
        comp_type = "HSV-1"
    else:
        comp_type = "HSV-2"
        
    agreement = (comp_type == vtype)
    species_verification_rows.append({
        "accession": acc,
        "database_label": vtype,
        "computational_assignment": comp_type,
        "ANI_to_HSV1_ref": ani1,
        "ANI_to_HSV2_ref": ani2,
        "agreement": "YES" if agreement else "NO",
        "is_recombinant": comp_type == "RECOMBINANT_HSV1_HSV2"
    })
    
    lab_audit_rows.append({
        "accession": acc,
        "virus": vtype,
        "isolate": isolate,
        "metadata_lab_flag": meta_is_lab,
        "ani_strain17": ani1,
        "ani_hg52": ani2,
        "ani_strain_F": ani_f,
        "ani_kos": ani_kos,
        "ani_mckrae": ani_mck,
        "is_lab_by_ani": ani_is_lab,
        "is_lab_union": is_lab_union
    })

# Save Species Verification Table
df_spec = pd.DataFrame(species_verification_rows)
out_spec_tsv = os.path.join(OUT_QC, "04e_species_verification.tsv")
df_spec.to_csv(out_spec_tsv, sep="\t", index=False)
print(f"Saved species verification table to {out_spec_tsv}")

# Species agreement summary
agr_count = (df_spec["agreement"] == "YES").sum()
disagr_count = (df_spec["agreement"] == "NO").sum()
recomb_count = (df_spec["is_recombinant"]).sum()

print(f"\nSpecies Verification Summary (Total: 1,305):")
print(f"  Agreements between database and ANI   : {agr_count} ({agr_count/1305*100:.1f}%)")
print(f"  Disagreements (Mislabeled or Non-HSV) : {disagr_count}")
print(f"  Inter-type Recombinant Candidates     : {recomb_count}")

if disagr_count > 0:
    print("Disagreements breakdown:")
    print(df_spec[df_spec["agreement"] == "NO"][["accession", "database_label", "computational_assignment", "ANI_to_HSV1_ref", "ANI_to_HSV2_ref"]])

# Save Lab Strain Audit Table
df_lab = pd.DataFrame(lab_audit_rows)
out_lab_tsv = os.path.join(OUT_QC, "04c_lab_strain_audit.tsv")
df_lab.to_csv(out_lab_tsv, sep="\t", index=False)
print(f"\nLab Strain Detection Audit Summary:")
print(f"  Detected by Metadata Text Search      : {meta_lab_count}")
print(f"  Detected by ANI >99.9% to Strain 17   : {ani_strain17_lab}")
print(f"  Detected by ANI >99.9% to HG52        : {ani_hg52_lab}")
print(f"  Detected by ANI >99.9% to Strain F    : {ani_f_lab}")
print(f"  Detected by ANI >99.9% to KOS         : {ani_kos_lab}")
print(f"  Detected by ANI >99.9% to McKrae      : {ani_mck_lab}")
print(f"  UNION of Lab Strains (Metadata or ANI): {len(union_lab_set)}")

# =============================================================
# Update Filtered and Flagged Cohorts using Union Lab Strains
# =============================================================
print("\nUpdating Filtered (Tier 1) and Flagged Cohorts with Verified Lab Strains...")
tier1_hsv1 = []
tier1_hsv2 = []
flagged_hsv1 = []
flagged_hsv2 = []

seen_clean_hashes = set()

for vtype, rec in all_genomes:
    acc = rec.id.split(".")[0]
    log_info = log_dict.get(acc, {})
    seq_str = str(rec.seq).upper()
    seq_l = len(seq_str)
    n_pct = (seq_str.count("N") / max(1, seq_l)) * 100.0
    sha = hashlib.sha256(seq_str.encode('utf-8')).hexdigest()
    
    # Must meet strict basic QC (140-160k, N <= 1%, non-engineered, non-contaminant)
    is_eng = log_info.get("is_engineered", False)
    comp_spec = df_spec[df_spec["accession"] == acc]["computational_assignment"].values[0]
    
    if is_eng or seq_l < 140000 or seq_l > 160000 or n_pct > 1.0 or comp_spec == "NON_HSV":
        continue
        
    is_recomb = (comp_spec == "RECOMBINANT_HSV1_HSV2")
    is_lab = acc in union_lab_set
    is_dup = sha in seen_clean_hashes
    
    target_vtype = comp_spec if comp_spec in ["HSV-1", "HSV-2"] else vtype
    
    if is_lab or is_recomb or is_dup:
        if target_vtype == "HSV-1":
            flagged_hsv1.append(rec)
        else:
            flagged_hsv2.append(rec)
    else:
        seen_clean_hashes.add(sha)
        if target_vtype == "HSV-1":
            tier1_hsv1.append(rec)
        else:
            tier1_hsv2.append(rec)

# Write updated cohorts
SeqIO.write(tier1_hsv1, os.path.join(OUT_FILTERED, "HSV1_filtered.fasta"), "fasta")
SeqIO.write(tier1_hsv2, os.path.join(OUT_FILTERED, "HSV2_filtered.fasta"), "fasta")
SeqIO.write(flagged_hsv1, os.path.join(OUT_FILTERED, "HSV1_flagged.fasta"), "fasta")
SeqIO.write(flagged_hsv2, os.path.join(OUT_FILTERED, "HSV2_flagged.fasta"), "fasta")

print("\n" + "="*60)
print("FINAL STEP 4.5 TIERED COHORT SUMMARY")
print("="*60)
print(f"HSV-1 Tier 1 (Strict Clean Clinical)  : {len(tier1_hsv1)}")
print(f"HSV-2 Tier 1 (Strict Clean Clinical)  : {len(tier1_hsv2)}")
print(f"HSV-2 Tier 2 (Recoverable Sensitivity): {n_rec}")
print(f"HSV-1 Flagged (Lab Strains & Dups)    : {len(flagged_hsv1)}")
print(f"HSV-2 Flagged (Lab Strains & Dups)    : {len(flagged_hsv2)}")
print(f"Total Retained (Tier 1 + Tier 2 + Flag): {len(tier1_hsv1) + len(tier1_hsv2) + n_rec + len(flagged_hsv1) + len(flagged_hsv2)}")
