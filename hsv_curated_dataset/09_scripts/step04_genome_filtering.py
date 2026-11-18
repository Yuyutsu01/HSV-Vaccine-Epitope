"""
Step 4: Genome-Level Filtering, Species Verification (ANI), and Exclusion Logging.
Filters all 1,305 complete genomes and produces filtered/flagged FASTA files and QC report.
"""

import os
import re
import math
import pandas as pd
from Bio import SeqIO
from collections import Counter, defaultdict

# Input paths
HSV1_RAW_FASTA = "hsv_curated_dataset/00_raw_downloads/HSV1/HSV1_all_genomes.fasta"
HSV2_RAW_FASTA = "hsv_curated_dataset/00_raw_downloads/HSV2/HSV2_all_genomes.fasta"
RAW_META_TSV = "hsv_curated_dataset/03_metadata/raw_metadata_all_genomes.tsv"
REF1_FASTA = "hsv_curated_dataset/01_reference_genomes/NC_001806.2.fasta"
REF2_FASTA = "hsv_curated_dataset/01_reference_genomes/NC_001798.2.fasta"

# Output paths
OUT_FILTERED_DIR = "hsv_curated_dataset/02_filtered_genomes"
OUT_QC_DIR = "hsv_curated_dataset/07_qc_reports"
os.makedirs(OUT_FILTERED_DIR, exist_ok=True)
os.makedirs(OUT_QC_DIR, exist_ok=True)

# -------------------------------------------------------------
# 1. Load Reference Genomes and Index K-mers for Fast ANI
# -------------------------------------------------------------
print("Loading reference genomes...")
ref1_record = SeqIO.read(REF1_FASTA, "fasta")
ref2_record = SeqIO.read(REF2_FASTA, "fasta")

ref1_seq = str(ref1_record.seq).upper()
ref2_seq = str(ref2_record.seq).upper()

K = 17
print(f"Indexing {K}-mers for reference genomes...")
kmers_ref1 = set(ref1_seq[i:i+K] for i in range(0, len(ref1_seq)-K+1, 2) if "N" not in ref1_seq[i:i+K])
kmers_ref2 = set(ref2_seq[i:i+K] for i in range(0, len(ref2_seq)-K+1, 2) if "N" not in ref2_seq[i:i+K])

def estimate_ani(query_seq, target_kmers):
    """
    Estimates Average Nucleotide Identity via k-mer containment and Jaccard transformation.
    Mash / containment distance: D = -1/k * ln(2*J / (1 + J)) or containment C = |A & B| / |A|.
    ANI ~ C^(1/k) * 100
    """
    q_kmers = set(query_seq[i:i+K] for i in range(0, len(query_seq)-K+1, 2) if "N" not in query_seq[i:i+K])
    if not q_kmers:
        return 0.0
    shared = len(q_kmers & target_kmers)
    containment = shared / len(q_kmers)
    if containment <= 0:
        return 0.0
    ani = (containment ** (1.0 / K)) * 100.0
    return min(100.0, ani)

# -------------------------------------------------------------
# 2. Screening Patterns
# -------------------------------------------------------------
ENGINEERED_PAT = re.compile(
    r"\b(BAC|oncolytic|vaccine vector|deletion mutant|recombinant construct|synthetic|engineered|patent|T-VEC|dl5-29|CJ9-gD|NV1020|vector)\b",
    re.I
)

LAB_STRAIN_PAT = re.compile(
    r"\b(strain 17\b|17\+|strain F\b|KOS|McKrae|333|HG52|strain MS\b|HFEM|MacIntyre|SC16|strain 186\b|strain RE\b)\b",
    re.I
)

def run_filtering():
    print("="*60)
    print("STARTING STEP 4: GENOME-LEVEL FILTERING")
    print("="*60)
    
    # Load raw metadata
    df_meta = pd.read_csv(RAW_META_TSV, sep="\t") if os.path.exists(RAW_META_TSV) else pd.DataFrame()
    meta_dict = {row["accession"]: row.to_dict() for _, row in df_meta.iterrows()} if not df_meta.empty else {}
    
    all_records = []
    
    # Read HSV-1 and HSV-2 raw genomes
    for vlabel, fasta_file in [("HSV-1", HSV1_RAW_FASTA), ("HSV-2", HSV2_RAW_FASTA)]:
        print(f"Reading {vlabel} genomes from {fasta_file}...")
        for r in SeqIO.parse(fasta_file, "fasta"):
            all_records.append((vlabel, r))
    print(f"Loaded {len(all_records)} total raw genomes.")
    
    # Tracking sets for duplicate detection
    seen_sequences = {}
    seen_isolates = defaultdict(list)
    
    filtering_results = []
    
    hsv1_filtered_records = []
    hsv2_filtered_records = []
    hsv1_flagged_records = []
    hsv2_flagged_records = []
    
    for idx, (db_type, record) in enumerate(all_records):
        acc = record.id.split(".")[0]
        full_id = record.id
        desc = record.description
        seq_str = str(record.seq).upper()
        seq_len = len(seq_str)
        
        # Calculate N count & percentage
        n_count = seq_str.count("N")
        n_pct = (n_count / max(1, seq_len)) * 100.0
        
        # Extract isolate/strain from description
        m_iso = re.search(r"\b(?:strain|isolate)\s+([A-Za-z0-9_\-\.\/]+)", desc, re.I)
        isolate = m_iso.group(1) if m_iso else "Unspecified"
        
        # Check engineered / synthetic
        meta_info = meta_dict.get(acc, {})
        meta_title = str(meta_info.get("title", "")) if pd.notna(meta_info.get("title", "")) else ""
        meta_subname = str(meta_info.get("subname", "")) if pd.notna(meta_info.get("subname", "")) else ""
        title_to_check = f"{desc} {meta_title} {meta_subname}"
        is_engineered = bool(ENGINEERED_PAT.search(title_to_check))
        
        # Check lab strain
        is_lab_strain = bool(LAB_STRAIN_PAT.search(title_to_check))
        
        # Compute ANI against both reference genomes
        ani_hsv1 = round(estimate_ani(seq_str, kmers_ref1), 2)
        ani_hsv2 = round(estimate_ani(seq_str, kmers_ref2), 2)
        
        # Computational species assignment
        if ani_hsv1 >= 95.0 and ani_hsv2 < 90.0:
            comp_type = "HSV-1"
        elif ani_hsv2 >= 95.0 and ani_hsv1 < 90.0:
            comp_type = "HSV-2"
        elif ani_hsv1 >= 85.0 and ani_hsv2 >= 85.0:
            comp_type = "RECOMBINANT_HSV1_HSV2"
        elif ani_hsv1 < 80.0 and ani_hsv2 < 80.0:
            comp_type = "NON_HSV"
        elif ani_hsv1 > ani_hsv2:
            comp_type = "HSV-1"
        else:
            comp_type = "HSV-2"
            
        is_recombinant = (comp_type == "RECOMBINANT_HSV1_HSV2")
        is_mislabeled = (comp_type in ["HSV-1", "HSV-2"] and comp_type != db_type)
        
        # Check duplicates
        # Sequence hash / prefix-suffix hash for speed
        seq_hash = (seq_len, seq_str[:500], seq_str[-500:], seq_str[seq_len//2:seq_len//2+500])
        is_exact_dup = seq_hash in seen_sequences
        dup_of = seen_sequences.get(seq_hash, None)
        
        # Decision logic
        status = "PASS"
        action = "FILTERED_KEEP"
        reasons = []
        
        # Exclusion criteria
        if is_engineered:
            status = "EXCLUDED"
            action = "EXCLUDED"
            reasons.append("Engineered, synthetic, oncolytic, or patent construct")
        elif seq_len < 140000 or seq_len > 160000:
            status = "EXCLUDED"
            action = "EXCLUDED"
            reasons.append(f"Length outside 140-160 kb ({seq_len} bp)")
        elif n_pct > 1.0:
            status = "EXCLUDED"
            action = "EXCLUDED"
            reasons.append(f"Ambiguity > 1% N ({n_pct:.2f}% N)")
        elif is_exact_dup:
            status = "EXCLUDED"
            action = "EXCLUDED"
            reasons.append(f"Exact sequence duplicate of {dup_of}")
        elif comp_type == "NON_HSV":
            status = "EXCLUDED"
            action = "EXCLUDED"
            reasons.append("Non-HSV contaminant based on ANI")
            
        # Flag but keep criteria (only if not already excluded)
        if status != "EXCLUDED":
            flag_reasons = []
            if is_lab_strain:
                flag_reasons.append(f"Laboratory passaged strain")
            if is_recombinant:
                flag_reasons.append(f"Recombinant HSV-1/HSV-2 genome")
            if is_mislabeled:
                flag_reasons.append(f"Mislabeled species (Database={db_type}, ANI={comp_type})")
            if isolate != "Unspecified" and isolate in seen_isolates[db_type]:
                flag_reasons.append(f"Multiple submission of isolate {isolate}")
                
            if flag_reasons:
                status = "FLAGGED"
                action = "FLAGGED_KEEP"
                reasons.extend(flag_reasons)
            else:
                status = "PASS"
                action = "FILTERED_KEEP"
                reasons.append("Passed all genome-level quality filters")
                seen_sequences[seq_hash] = acc
                if isolate != "Unspecified":
                    seen_isolates[db_type].append(isolate)

        # Store to output bins
        assigned_target_virus = comp_type if comp_type in ["HSV-1", "HSV-2"] else db_type
        
        if action == "FILTERED_KEEP":
            if assigned_target_virus == "HSV-1":
                hsv1_filtered_records.append(record)
            else:
                hsv2_filtered_records.append(record)
        elif action == "FLAGGED_KEEP":
            if assigned_target_virus == "HSV-1":
                hsv1_flagged_records.append(record)
            else:
                hsv2_flagged_records.append(record)

        filtering_results.append({
            "accession": acc,
            "full_header": desc,
            "database_type": db_type,
            "computational_type": comp_type,
            "ani_to_hsv1_ref": ani_hsv1,
            "ani_to_hsv2_ref": ani_hsv2,
            "genome_length": seq_len,
            "n_count": n_count,
            "n_pct": round(n_pct, 2),
            "isolate": isolate,
            "is_lab_strain": is_lab_strain,
            "is_engineered": is_engineered,
            "is_recombinant": is_recombinant,
            "is_mislabeled": is_mislabeled,
            "is_exact_dup": is_exact_dup,
            "filter_status": status,
            "filter_action": action,
            "decision_reasons": "; ".join(reasons)
        })

    # Save output FASTA files
    SeqIO.write(hsv1_filtered_records, os.path.join(OUT_FILTERED_DIR, "HSV1_filtered.fasta"), "fasta")
    SeqIO.write(hsv2_filtered_records, os.path.join(OUT_FILTERED_DIR, "HSV2_filtered.fasta"), "fasta")
    SeqIO.write(hsv1_flagged_records, os.path.join(OUT_FILTERED_DIR, "HSV1_flagged.fasta"), "fasta")
    SeqIO.write(hsv2_flagged_records, os.path.join(OUT_FILTERED_DIR, "HSV2_flagged.fasta"), "fasta")
    
    # Save exclusion and filtering log
    df_filter = pd.DataFrame(filtering_results)
    log_path = os.path.join(OUT_QC_DIR, "04_genome_filtering_log.tsv")
    df_filter.to_csv(log_path, sep="\t", index=False)
    print(f"Saved filtering log to {log_path} with {len(df_filter)} entries.")
    
    # -------------------------------------------------------------
    # Summary Reporting
    # -------------------------------------------------------------
    print("\n" + "="*60)
    print("STEP 4 FILTERING SUMMARY REPORT")
    print("="*60)
    
    for vlabel in ["HSV-1", "HSV-2"]:
        sub = df_filter[df_filter["database_type"] == vlabel]
        tot = len(sub)
        n_pass = len(hsv1_filtered_records if vlabel == "HSV-1" else hsv2_filtered_records)
        n_flag = len(hsv1_flagged_records if vlabel == "HSV-1" else hsv2_flagged_records)
        n_excl = len(sub[sub["filter_status"] == "EXCLUDED"])
        
        print(f"\n{vlabel} GENOMES:")
        print(f"  Total Downloaded Genomes         : {tot}")
        print(f"  Passed High-Quality Filter (KEEP): {n_pass} ({n_pass/tot*100:.1f}%)")
        print(f"  Flagged (Lab/Recomb/Patient)     : {n_flag} ({n_flag/tot*100:.1f}%)")
        print(f"  Excluded (Failed Length/N%/Eng)  : {n_excl} ({n_excl/tot*100:.1f}%)")
        print(f"  Total Retained (Clean + Flagged) : {n_pass + n_flag}")

    # Exclusions breakdown
    print("\nEXCLUSION REASONS BREAKDOWN:")
    excl_df = df_filter[df_filter["filter_status"] == "EXCLUDED"]
    excl_reasons = Counter()
    for r in excl_df["decision_reasons"]:
        for item in r.split("; "):
            excl_reasons[item] += 1
    for reason, count in excl_reasons.most_common():
        print(f"  {count:4d} : {reason}")

if __name__ == "__main__":
    run_filtering()
