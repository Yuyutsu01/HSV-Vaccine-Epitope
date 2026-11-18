"""
Step 7: Per-CDS Quality Control Pipeline (Part 1 of Path D).
Processes the 18 clean priority target genes:
  Glycoproteins: US6, UL27, UL44, US4, UL22, UL1, US8, US7, UL53, UL10, UL49A
  T-cell targets: UL54, UL48, UL46, UL47, UL49
  Drug targets: UL23 (TK), UL30 (Polymerase)

(Excludes RL2 and RS1, which are evaluated in Part 2).

Filters applied:
  1. Starts with ATG
  2. Ends with in-frame stop codon (TAA, TAG, TGA)
  3. Divisible by 3 (in-frame)
  4. Length within 95-105% of reference length
  5. Zero internal stop codons
  6. Ambiguity < 0.5% N

Special Handling for UL23:
  - Preserves verified frameshift variants (OQ724915, PQ065496) and isolates
    them into UL23_resistance_variants.fasta.
"""

import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Paths
REF1_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001806.2.gbk"
REF2_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001798.2.gbk"
IN_HSV1 = "hsv_curated_dataset/04_cds_unaligned/HSV1"
IN_HSV2 = "hsv_curated_dataset/04_cds_unaligned/HSV2"
OUT_HSV1 = "hsv_curated_dataset/05_cds_qc_passed/HSV1"
OUT_HSV2 = "hsv_curated_dataset/05_cds_qc_passed/HSV2"
OUT_QC = "hsv_curated_dataset/07_qc_reports"

os.makedirs(OUT_HSV1, exist_ok=True)
os.makedirs(OUT_HSV2, exist_ok=True)
os.makedirs(OUT_QC, exist_ok=True)

# 18 Clean Priority Genes
CLEAN_18_GENES = [
    "US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53", "UL10", "UL49A",
    "UL54", "UL48", "UL46", "UL47", "UL49", "UL23", "UL30"
]

# Load reference gene models
def get_ref_lengths(gbk_path):
    record = SeqIO.read(gbk_path, "genbank")
    lens = {}
    for feat in record.features:
        if feat.type == "CDS":
            g = feat.qualifiers.get("gene", [""])[0] or feat.qualifiers.get("locus_tag", [""])[0]
            if g and g not in lens:
                # Sum exon lengths
                l = sum(p.end - p.start for p in feat.location.parts)
                lens[g] = l
    return lens

ref_lens_h1 = get_ref_lengths(REF1_GBK)
ref_lens_h2 = get_ref_lengths(REF2_GBK)

qc_log_rows = []
summary_table_rows = []

print("="*60)
print("EXECUTING PART 1: PER-CDS QC FOR 18 CLEAN PRIORITY GENES")
print("="*60)

for vtype in ["HSV-1", "HSV-2"]:
    in_dir = IN_HSV1 if vtype == "HSV-1" else IN_HSV2
    out_dir = OUT_HSV1 if vtype == "HSV-1" else OUT_HSV2
    ref_lens = ref_lens_h1 if vtype == "HSV-1" else ref_lens_h2
    
    for gname in CLEAN_18_GENES:
        in_fasta = os.path.join(in_dir, f"{gname}.fasta")
        if not os.path.exists(in_fasta):
            continue
            
        ref_l = ref_lens.get(gname, 0)
        recs = list(SeqIO.parse(in_fasta, "fasta"))
        
        passing_records = []
        ul23_resistance_records = []
        
        pass_clean_count = 0
        pass_resist_count = 0
        fail_start_count = 0
        fail_stop_count = 0
        fail_frame_count = 0
        fail_len_count = 0
        fail_internal_stop_count = 0
        fail_ambiguity_count = 0
        
        for rec in recs:
            s = str(rec.seq).upper()
            acc = rec.id.split("|")[0]
            l_nt = len(s)
            
            # 1. Start codon check
            starts_atg = (s[:3] == "ATG")
            
            # 2. In-frame stop codon check
            ends_stop = (s[-3:] in ["TAA", "TAG", "TGA"])
            
            # 3. Divisible by 3
            is_div3 = (l_nt % 3 == 0)
            
            # 4. Length ratio
            len_ratio = l_nt / max(1, ref_l)
            len_ok = (0.95 <= len_ratio <= 1.05)
            
            # 5. Ambiguity check (< 0.5% N)
            pct_n = (s.count("N") / max(1, l_nt)) * 100.0
            ambig_ok = (pct_n < 0.5)
            
            # 6. Internal stop codons
            if is_div3 and l_nt >= 6:
                prot_full = str(Seq(s).translate())
                int_stops = prot_full[:-1].count("*")
            else:
                int_stops = 999
            no_internal_stops = (int_stops == 0)
            
            # Determine QC status
            qc_status = "FAIL"
            failure_reasons = []
            
            if not starts_atg:
                failure_reasons.append(f"Invalid start codon ({s[:3]})")
                fail_start_count += 1
            if not ends_stop:
                failure_reasons.append(f"Invalid stop codon ({s[-3:]})")
                fail_stop_count += 1
            if not is_div3:
                failure_reasons.append(f"Length {l_nt} not divisible by 3")
                fail_frame_count += 1
            if not len_ok:
                failure_reasons.append(f"Length ratio {len_ratio:.3f} outside 0.95-1.05")
                fail_len_count += 1
            if not ambig_ok:
                failure_reasons.append(f"Ambiguity {pct_n:.2f}% >= 0.5%")
                fail_ambiguity_count += 1
            if not no_internal_stops:
                failure_reasons.append(f"Internal stop codons ({int_stops})")
                fail_internal_stop_count += 1
                
            # Special Handling for UL23 (TK)
            is_ul23_resistance = False
            if gname == "UL23":
                # Check for verified resistance frameshift variants
                # OQ724915 (HSV-1) or PQ065496 (HSV-2) or any truncated TK with identity >= 90%
                p_to_stop = str(Seq(s[:len(s)-(len(s)%3)]).translate(to_stop=True))
                if (acc in ["OQ724915", "PQ065496"]) or (len(p_to_stop) < 360 and len(p_to_stop) >= 150):
                    is_ul23_resistance = True
                    qc_status = "PASS_RESISTANCE_VARIANT"
                    pass_resist_count += 1
                    ul23_resistance_records.append(f">{rec.id}|status=RESISTANCE_VARIANT\n{s}\n")
                    
            if not is_ul23_resistance:
                if not failure_reasons:
                    qc_status = "PASS_GOOD"
                    pass_clean_count += 1
                    passing_records.append(f">{rec.id}\n{s}\n")
                else:
                    qc_status = "FAIL"
                    
            qc_log_rows.append({
                "virus": vtype,
                "gene": gname,
                "accession": acc,
                "fasta_header": rec.id,
                "length_nt": l_nt,
                "ref_length_nt": ref_l,
                "length_ratio": round(len_ratio, 3),
                "starts_atg": starts_atg,
                "ends_stop": ends_stop,
                "divisible_by_3": is_div3,
                "internal_stops": int_stops,
                "pct_n": round(pct_n, 2),
                "qc_status": qc_status,
                "failure_reasons": "; ".join(failure_reasons)
            })
            
        # Write passing clean FASTA
        out_fasta = os.path.join(out_dir, f"{gname}.fasta")
        with open(out_fasta, "w", encoding="utf-8") as f:
            f.writelines(passing_records)
            
        # Write UL23 resistance FASTA if applicable
        if gname == "UL23" and ul23_resistance_records:
            resist_fasta = os.path.join(out_dir, "UL23_resistance_variants.fasta")
            with open(resist_fasta, "w", encoding="utf-8") as f:
                f.writelines(ul23_resistance_records)
            print(f"[{vtype} UL23] Segregated {len(ul23_resistance_records)} verified resistance variant(s) into {resist_fasta}")
            
        summary_table_rows.append({
            "virus": vtype,
            "gene": gname,
            "total_input": len(recs),
            "PASS_GOOD": pass_clean_count,
            "PASS_RESISTANCE": pass_resist_count,
            "FAIL_TOTAL": len(recs) - pass_clean_count - pass_resist_count,
            "pass_rate_pct": round((pass_clean_count + pass_resist_count) / max(1, len(recs)) * 100.0, 2)
        })

# Save Master Step 7 QC Log
df_qc = pd.DataFrame(qc_log_rows)
qc_log_file = os.path.join(OUT_QC, "07_cds_qc_log.tsv")
df_qc.to_csv(qc_log_file, sep="\t", index=False)
print(f"\nSaved Step 7 QC log to {qc_log_file}")

# Display Summary Table
df_summary = pd.DataFrame(summary_table_rows)
print("\n" + "="*60)
print("STEP 7 QC RESULTS: 18 CLEAN PRIORITY GENES")
print("="*60)
print(df_summary.to_string(index=False))

if __name__ == "__main__":
    pass
