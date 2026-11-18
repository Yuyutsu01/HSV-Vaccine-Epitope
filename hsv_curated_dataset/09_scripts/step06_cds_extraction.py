"""
Step 6: Exon-Aware, Annotation-Verified CDS Extraction Pipeline.
Extracts all canonical CDS from Tier 1 (HSV-1, HSV-2) and Tier 2 (HSV-2 non-repeat) cohorts.
Uses high-performance sliding terminal k-mer anchor matching with codon boundary snapping,
with local dynamic programming fallback for divergent isolates.
Handles spliced genes (RL2, UL15) via exon-aware concatenation.
Enforces strictly single-copy extraction (Copy 1) for diploid repeat genes (RL1, RL2, RS1).
"""

import os
import re
import time
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

# -------------------------------------------------------------
# Input and Output File Paths
# -------------------------------------------------------------
REF1_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001806.2.gbk"
REF2_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001798.2.gbk"
TIER1_HSV1 = "hsv_curated_dataset/02_filtered_genomes/HSV1_filtered.fasta"
TIER1_HSV2 = "hsv_curated_dataset/02_filtered_genomes/HSV2_filtered.fasta"
TIER2_HSV2 = "hsv_curated_dataset/02_filtered_genomes/HSV2_filtered_tier2.fasta"
CURATED_META = "hsv_curated_dataset/03_metadata/curated_metadata.tsv"

OUT_CDS_HSV1 = "hsv_curated_dataset/04_cds_unaligned/HSV1"
OUT_CDS_HSV2 = "hsv_curated_dataset/04_cds_unaligned/HSV2"
OUT_QC = "hsv_curated_dataset/07_qc_reports"
os.makedirs(OUT_CDS_HSV1, exist_ok=True)
os.makedirs(OUT_CDS_HSV2, exist_ok=True)
os.makedirs(OUT_QC, exist_ok=True)

# -------------------------------------------------------------
# Priority Target Genes (20 genes)
# -------------------------------------------------------------
PRIORITY_GENES = [
    "US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53", "UL10", "UL49A",
    "RL2", "RS1", "UL54", "UL48", "UL46", "UL47", "UL49", "UL23", "UL30"
]

# -------------------------------------------------------------
# 1. Parse Reference Gene Models
# -------------------------------------------------------------
def load_ref_gene_models(gbk_path):
    """
    Parses GenBank reference record to extract exact exon coordinates,
    splicing structure, strand, and reference nucleotide sequences.
    For diploid inverted repeat genes (RL1, RL2, RS1), records Copy 1 only.
    """
    record = SeqIO.read(gbk_path, "genbank")
    full_seq = str(record.seq).upper()
    
    gene_models = {}
    seen_genes = set()
    
    for feat in record.features:
        if feat.type == "CDS":
            g = feat.qualifiers.get("gene", [""])[0]
            locus = feat.qualifiers.get("locus_tag", [""])[0]
            gname = g or locus
            
            # Enforce single copy extraction for repeat genes (Copy 1)
            if gname in seen_genes:
                continue
            seen_genes.add(gname)
            
            strand = "+" if feat.location.strand == 1 else "-"
            parts = []
            for part in feat.location.parts:
                p_start = int(part.start)
                p_end = int(part.end)
                p_seq = full_seq[p_start:p_end]
                if strand == "-":
                    p_seq = str(Seq(p_seq).reverse_complement())
                parts.append({
                    "start": p_start,
                    "end": p_end,
                    "seq": p_seq
                })
                
            full_cds = "".join(p["seq"] for p in parts)
            
            gene_models[gname] = {
                "gene": gname,
                "strand": strand,
                "is_spliced": len(parts) > 1,
                "parts": parts,
                "full_cds": full_cds,
                "ref_length": len(full_cds),
                "is_repeat": gname in ["RL1", "RL2", "RS1"]
            }
            
    return gene_models

# Fallback PairwiseAligner for narrow window alignment
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -3
aligner.extend_gap_score = -1

# -------------------------------------------------------------
# 2. Exon Extraction Function
# -------------------------------------------------------------
def extract_single_exon(genome_seq, ref_exon_seq, p_start, p_end, strand, is_start_exon, is_end_exon):
    """
    Extracts an individual exon from target genome using sliding terminal anchors,
    falling back to local dynamic programming on narrow boundary windows.
    """
    seq_l = len(genome_seq)
    pad = 2500  # Wide enough to accommodate any cumulative genome coordinate drift
    w_start = max(0, p_start - pad)
    w_end = min(seq_l, p_end + pad)
    win = genome_seq[w_start:w_end]
    
    # Orient window to 5' -> 3' coding direction
    if strand == "-":
        win = str(Seq(win).reverse_complement())
        
    exon_l = len(ref_exon_seq)
    
    # 1. High-speed sliding terminal anchor matching (18-mers)
    pos_s = -1
    pos_e = -1
    
    # Start anchor
    if is_start_exon:
        # Must snap to ATG
        for off in range(0, min(90, exon_l - 18), 3):
            k = ref_exon_seq[off:off+18]
            idx = win.find(k)
            if idx != -1:
                cand = idx - off
                if 0 <= cand < len(win) and win[cand:cand+3] == "ATG":
                    pos_s = cand
                    break
    else:
        for off in range(0, min(90, exon_l - 18), 3):
            k = ref_exon_seq[off:off+18]
            idx = win.find(k)
            if idx != -1:
                cand = idx - off
                if 0 <= cand < len(win):
                    pos_s = cand
                    break
                    
    # End anchor
    if is_end_exon:
        # Must snap to in-frame stop codon
        for off in range(0, min(90, exon_l - 18), 3):
            k = ref_exon_seq[exon_l - 18 - off : exon_l - off]
            idx = win.find(k)
            if idx != -1:
                cand = idx + 18 + off
                if cand <= len(win) and win[cand-3:cand] in ["TAA", "TAG", "TGA"]:
                    pos_e = cand
                    break
    else:
        for off in range(0, min(90, exon_l - 18), 3):
            k = ref_exon_seq[exon_l - 18 - off : exon_l - off]
            idx = win.find(k)
            if idx != -1:
                cand = idx + 18 + off
                if cand <= len(win):
                    pos_e = cand
                    break
                    
    # Validate anchor-derived exon
    if pos_s != -1 and pos_e != -1 and pos_e > pos_s:
        cand_exon = win[pos_s:pos_e]
        if abs(len(cand_exon) - exon_l) <= 120:
            return cand_exon
            
    # 2. Fallback: Local alignment on candidate sub-window if terminal anchor missed
    # Find rough location via internal conserved 24-mer
    mid_k = ref_exon_seq[max(0, exon_l//2 - 12) : max(0, exon_l//2 + 12)]
    idx_mid = win.find(mid_k)
    if idx_mid != -1:
        # Target sub-region around predicted bounds
        approx_s = max(0, idx_mid - (exon_l // 2) - 150)
        approx_e = min(len(win), idx_mid + (exon_l // 2) + 150)
        sub_win = win[approx_s:approx_e]
        
        alns = aligner.align(sub_win, ref_exon_seq)
        top_aln = next(iter(alns), None)
        if top_aln is not None:
            qs, qe = top_aln.aligned[0][0][0], top_aln.aligned[0][-1][1]
            cand = sub_win[qs:qe]
            
            # Snap start to ATG
            if is_start_exon and cand[:3] != "ATG":
                sub_atg = sub_win[max(0, qs-12):min(len(sub_win), qs+12)]
                atg_idx = sub_atg.find("ATG")
                if atg_idx != -1:
                    qs = max(0, qs-12) + atg_idx
                    cand = sub_win[qs:qe]
                    
            # Snap stop
            if is_end_exon and cand[-3:] not in ["TAA", "TAG", "TGA"]:
                for off in [0, 3, -3, 6, -6, 9, -9]:
                    cand_off = sub_win[qs:qe+off]
                    if len(cand_off) % 3 == 0 and cand_off[-3:] in ["TAA", "TAG", "TGA"]:
                        cand = cand_off
                        break
            return cand
            
    return ""

# -------------------------------------------------------------
# 3. Main Extraction Pipeline Function
# -------------------------------------------------------------
def run_extraction_pipeline():
    print("="*60)
    print("STARTING STEP 6: EXON-AWARE CDS EXTRACTION")
    print("="*60)
    
    t_start = time.time()
    
    ref_models_hsv1 = load_ref_gene_models(REF1_GBK)
    ref_models_hsv2 = load_ref_gene_models(REF2_GBK)
    print(f"Loaded {len(ref_models_hsv1)} HSV-1 and {len(ref_models_hsv2)} HSV-2 reference gene models.")
    
    # Curated metadata dictionary for header construction
    df_meta = pd.read_csv(CURATED_META, sep="\t") if os.path.exists(CURATED_META) else pd.DataFrame()
    meta_dict = {}
    if not df_meta.empty:
        for _, row in df_meta.iterrows():
            meta_dict[str(row["accession"]).split(".")[0]] = row.to_dict()
            
    # Load all target genomes
    datasets = [
        ("HSV-1", "Tier 1", TIER1_HSV1, ref_models_hsv1),
        ("HSV-2", "Tier 1", TIER1_HSV2, ref_models_hsv2),
        ("HSV-2", "Tier 2", TIER2_HSV2, ref_models_hsv2)
    ]
    
    gene_fasta_dict = defaultdict(lambda: defaultdict(list))
    log_records = []
    
    total_genomes_processed = 0
    total_extractions = 0
    total_success = 0
    
    for vtype, tier, fasta_path, ref_models in datasets:
        print(f"\nProcessing {vtype} {tier} genomes from {fasta_path}...")
        records = list(SeqIO.parse(fasta_path, "fasta"))
        print(f"Loaded {len(records)} genomes.")
        
        for rec in records:
            total_genomes_processed += 1
            seq_str = str(rec.seq).upper()
            acc = rec.id.split(".")[0]
            meta = meta_dict.get(acc, {})
            
            for gname, model in ref_models.items():
                total_extractions += 1
                
                # Tier 2 HSV-2 genomes: Exclude repeat genes
                if tier == "Tier 2" and model["is_repeat"]:
                    continue
                    
                strand = model["strand"]
                is_spliced = model["is_spliced"]
                ref_len = model["ref_length"]
                parts = model["parts"]
                
                status = "FAILED"
                failure_reason = ""
                extracted_cds = ""
                
                exon_pieces = []
                all_exons_ok = True
                
                for i, part in enumerate(parts):
                    is_start = (i == 0)
                    is_end = (i == len(parts) - 1)
                    ex_seq = extract_single_exon(
                        seq_str, part["seq"], part["start"], part["end"], strand, is_start, is_end
                    )
                    if not ex_seq:
                        all_exons_ok = False
                        failure_reason = f"Exon {i+1} extraction failed"
                        break
                    exon_pieces.append(ex_seq)
                    
                if all_exons_ok:
                    extracted_cds = "".join(exon_pieces)
                    ext_len = len(extracted_cds)
                    
                    # Evaluation
                    if ext_len >= 0.8 * ref_len and ext_len <= 1.2 * ref_len:
                        status = "SUCCESS"
                        total_success += 1
                    else:
                        status = "FAILED"
                        failure_reason = f"Extracted length {ext_len} outside 80-120% of ref length {ref_len}"
                else:
                    status = "FAILED"
                    
                # Build standardized FASTA header
                # Format: >accession|gene|isolate|country|collection_year|BioSample|lab_strain_flag|recombinant_flag
                isolate = str(meta.get("isolate", "")) if pd.notna(meta.get("isolate", "")) else ""
                country = str(meta.get("country", "")) if pd.notna(meta.get("country", "")) else ""
                year = str(meta.get("collection_year", "")) if pd.notna(meta.get("collection_year", "")) else ""
                biosample = str(meta.get("BioSample", "")) if pd.notna(meta.get("BioSample", "")) else ""
                lab_flag = meta.get("lab_strain_flag", False)
                recomb_flag = meta.get("recombinant_flag", False)
                
                fasta_header = f"{acc}|{gname}|{isolate}|{country}|{year}|{biosample}|{lab_flag}|{recomb_flag}"
                
                if status == "SUCCESS" and len(extracted_cds) > 0:
                    gene_fasta_dict[vtype][gname].append(f">{fasta_header}\n{extracted_cds}\n")
                    
                # Append to TSV log record
                log_records.append({
                    "accession": acc,
                    "virus": vtype,
                    "cohort_tier": tier,
                    "gene": gname,
                    "is_priority": gname in PRIORITY_GENES,
                    "is_spliced": is_spliced,
                    "is_repeat": model["is_repeat"],
                    "copy_extracted": 1,
                    "extraction_status": status,
                    "extracted_length": len(extracted_cds),
                    "ref_length": ref_len,
                    "start_codon": extracted_cds[:3] if len(extracted_cds) >= 3 else "",
                    "stop_codon": extracted_cds[-3:] if len(extracted_cds) >= 3 else "",
                    "in_frame": (len(extracted_cds) % 3 == 0) if len(extracted_cds) > 0 else False,
                    "pct_length": round(len(extracted_cds) / max(1, ref_len) * 100.0, 2),
                    "failure_reason": failure_reason,
                    "fasta_header": fasta_header
                })

    # Write per-gene unaligned FASTA files
    print("\nWriting per-gene nucleotide FASTA files...")
    for vtype, genes_dict in gene_fasta_dict.items():
        out_dir = OUT_CDS_HSV1 if vtype == "HSV-1" else OUT_CDS_HSV2
        for gname, lines in genes_dict.items():
            fpath = os.path.join(out_dir, f"{gname}.fasta")
            with open(fpath, "w", encoding="utf-8") as f:
                f.writelines(lines)
                
    # Save complete extraction log
    df_log = pd.DataFrame(log_records)
    log_tsv = os.path.join(OUT_QC, "06_cds_extraction_log.tsv")
    df_log.to_csv(log_tsv, sep="\t", index=False)
    print(f"Saved extraction log to {log_tsv}")
    
    t_end = time.time()
    elapsed = round(t_end - t_start, 2)
    print(f"\nCompleted Step 6 extraction in {elapsed} seconds ({elapsed/60.0:.2f} minutes)!")
    print(f"Total genomes processed: {total_genomes_processed}")
    print(f"Total CDS extraction attempts: {total_extractions}")
    print(f"Successful extractions: {total_success} ({total_success/total_extractions*100:.2f}%)")
    
    # -------------------------------------------------------------
    # Summary of Priority Target Genes
    # -------------------------------------------------------------
    print("\n" + "="*60)
    print("STEP 6 EXTRACTION SUMMARY: PRIORITY GENES")
    print("="*60)
    
    prio_summary = []
    for g in PRIORITY_GENES:
        h1_sub = df_log[(df_log["virus"] == "HSV-1") & (df_log["gene"] == g) & (df_log["extraction_status"] == "SUCCESS")]
        h2_t1 = df_log[(df_log["virus"] == "HSV-2") & (df_log["cohort_tier"] == "Tier 1") & (df_log["gene"] == g) & (df_log["extraction_status"] == "SUCCESS")]
        h2_t2 = df_log[(df_log["virus"] == "HSV-2") & (df_log["cohort_tier"] == "Tier 2") & (df_log["gene"] == g) & (df_log["extraction_status"] == "SUCCESS")]
        
        c_h1 = len(h1_sub)
        c_h2_t1 = len(h2_t1)
        c_h2_t2 = len(h2_t2)
        
        prio_summary.append({
            "Gene": g,
            "HSV1_Tier1_Clean": c_h1,
            "HSV2_Tier1_Clean": c_h2_t1,
            "HSV2_Tier2_Sensitivity": c_h2_t2,
            "Combined_Primary": c_h1 + c_h2_t1,
            "Combined_Total_With_Tier2": c_h1 + c_h2_t1 + c_h2_t2
        })
        
    df_prio = pd.DataFrame(prio_summary)
    print(df_prio.to_string(index=False))

if __name__ == "__main__":
    run_extraction_pipeline()
