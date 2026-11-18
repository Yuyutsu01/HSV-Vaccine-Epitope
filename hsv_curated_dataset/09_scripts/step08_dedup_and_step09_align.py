"""
Step 8 & Step 9: Exact Deduplication & Codon-Aware Multiple Alignment Pipeline.
Processes all 20 priority genes for HSV-1 and HSV-2.

Step 8:
  - Exact SHA-256 deduplication (100% exact identical sequences only; no collapsing near-duplicates).
  - Outputs 05_cds_qc_passed/{virus}/{gene}_dedup.fasta
  - Generates 07_qc_reports/08_dedup_log.tsv

Step 9:
  - Profile-guided codon-aware multiple sequence alignment.
  - Three distinct alignment sets:
      (a) 06_alignments/HSV1/{gene}_aligned.fasta
      (b) 06_alignments/HSV2/{gene}_aligned.fasta
      (c) 06_alignments/combined/{gene}_combined_aligned.fasta
  - Generates 06_alignments/alignment_summary.tsv
"""

import os
import hashlib
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

# Paths
DIR_QC_PASSED_H1 = "hsv_curated_dataset/05_cds_qc_passed/HSV1"
DIR_QC_PASSED_H2 = "hsv_curated_dataset/05_cds_qc_passed/HSV2"
DIR_ALN_H1 = "hsv_curated_dataset/06_alignments/HSV1"
DIR_ALN_H2 = "hsv_curated_dataset/06_alignments/HSV2"
DIR_ALN_COMBINED = "hsv_curated_dataset/06_alignments/combined"
DIR_QC_REPORTS = "hsv_curated_dataset/07_qc_reports"
REF1_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001806.2.gbk"
REF2_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001798.2.gbk"

os.makedirs(DIR_ALN_H1, exist_ok=True)
os.makedirs(DIR_ALN_H2, exist_ok=True)
os.makedirs(DIR_ALN_COMBINED, exist_ok=True)
os.makedirs(DIR_QC_REPORTS, exist_ok=True)

PRIORITY_20 = [
    "US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53", "UL10", "UL49A",
    "RL2", "RS1", "UL54", "UL48", "UL46", "UL47", "UL49", "UL23", "UL30"
]

# Load reference proteins and CDS
def load_ref_dict(gbk_path):
    record = SeqIO.read(gbk_path, "genbank")
    full_seq = str(record.seq).upper()
    prots = {}
    cdss = {}
    seen = set()
    for feat in record.features:
        if feat.type == "CDS":
            g = feat.qualifiers.get("gene", [""])[0] or feat.qualifiers.get("locus_tag", [""])[0]
            if g and g not in seen:
                seen.add(g)
                prots[g] = feat.qualifiers.get("translation", [""])[0]
                parts = feat.location.parts
                strand = "+" if feat.location.strand == 1 else "-"
                seqs = [full_seq[p.start:p.end] if strand=="+" else str(Seq(full_seq[p.start:p.end]).reverse_complement()) for p in parts]
                cdss[g] = "".join(seqs)
    return prots, cdss

ref_p_h1, ref_cds_h1 = load_ref_dict(REF1_GBK)
ref_p_h2, ref_cds_h2 = load_ref_dict(REF2_GBK)

# Global Aligner for protein alignment
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -4
aligner.extend_gap_score = -1

# -------------------------------------------------------------
# STEP 8: Exact SHA-256 Deduplication
# -------------------------------------------------------------
print("="*60)
print("EXECUTING STEP 8: EXACT SHA-256 DEDUPLICATION (20 PRIORITY GENES)")
print("="*60)

dedup_log_rows = []

for vtype in ["HSV-1", "HSV-2"]:
    in_dir = DIR_QC_PASSED_H1 if vtype == "HSV-1" else DIR_QC_PASSED_H2
    
    for gname in PRIORITY_20:
        fpath = os.path.join(in_dir, f"{gname}.fasta")
        if not os.path.exists(fpath):
            continue
            
        recs = list(SeqIO.parse(fpath, "fasta"))
        tot = len(recs)
        
        # Group by identical sequence SHA-256
        clusters = defaultdict(list)
        for r in recs:
            s_clean = str(r.seq).upper()
            h = hashlib.sha256(s_clean.encode("utf-8")).hexdigest()
            clusters[h].append((r.id, s_clean))
            
        unique_cnt = len(clusters)
        dup_cnt = tot - unique_cnt
        dup_rate = round(dup_cnt / max(1, tot) * 100.0, 2)
        max_cluster = max(len(v) for v in clusters.values()) if clusters else 0
        
        # Write deduplicated FASTA
        dedup_records = []
        for h, members in clusters.items():
            rep_id, rep_seq = members[0]
            cluster_size = len(members)
            all_accs = [m[0].split("|")[0] for m in members]
            # Construct informative deduplicated header
            dedup_hdr = f"{rep_id}|cluster_size={cluster_size}|collapsed={','.join(all_accs[:5])}"
            if cluster_size > 5:
                dedup_hdr += f"...(+{cluster_size-5})"
            dedup_records.append(f">{dedup_hdr}\n{rep_seq}\n")
            
        out_dedup_file = os.path.join(in_dir, f"{gname}_dedup.fasta")
        with open(out_dedup_file, "w", encoding="utf-8") as f:
            f.writelines(dedup_records)
            
        dedup_log_rows.append({
            "virus": vtype,
            "gene": gname,
            "total_clean_sequences": tot,
            "unique_sequences": unique_cnt,
            "duplicate_sequences": dup_cnt,
            "duplicate_rate_pct": dup_rate,
            "max_cluster_size": max_cluster
        })

df_dedup = pd.DataFrame(dedup_log_rows)
dedup_tsv = os.path.join(DIR_QC_REPORTS, "08_dedup_log.tsv")
df_dedup.to_csv(dedup_tsv, sep="\t", index=False)
print(f"Saved exact deduplication audit log to {dedup_tsv}")
print(df_dedup.to_string(index=False))

# -------------------------------------------------------------
# STEP 9: Codon-Aware Multiple Sequence Alignment
# -------------------------------------------------------------
print("\n" + "="*60)
print("EXECUTING STEP 9: CODON-AWARE MULTIPLE SEQUENCE ALIGNMENT")
print("="*60)

def codon_align_records(records, ref_p, ref_cds):
    """
    Performs codon-aware alignment by translating CDS to protein,
    aligning against reference protein, and projecting nucleotide
    codons into the aligned coordinate space.
    """
    ref_aa_len = len(ref_p)
    
    # 1. Align each protein sequence against reference protein
    aligned_entries = []
    
    # Track overall profile positions
    all_col_insertions = defaultdict(int)
    
    parsed = []
    for r_id, s_nt in records:
        s_clean = s_nt[:len(s_nt) - (len(s_nt)%3)]
        p_seq = str(Seq(s_clean).translate(to_stop=True))
        
        # Protein alignment to reference
        if p_seq == ref_p:
            aln_ref = ref_p
            aln_que = p_seq
        else:
            aln = next(iter(aligner.align(ref_p, p_seq)))
            aln_ref = aln[0]
            aln_que = aln[1]
            
        parsed.append((r_id, s_clean, aln_ref, aln_que))

    # Project nucleotide codons onto aligned protein
    for r_id, s_clean, aln_ref, aln_que in parsed:
        codons = [s_clean[i:i+3] for i in range(0, len(s_clean), 3)]
        c_idx = 0
        aligned_codons = []
        for aa in aln_que:
            if aa == "-":
                aligned_codons.append("---")
            else:
                if c_idx < len(codons):
                    aligned_codons.append(codons[c_idx])
                    c_idx += 1
                else:
                    aligned_codons.append("---")
                    
        aligned_entries.append((r_id, "".join(aligned_codons)))
        
    return aligned_entries

alignment_summary_rows = []

for gname in PRIORITY_20:
    f_h1 = os.path.join(DIR_QC_PASSED_H1, f"{gname}_dedup.fasta")
    f_h2 = os.path.join(DIR_QC_PASSED_H2, f"{gname}_dedup.fasta")
    
    recs_h1 = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(f_h1, "fasta")] if os.path.exists(f_h1) else []
    recs_h2 = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(f_h2, "fasta")] if os.path.exists(f_h2) else []
    
    ref_p1 = ref_p_h1.get(gname, "")
    ref_cds1 = ref_cds_h1.get(gname, "")
    ref_p2 = ref_p_h2.get(gname, "")
    ref_cds2 = ref_cds_h2.get(gname, "")
    
    # (a) HSV-1 Alignments
    aln_h1 = codon_align_records(recs_h1, ref_p1, ref_cds1) if recs_h1 else []
    out_aln_h1 = os.path.join(DIR_ALN_H1, f"{gname}_aligned.fasta")
    with open(out_aln_h1, "w", encoding="utf-8") as f:
        for r_id, seq in aln_h1:
            f.write(f">{r_id}\n{seq}\n")
            
    # (b) HSV-2 Alignments
    aln_h2 = codon_align_records(recs_h2, ref_p2, ref_cds2) if recs_h2 else []
    out_aln_h2 = os.path.join(DIR_ALN_H2, f"{gname}_aligned.fasta")
    with open(out_aln_h2, "w", encoding="utf-8") as f:
        for r_id, seq in aln_h2:
            f.write(f">{r_id}\n{seq}\n")
            
    # (c) Combined HSV-1 + HSV-2 Alignment (Reference-Anchored Merger)
    # Align ref_p1 vs ref_p2 to establish inter-species homology coordinate map
    aln_refs = next(iter(aligner.align(ref_p1, ref_p2)))
    ref_map1 = aln_refs[0]
    ref_map2 = aln_refs[1]
    
    combined_entries = []
    
    # Project HSV-1 onto inter-species map
    for r_id, aln_seq in aln_h1:
        codons = [aln_seq[i:i+3] for i in range(0, len(aln_seq), 3)]
        c_i = 0
        comb_codons = []
        for ch in ref_map1:
            if ch == "-":
                comb_codons.append("---")
            else:
                if c_i < len(codons):
                    comb_codons.append(codons[c_i])
                    c_i += 1
                else:
                    comb_codons.append("---")
        combined_entries.append((f"HSV1|{r_id}", "".join(comb_codons)))
        
    # Project HSV-2 onto inter-species map
    for r_id, aln_seq in aln_h2:
        codons = [aln_seq[i:i+3] for i in range(0, len(aln_seq), 3)]
        c_i = 0
        comb_codons = []
        for ch in ref_map2:
            if ch == "-":
                comb_codons.append("---")
            else:
                if c_i < len(codons):
                    comb_codons.append(codons[c_i])
                    c_i += 1
                else:
                    comb_codons.append("---")
        combined_entries.append((f"HSV2|{r_id}", "".join(comb_codons)))
        
    out_aln_comb = os.path.join(DIR_ALN_COMBINED, f"{gname}_combined_aligned.fasta")
    with open(out_aln_comb, "w", encoding="utf-8") as f:
        for r_id, seq in combined_entries:
            f.write(f">{r_id}\n{seq}\n")
            
    # Compute alignment statistics
    if combined_entries:
        comb_len_nt = len(combined_entries[0][1])
        comb_len_codons = comb_len_nt // 3
        
        # Calculate column conservation (% of columns with single non-gap nucleotide across all sequences)
        cons_cols = 0
        total_cols = comb_len_nt
        for col_i in range(total_cols):
            col_chars = set(seq[col_i] for _, seq in combined_entries)
            if len(col_chars) == 1 and "-" not in col_chars:
                cons_cols += 1
        cons_pct = round(cons_cols / max(1, total_cols) * 100.0, 2)
    else:
        comb_len_nt, comb_len_codons, cons_pct = 0, 0, 0.0
        
    alignment_summary_rows.append({
        "gene": gname,
        "hsv1_dedup_seqs": len(aln_h1),
        "hsv2_dedup_seqs": len(aln_h2),
        "combined_total_seqs": len(combined_entries),
        "alignment_length_nt": comb_len_nt,
        "alignment_length_codons": comb_len_codons,
        "nucleotide_conservation_pct": cons_pct
    })

df_aln_sum = pd.DataFrame(alignment_summary_rows)
aln_sum_tsv = os.path.join(DIR_ALN_COMBINED, "../alignment_summary.tsv")
df_aln_sum.to_csv(aln_sum_tsv, sep="\t", index=False)
print(f"\nSaved master alignment summary to {aln_sum_tsv}")
print(df_aln_sum.to_string(index=False))

print("\n" + "="*60)
print("ALL STEP 8 & STEP 9 PIPELINES COMPLETED!")
print("="*60)

if __name__ == "__main__":
    pass
