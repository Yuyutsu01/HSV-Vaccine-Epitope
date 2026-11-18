"""
Step 6d: Nucleotide Diagnostic for RL2 and RS1 Failures (Part 2 of Path D).
Evaluates whether translation failures in RL2 and RS1 are true clinical/assembly
indels/frameshifts (REAL_INDEL) or coordinate extraction artifacts (EXTRACTION_ERROR).
"""

import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

# File paths
REF1_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001806.2.gbk"
REF2_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001798.2.gbk"
VAL_TSV = "hsv_curated_dataset/07_qc_reports/06b_translation_validation_per_gene.tsv"
OUT_TSV = "hsv_curated_dataset/07_qc_reports/06d_rl2_rs1_nucleotide_diagnostic.tsv"

# 1. Load Reference Nucleotide CDS and Proteins
def get_reference_models(gbk_path):
    record = SeqIO.read(gbk_path, "genbank")
    full_seq = str(record.seq).upper()
    cds_dict = {}
    prot_dict = {}
    seen = set()
    for feat in record.features:
        if feat.type == "CDS":
            g = feat.qualifiers.get("gene", [""])[0] or feat.qualifiers.get("locus_tag", [""])[0]
            if not g or g in seen:
                continue
            seen.add(g)
            t = feat.qualifiers.get("translation", [""])[0]
            prot_dict[g] = t
            parts = feat.location.parts
            strand = "+" if feat.location.strand == 1 else "-"
            seqs = []
            for p in parts:
                p_s = full_seq[p.start:p.end]
                if strand == "-":
                    p_s = str(Seq(p_s).reverse_complement())
                seqs.append(p_s)
            cds_dict[g] = "".join(seqs)
    return cds_dict, prot_dict

ref_cds_h1, ref_prots_h1 = get_reference_models(REF1_GBK)
ref_cds_h2, ref_prots_h2 = get_reference_models(REF2_GBK)

# Protein Aligner for calculating aa_identity
prot_aligner = PairwiseAligner()
prot_aligner.mode = "global"
prot_aligner.match_score = 2
prot_aligner.mismatch_score = -1
prot_aligner.open_gap_score = -4
prot_aligner.extend_gap_score = -1

# Local Nucleotide Aligner
nt_aligner = PairwiseAligner()
nt_aligner.mode = "local"
nt_aligner.match_score = 2
nt_aligner.mismatch_score = -1
nt_aligner.open_gap_score = -3
nt_aligner.extend_gap_score = -1

# Target sets: (Virus, Gene, FASTA file, Ref CDS, Ref Prot)
targets = [
    ("HSV-1", "RS1", "hsv_curated_dataset/04_cds_unaligned/HSV1/RS1.fasta", ref_cds_h1["RS1"], ref_prots_h1["RS1"]),
    ("HSV-2", "RL2", "hsv_curated_dataset/04_cds_unaligned/HSV2/RL2.fasta", ref_cds_h2["RL2"], ref_prots_h2["RL2"]),
    ("HSV-2", "RS1", "hsv_curated_dataset/04_cds_unaligned/HSV2/RS1.fasta", ref_cds_h2["RS1"], ref_prots_h2["RS1"])
]

diagnostic_rows = []

print("="*60)
print("EXECUTING PART 2: RL2 / RS1 NUCLEOTIDE DIAGNOSTIC")
print("="*60)

for vtype, gname, fasta_path, ref_nt, ref_aa in targets:
    print(f"\nAnalyzing {vtype} {gname} FAIL_ANCHOR_ERROR sequences from {fasta_path}...")
    ref_nt_len = len(ref_nt)
    ref_aa_len = len(ref_aa)
    
    recs = list(SeqIO.parse(fasta_path, "fasta"))
    
    for rec in recs:
        acc = rec.id.split("|")[0]
        s_cds = str(rec.seq).upper()
        
        # 1. Recompute Step 6.5 classification to isolate FAIL_ANCHOR_ERROR
        valid_len = len(s_cds) - (len(s_cds) % 3)
        s_clean = s_cds[:valid_len]
        p_to_stop = str(Seq(s_clean).translate(to_stop=True))
        p_full = str(Seq(s_clean).translate())
        internal_stops = p_full[:-1].count("*") if len(p_full) > 1 else 0
        trans_len = len(p_to_stop)
        length_ratio = trans_len / max(1, ref_aa_len)
        
        # AA identity
        if p_to_stop == ref_aa:
            aa_id = 100.0
        else:
            aln_p = next(iter(prot_aligner.align(p_to_stop, ref_aa)))
            matches_p = sum(c1 == c2 for c1, c2 in zip(aln_p[0], aln_p[1]))
            aa_id = round(matches_p / max(1, len(aln_p[0])) * 100.0, 2)
            
        is_fail_anchor = (aa_id < 80.0) or (length_ratio < 0.5) or (length_ratio > 1.2) or (internal_stops > 1)
        
        if not is_fail_anchor:
            continue
            
        # 2. Local Nucleotide Alignment to Reference CDS
        aln_nt = next(iter(nt_aligner.align(s_cds, ref_nt)))
        matches_nt = sum(c1 == c2 for c1, c2 in zip(aln_nt[0], aln_nt[1]))
        aln_nt_len = len(aln_nt[0])
        nt_identity = round(matches_nt / max(1, aln_nt_len) * 100.0, 2)
        
        # 3. Classify:
        # REAL_INDEL: nt_identity >= 95%
        # EXTRACTION_ERROR: nt_identity < 90%
        # AMBIGUOUS: 90% <= nt_identity < 95%
        if nt_identity >= 95.0:
            classification = "REAL_INDEL"
        elif nt_identity < 90.0:
            classification = "EXTRACTION_ERROR"
        else:
            classification = "AMBIGUOUS"
            
        diagnostic_rows.append({
            "virus": vtype,
            "gene": gname,
            "accession": acc,
            "cds_length": len(s_cds),
            "ref_cds_length": ref_nt_len,
            "nt_identity": nt_identity,
            "aa_identity_from_6.5": aa_id,
            "classification": classification
        })

df_diag = pd.DataFrame(diagnostic_rows)
df_diag.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\nSaved nucleotide diagnostic to {OUT_TSV}")
print(f"Total FAIL_ANCHOR_ERROR sequences evaluated: {len(df_diag)}")

# Print Summary Table
print("\n" + "="*60)
print("PART 2 SUMMARY: NUCLEOTIDE DIAGNOSTIC PER GENE")
print("="*60)

summary_rows = []
for (vtype, gname), grp in df_diag.groupby(["virus", "gene"]):
    tot = len(grp)
    c_real = sum(grp["classification"] == "REAL_INDEL")
    c_ext = sum(grp["classification"] == "EXTRACTION_ERROR")
    c_amb = sum(grp["classification"] == "AMBIGUOUS")
    pct_real = round(c_real / tot * 100.0, 2)
    pct_ext = round(c_ext / tot * 100.0, 2)
    pct_amb = round(c_amb / tot * 100.0, 2)
    summary_rows.append({
        "virus": vtype,
        "gene": gname,
        "total_failures": tot,
        "REAL_INDEL": c_real,
        "REAL_INDEL_pct": pct_real,
        "EXTRACTION_ERROR": c_ext,
        "EXTRACTION_ERROR_pct": pct_ext,
        "AMBIGUOUS": c_amb,
        "AMBIGUOUS_pct": pct_amb
    })

df_sum = pd.DataFrame(summary_rows)
print(df_sum.to_string(index=False))

if __name__ == "__main__":
    pass
