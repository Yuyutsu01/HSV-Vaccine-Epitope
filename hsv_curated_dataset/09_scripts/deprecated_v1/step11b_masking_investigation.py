"""
Step 11b: Tier B Masking Investigation (RL2 and RS1).
Investigates low-complexity regions, clinical repeat-length polymorphism,
generates codon-masked alignments, and validates within-type AA identity.

Outputs:
  - 07_qc_reports/11b_tier_b_masking_investigation.tsv
  - 06_alignments/masked/HSV1/{gene}_masked.fasta
  - 06_alignments/masked/HSV2/{gene}_masked.fasta
  - 06_alignments/masked/combined/{gene}_combined_masked.fasta
"""

import os
import math
import random
from collections import Counter
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

# Ensure reproducible sampling
random.seed(42)

# File paths
DIR_QC_PASSED_H1 = "hsv_curated_dataset/05_cds_qc_passed/HSV1"
DIR_QC_PASSED_H2 = "hsv_curated_dataset/05_cds_qc_passed/HSV2"
DIR_ALN_H1 = "hsv_curated_dataset/06_alignments/HSV1"
DIR_ALN_H2 = "hsv_curated_dataset/06_alignments/HSV2"
DIR_ALN_COMB = "hsv_curated_dataset/06_alignments/combined"
DIR_MASKED_H1 = "hsv_curated_dataset/06_alignments/masked/HSV1"
DIR_MASKED_H2 = "hsv_curated_dataset/06_alignments/masked/HSV2"
DIR_MASKED_COMB = "hsv_curated_dataset/06_alignments/masked/combined"
DIR_REPORTS = "hsv_curated_dataset/07_qc_reports"
REF1_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001806.2.gbk"
REF2_GBK = "hsv_curated_dataset/01_reference_genomes/NC_001798.2.gbk"

os.makedirs(DIR_MASKED_H1, exist_ok=True)
os.makedirs(DIR_MASKED_H2, exist_ok=True)
os.makedirs(DIR_MASKED_COMB, exist_ok=True)
os.makedirs(DIR_REPORTS, exist_ok=True)

# Pairwise aligner for global reference-anchored projection
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -4
aligner.extend_gap_score = -1


def shannon_entropy(window: str) -> float:
    """Computes Shannon entropy (in bits) over a window of amino acids."""
    counts = Counter(window)
    length = len(window)
    return -sum((c / length) * math.log2(c / length) for c in counts.values())


def find_low_complexity_regions(seq: str, window_size: int = 20, threshold: float = 2.2):
    """
    Scans a protein sequence with a sliding window to detect low-complexity tracts.
    Positions in any window with Shannon entropy < threshold bits are flagged.
    Returns:
      coords: list of (start_aa, end_aa) tuples (1-based, inclusive)
      fraction: fraction of total protein length classified as low complexity
      is_lc_mask: boolean list of length len(seq)
    """
    is_lc = [False] * len(seq)
    for i in range(len(seq) - window_size + 1):
        win = seq[i : i + window_size]
        if shannon_entropy(win) < threshold:
            for j in range(i, i + window_size):
                is_lc[j] = True

    coords = []
    in_region = False
    start = 0
    for idx, flag in enumerate(is_lc):
        if flag and not in_region:
            in_region = True
            start = idx + 1
        elif not flag and in_region:
            in_region = False
            coords.append((start, idx))
    if in_region:
        coords.append((start, len(seq)))

    fraction = sum(is_lc) / max(1, len(seq))
    return coords, fraction, is_lc


def get_ref_protein(gbk_path: str, gene: str) -> str:
    """Extracts the reference translation for a given gene from a GenBank record."""
    rec = SeqIO.read(gbk_path, "genbank")
    for feat in rec.features:
        if feat.type == "CDS":
            g = feat.qualifiers.get("gene", [""])[0] or feat.qualifiers.get("locus_tag", [""])[0]
            if g == gene:
                return feat.qualifiers.get("translation", [""])[0]
    return ""


def compute_pairwise_identity(seqs: list, n_pairs: int = 500) -> float:
    """
    Computes mean pairwise amino acid identity across aligned sequences.
    Ignores positions where either sequence contains a gap or mask character.
    """
    if len(seqs) < 2:
        return 100.0
    n = len(seqs)
    pairs = []
    if n <= 50:
        for i in range(n):
            for j in range(i + 1, n):
                pairs.append((seqs[i], seqs[j]))
    else:
        for _ in range(n_pairs):
            i, j = random.sample(range(n), 2)
            pairs.append((seqs[i], seqs[j]))

    identities = []
    for s1, s2 in pairs:
        matches = 0
        valid = 0
        for c1, c2 in zip(s1, s2):
            if c1 != "-" and c2 != "-" and c1 != "X" and c2 != "X":
                valid += 1
                if c1 == c2:
                    matches += 1
        if valid > 0:
            identities.append(matches / valid * 100.0)

    return round(float(pd.Series(identities).mean()), 2) if identities else 100.0


print("=" * 70)
print("EXECUTING TIER B MASKING INVESTIGATION (RL2 & RS1)")
print("=" * 70)

# Baseline unmasked values from Step 10 diagnostic (09_conservation_metric_diagnostic.tsv)
BASELINE_UNMASKED = {
    ("HSV-1", "RS1"): 55.72,
    ("HSV-2", "RS1"): 30.76,
    ("HSV-1", "RL2"): 64.18,
    ("HSV-2", "RL2"): 85.15,
}

report_rows = []

for gene in ["RL2", "RS1"]:
    # 1. Reference protein analysis
    ref_h1 = get_ref_protein(REF1_GBK, gene)
    ref_h2 = get_ref_protein(REF2_GBK, gene)

    coords_h1, frac_h1, mask_h1 = find_low_complexity_regions(ref_h1)
    coords_h2, frac_h2, mask_h2 = find_low_complexity_regions(ref_h2)

    # 2. Inter-reference alignment for combined coordinate mapping
    aln_ref = next(iter(aligner.align(ref_h1, ref_h2)))
    ref_map1, ref_map2 = aln_ref[0], aln_ref[1]
    comb_map_len = len(ref_map1)

    # Combined alignment low-complexity column mask
    comb_col_lc = [False] * comb_map_len
    idx1, idx2 = 0, 0
    for col in range(comb_map_len):
        c1, c2 = ref_map1[col], ref_map2[col]
        is_lc = False
        if c1 != "-":
            if idx1 < len(mask_h1) and mask_h1[idx1]:
                is_lc = True
            idx1 += 1
        if c2 != "-":
            if idx2 < len(mask_h2) and mask_h2[idx2]:
                is_lc = True
            idx2 += 1
        if is_lc:
            comb_col_lc[col] = True

    # Process each virus type
    for vtype, ref_p, coords, frac, mask_ref, in_dir, out_dir in [
        ("HSV-1", ref_h1, coords_h1, frac_h1, mask_h1, DIR_QC_PASSED_H1, DIR_MASKED_H1),
        ("HSV-2", ref_h2, coords_h2, frac_h2, mask_h2, DIR_QC_PASSED_H2, DIR_MASKED_H2),
    ]:
        fasta_file = os.path.join(in_dir, f"{gene}.fasta")
        dedup_file = os.path.join(in_dir, f"{gene}_dedup.fasta")
        active_fasta = dedup_file if os.path.exists(dedup_file) else fasta_file

        recs = list(SeqIO.parse(active_fasta, "fasta"))
        n_seqs = len(recs)

        # Amino acid lengths of clinical sequences
        aa_lens = [len(str(r.seq)) // 3 for r in recs]
        q_min = int(min(aa_lens)) if aa_lens else 0
        q25, q_med, q75 = [round(float(x), 1) for x in np.percentile(aa_lens, [25, 50, 75])] if aa_lens else (0, 0, 0)
        q_max = int(max(aa_lens)) if aa_lens else 0

        # Construct reference-anchored per-type codon alignment
        # This properly aligns each sequence to the reference coordinates,
        # preventing repeat expansions from phase-shifting downstream codons.
        aligned_per_type_codons = []
        masked_per_type_codons = []
        translated_clean_unmasked = []
        translated_masked_prots = []

        for r in recs:
            s_clean = str(r.seq).upper()
            s_clean = s_clean[: len(s_clean) - (len(s_clean) % 3)]
            codons = [s_clean[i : i + 3] for i in range(0, len(s_clean), 3)]
            p_seq = str(Seq(s_clean).translate(to_stop=True))

            # Align query protein to reference protein
            if p_seq == ref_p:
                seq_codons = codons[:len(ref_p)]
                seq_masked_codons = ["---" if (idx < len(mask_ref) and mask_ref[idx]) else c for idx, c in enumerate(seq_codons)]
                p_unmasked = ref_p
                p_masked = "".join(["-" if (idx < len(mask_ref) and mask_ref[idx]) else ref_p[idx] for idx in range(len(ref_p))])
            else:
                aln = next(iter(aligner.align(ref_p, p_seq)))
                ref_aln_str, que_aln_str = aln[0], aln[1]

                # Map query codons and amino acids strictly to reference positions
                ref_codons = []
                ref_masked_codons = []
                ref_aa_chars = []
                ref_masked_aa_chars = []
                que_c_idx = 0

                for r_ch, q_ch in zip(ref_aln_str, que_aln_str):
                    if q_ch != "-":
                        curr_c = codons[que_c_idx] if que_c_idx < len(codons) else "---"
                        curr_aa = q_ch
                        que_c_idx += 1
                    else:
                        curr_c = "---"
                        curr_aa = "-"

                    if r_ch != "-":  # Reference position
                        pos_idx = len(ref_codons)
                        is_lc = mask_ref[pos_idx] if pos_idx < len(mask_ref) else False

                        ref_codons.append(curr_c)
                        ref_aa_chars.append(curr_aa)

                        if is_lc:
                            ref_masked_codons.append("---")
                            ref_masked_aa_chars.append("-")
                        else:
                            ref_masked_codons.append(curr_c)
                            ref_masked_aa_chars.append(curr_aa)

                seq_codons = ref_codons
                seq_masked_codons = ref_masked_codons
                p_unmasked = "".join(ref_aa_chars)
                p_masked = "".join(ref_masked_aa_chars)

            aligned_per_type_codons.append((r.id, "".join(seq_codons)))
            masked_per_type_codons.append((r.id, "".join(seq_masked_codons)))
            translated_clean_unmasked.append(p_unmasked)
            translated_masked_prots.append(p_masked)

        # Write per-type masked FASTA
        out_masked_fasta = os.path.join(out_dir, f"{gene}_masked.fasta")
        with open(out_masked_fasta, "w", encoding="utf-8") as f:
            for r_id, seq_str in masked_per_type_codons:
                f.write(f">{r_id}\n{seq_str}\n")

        # Compute within-type AA identities
        masked_id = compute_pairwise_identity(translated_masked_prots)
        clean_unmasked_id = compute_pairwise_identity(translated_clean_unmasked)
        baseline_id = BASELINE_UNMASKED.get((vtype, gene), 0.0)

        # Formatted coordinates string
        coords_str = "; ".join([f"{s}-{e}" for s, e in coords])

        unmasked_aln_len = len(ref_p) * 3
        masked_aln_len = (len(ref_p) - sum(mask_ref)) * 3

        report_rows.append({
            "virus": vtype,
            "gene": gene,
            "ref_length_aa": len(ref_p),
            "low_complexity_coords": coords_str,
            "low_complexity_fraction_pct": round(frac * 100.0, 2),
            "num_clinical_sequences": n_seqs,
            "aa_min": q_min,
            "aa_25th": q25,
            "aa_median": q_med,
            "aa_75th": q75,
            "aa_max": q_max,
            "unmasked_aln_len_codons": len(ref_p),
            "masked_aln_len_codons": len(ref_p) - sum(mask_ref),
            "baseline_step10_unmasked_id_pct": baseline_id,
            "clean_ref_aligned_unmasked_id_pct": clean_unmasked_id,
            "masked_within_type_id_pct": masked_id,
            "artifact_confirmed": (masked_id >= 95.0),
            "verdict": "CONFIRMED_REPEAT_ARTIFACT (>=95% identity)" if masked_id >= 95.0 else "UNRESOLVED_DIVERSITY"
        })

    # 3. Generate the Combined Masked Alignment for downstream HyPhy evaluation
    comb_unmasked_path = os.path.join(DIR_ALN_COMB, f"{gene}_combined_aligned.fasta")
    comb_masked_path = os.path.join(DIR_MASKED_COMB, f"{gene}_combined_masked.fasta")
    if os.path.exists(comb_unmasked_path):
        comb_recs = list(SeqIO.parse(comb_unmasked_path, "fasta"))
        comb_masked_entries = []
        for r in comb_recs:
            s_nt = str(r.seq).upper()
            codons = [s_nt[i : i + 3] for i in range(0, len(s_nt), 3)]
            masked_c = []
            for col_i, c in enumerate(codons):
                if col_i < len(comb_col_lc) and comb_col_lc[col_i]:
                    masked_c.append("---")
                else:
                    masked_c.append(c)
            comb_masked_entries.append((r.id, "".join(masked_c)))

        with open(comb_masked_path, "w", encoding="utf-8") as f:
            for r_id, seq_str in comb_masked_entries:
                f.write(f">{r_id}\n{seq_str}\n")
        print(f"Generated combined masked alignment: {comb_masked_path}")

# Save Tier B investigation table
df_report = pd.DataFrame(report_rows)
out_tsv = os.path.join(DIR_REPORTS, "11b_tier_b_masking_investigation.tsv")
df_report.to_csv(out_tsv, sep="\t", index=False)
print(f"\nSaved Tier B Masking Investigation Report to {out_tsv}")
print(df_report.to_string(index=False))

print("\n" + "=" * 70)
print("TIER B MASKING INVESTIGATION COMPLETED SUCCESSFULLY!")
print("=" * 70)
