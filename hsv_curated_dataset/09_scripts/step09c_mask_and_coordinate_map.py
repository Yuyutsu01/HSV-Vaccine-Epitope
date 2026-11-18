"""
Step 09c (CORRECTION): Alignment QC, principled codon masking, reference coordinate maps.

Supersedes step11_prereq3_clean_masked.py / step11b_masking_investigation.py, which masked
low-complexity tracts in order to hide the out-of-register columns produced by the old star
alignment (CORRECTIONS_AND_PENDING_WORK.md, A5). With a genuine MSA in place, masking is now
needed only where a repeat expansion is genuinely sparse across isolates.

Per alignment:
  1. Codon-column gap fractions, computed over the clinical sequences (reference excluded).
  2. Mask (drop) whole codon columns present in fewer than MIN_OCCUPANCY of clinical sequences.
     Dropping whole codons keeps the reading frame intact for HyPhy.
  3. Flag outlier sequences by distance to the column-wise consensus (median + k*MAD),
     for investigation -- NOT auto-removal, since a true recombinant is biology.
  4. Emit a column -> reference-residue map, so every downstream selection site and
     predicted epitope can be reported in reference numbering.

Outputs:
  - 06_alignments_v2/masked/{HSV1,HSV2,combined}/<gene>[_combined]_masked.fasta
  - 06_alignments_v2/coordinate_maps/<cohort>/<gene>_colmap.tsv
  - 06_alignments_v2/masking_report_v2.tsv
  - 06_alignments_v2/outlier_sequences_v2.tsv
  - 06_alignments_v2/masking.log
"""

import os
import sys
import glob
import random
import itertools
from datetime import datetime

import pandas as pd

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_V2 = os.path.join(DIR_BASE, "06_alignments_v2")
DIR_MASKED = os.path.join(DIR_V2, "masked")
DIR_MAPS = os.path.join(DIR_V2, "coordinate_maps")

# A codon column is kept when at least this fraction of clinical sequences occupy it.
MIN_OCCUPANCY = 0.50
# Outlier flag threshold, in MADs above the median consensus distance.
OUTLIER_K = 6.0

LOG_PATH = os.path.join(DIR_V2, "masking.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def read_fasta(path):
    ids, seqs, cur = [], [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    seqs.append(cur)
                ids.append(line[1:])
                cur = ""
            elif cur is not None:
                cur += line
    if cur is not None:
        seqs.append(cur)
    return ids, seqs


def write_fasta(ids, seqs, path):
    with open(path, "w", encoding="utf-8") as fh:
        for i, s in zip(ids, seqs):
            fh.write(f">{i}\n{s}\n")


def codonise(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


def ranges_from_indices(idx):
    """Compact a sorted index list into 1-based inclusive ranges."""
    out = []
    for i in idx:
        if out and i == out[-1][1]:
            out[-1][1] = i + 1
        else:
            out.append([i, i + 1])
    return ";".join(f"{a + 1}-{b}" for a, b in out)


def mean_pairwise_identity(seqs, sample=60, seed=0):
    if len(seqs) < 2:
        return float("nan")
    random.seed(seed)
    idx = random.sample(range(len(seqs)), min(sample, len(seqs)))
    tot, cnt = 0.0, 0
    for a, b in itertools.combinations(idx, 2):
        A, B = seqs[a], seqs[b]
        m = sum(1 for x, y in zip(A, B) if x == y and x != "-")
        d = sum(1 for x, y in zip(A, B) if x != "-" and y != "-")
        if d:
            tot += m / d
            cnt += 1
    return round(100 * tot / cnt, 2) if cnt else float("nan")


def consensus_outliers(ids, seqs):
    """Distance of each sequence to the column-wise consensus; flag by MAD."""
    L = len(seqs[0])
    cons = []
    for i in range(L):
        col = [s[i] for s in seqs if s[i] != "-"]
        cons.append(max(set(col), key=col.count) if col else "-")
    cons = "".join(cons)

    dist = []
    for s in seqs:
        d = sum(1 for x, y in zip(s, cons) if x != y and x != "-" and y != "-")
        n = sum(1 for x, y in zip(s, cons) if x != "-" and y != "-")
        dist.append(d / n if n else 0.0)

    srt = sorted(dist)
    med = srt[len(srt) // 2]
    absdev = sorted(abs(d - med) for d in dist)
    mad = absdev[len(absdev) // 2] or 1e-9
    cutoff = med + OUTLIER_K * mad

    return [(i, round(d, 5), round(cutoff, 5))
            for i, d in zip(ids, dist) if d > cutoff], round(med, 5), round(cutoff, 5)


def process(path, cohort, gene):
    ids, seqs = read_fasta(path)
    if len(seqs) < 2:
        return None, None, None

    cod = [codonise(s) for s in seqs]
    ncol = len(cod[0])

    ref_rows = [k for k, i in enumerate(ids) if i.startswith("REF_")]
    clinical = [k for k in range(len(ids)) if k not in ref_rows]

    # Occupancy over clinical sequences only.
    keep, drop = [], []
    for c in range(ncol):
        occupied = sum(1 for k in clinical if cod[k][c] != "---")
        (keep if occupied / len(clinical) >= MIN_OCCUPANCY else drop).append(c)

    masked = ["".join(row[c] for c in keep) for row in cod]

    # Column -> reference residue map (on the unmasked alignment).
    map_rows = []
    if ref_rows:
        rcod = cod[ref_rows[0]]
        ref_res = 0
        for c in range(ncol):
            present = rcod[c] != "---"
            if present:
                ref_res += 1
            map_rows.append(dict(
                alignment_codon_col=c + 1,
                ref_residue=ref_res if present else pd.NA,
                ref_present=present,
                clinical_occupancy=round(
                    sum(1 for k in clinical if cod[k][c] != "---") / len(clinical), 4),
                retained_after_masking=c in keep,
                masked_codon_col=(keep.index(c) + 1) if c in keep else pd.NA,
            ))

    outliers, med, cutoff = consensus_outliers(ids, masked)

    report = dict(
        gene=gene, cohort=cohort, n_sequences=len(ids), n_clinical=len(clinical),
        unmasked_codons=ncol, masked_codons=len(keep), codons_dropped=len(drop),
        pct_codons_dropped=round(100 * len(drop) / ncol, 2),
        dropped_codon_ranges=ranges_from_indices(drop) if drop else "None",
        identity_before=mean_pairwise_identity(seqs),
        identity_after=mean_pairwise_identity(masked),
        consensus_dist_median=med, outlier_cutoff=cutoff, n_outliers=len(outliers),
    )
    out_rows = [dict(gene=gene, cohort=cohort, sequence_id=i,
                     consensus_distance=d, cutoff=c) for i, d, c in outliers]
    return (ids, masked), report, (report, map_rows, out_rows)


def main():
    reports, outlier_rows = [], []

    for cohort in ("HSV1", "HSV2", "combined"):
        os.makedirs(os.path.join(DIR_MASKED, cohort), exist_ok=True)
        os.makedirs(os.path.join(DIR_MAPS, cohort), exist_ok=True)

    log("=" * 78)
    log(f"STEP 09c: masking (min clinical occupancy {MIN_OCCUPANCY:.0%}) "
        f"and reference coordinate maps")

    for cohort in ("HSV1", "HSV2", "combined"):
        pattern = os.path.join(DIR_V2, cohort, "*_aligned.fasta")
        for path in sorted(glob.glob(pattern)):
            base = os.path.basename(path)
            gene = base.replace("_combined_aligned.fasta", "").replace("_aligned.fasta", "")

            res, report, extras = process(path, cohort, gene)
            if res is None:
                log(f"{gene} {cohort}: too few sequences, skipped")
                continue
            ids, masked = res
            _r, map_rows, out_rows = extras

            suffix = "_combined_masked.fasta" if cohort == "combined" else "_masked.fasta"
            write_fasta(ids, masked, os.path.join(DIR_MASKED, cohort, gene + suffix))

            if map_rows:
                pd.DataFrame(map_rows).to_csv(
                    os.path.join(DIR_MAPS, cohort, f"{gene}_colmap.tsv"),
                    sep="\t", index=False)

            reports.append(report)
            outlier_rows.extend(out_rows)

            log(f"{gene} {cohort}: {report['unmasked_codons']} -> {report['masked_codons']} codons "
                f"({report['pct_codons_dropped']}% dropped), "
                f"ident {report['identity_before']}% -> {report['identity_after']}%, "
                f"{report['n_outliers']} outlier(s)")

    pd.DataFrame(reports).to_csv(
        os.path.join(DIR_V2, "masking_report_v2.tsv"), sep="\t", index=False)
    pd.DataFrame(outlier_rows).to_csv(
        os.path.join(DIR_V2, "outlier_sequences_v2.tsv"), sep="\t", index=False)

    log(f"Wrote masking_report_v2.tsv ({len(reports)} alignments) and "
        f"outlier_sequences_v2.tsv ({len(outlier_rows)} flagged sequences)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
