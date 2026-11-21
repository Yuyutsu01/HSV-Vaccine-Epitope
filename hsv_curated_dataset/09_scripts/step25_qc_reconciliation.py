"""
Step 25: complete, authoritative CDS QC record (verification pass V2/V3, 2026-11-21).

Why this exists
---------------
`07_qc_reports/07_cds_qc_log.tsv` holds 11,504 assessments, but 12,054 CDS were actually
carried forward into `05_cds_qc_passed/` and from there into deduplication, alignment and
selection. 697 of the carried-forward sequences have no row in that log, because a later
extraction top-up updated the FASTAs and `06_cds_extraction_log.tsv` but not the QC log;
19 more carry a stale `FAIL`. The consequence was that no single file described the QC
state of the data that was actually analysed, and the supplements quoted 11,340/164 as if
it did.

This script re-assesses **every sequence present in the carried-forward FASTAs** against
the same criteria `step07_part1_clean_qc.py` applies, and writes one complete log.

It is deliberately read-only with respect to the analysis
--------------------------------------------------------
It does not write to `05_cds_qc_passed/`, the alignments, the trees or the selection
outputs. Re-filtering the cohort at this point would change the alignments, which would
invalidate every GARD partition hash, every HyPhy `.done` marker and the whole selection
arm. The cohort is pinned by accession in `19_supplements/S2_final_accessions.tsv`; this
step documents it, it does not redefine it. Any sequence it reports as FAIL is therefore a
finding to report, not a sequence that gets dropped here.

Status vocabulary
-----------------
  PASS_GOOD                  meets every criterion
  PASS_LENGTH_POLYMORPHISM   meets every correctness criterion (ATG start, valid stop,
                             in frame, no internal stop, ambiguity under threshold) and
                             fails ONLY the 0.95-1.05 length-ratio window. That window is
                             a heuristic for mis-annotation, not a correctness test, and
                             `07c_length_filter_audit.tsv` / `07b_us7_failure_breakdown.tsv`
                             establish the US7/gI cases as a real in-frame +84 nt indel.
                             Named rather than silently folded into PASS.
  PASS_RESISTANCE_VARIANT    carried over from the original log (UL23 resistance alleles)
  FAIL                       fails a correctness criterion

Outputs
-------
  07_qc_reports/07d_cds_qc_log_complete.tsv   one row per carried-forward CDS
  07_qc_reports/07e_qc_reconciliation.tsv     summary of the old log vs this one
"""

import os
import sys
import glob
import re
from datetime import datetime

import pandas as pd

DIR_BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DIR_QC = os.path.join(DIR_BASE, "07_qc_reports")
DIR_PASSED = os.path.join(DIR_BASE, "05_cds_qc_passed")
DIR_READY = os.path.join(DIR_BASE, "06_alignments_v2", "hyphy_ready")
DIR_GENOMES = os.path.join(DIR_BASE, "02_filtered_genomes")
LOG_PATH = os.path.join(DIR_QC, "qc_reconciliation.log")

STOPS = ("TAA", "TAG", "TGA")
NORM = {"HSV-1": "HSV1", "HSV-2": "HSV2"}


def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def strip_ver(acc):
    return re.sub(r"\.\d+$", "", str(acc))


def read_fasta(path):
    """Yield (header, sequence). Headers are kept whole: the tier field contains a space,
    so splitting on whitespace is exactly the bug that truncated `cohort_tier=Tier 1`."""
    header, buf = None, []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header, buf = line[1:], []
            elif header is not None:
                buf.append(line.strip())
    if header is not None:
        yield header, "".join(buf)


def assess(seq, ref_len):
    """Apply the step07 criteria. Returns (status, reasons, metrics)."""
    s = seq.replace("-", "").upper()
    n = len(s)
    starts_atg = s[:3] == "ATG"
    ends_stop = s[-3:] in STOPS
    div3 = (n % 3) == 0
    ratio = (n / ref_len) if ref_len else float("nan")
    len_ok = (0.95 <= ratio <= 1.05) if ref_len else True
    pct_n = (s.count("N") / max(1, n)) * 100.0
    ambig_ok = pct_n < 0.5

    internal = 0
    if div3:
        for i in range(0, n - 3, 3):
            if s[i:i + 3] in STOPS:
                internal += 1

    reasons = []
    if not starts_atg:
        reasons.append(f"Invalid start codon ({s[:3]})")
    if not ends_stop:
        reasons.append(f"Invalid stop codon ({s[-3:]})")
    if not div3:
        reasons.append(f"Length {n} not divisible by 3")
    if not len_ok:
        reasons.append(f"Length ratio {ratio:.3f} outside 0.95-1.05")
    if not ambig_ok:
        reasons.append(f"Ambiguity {pct_n:.2f}% >= 0.5%")
    if internal:
        reasons.append(f"Internal stop codons ({internal})")

    correctness_ok = starts_atg and ends_stop and div3 and ambig_ok and internal == 0
    if correctness_ok and len_ok:
        status = "PASS_GOOD"
    elif correctness_ok:
        status = "PASS_LENGTH_POLYMORPHISM"
    else:
        status = "FAIL"

    return status, reasons, dict(length_nt=n, ref_length_nt=ref_len,
                                 length_ratio=round(ratio, 3) if ref_len else None,
                                 starts_atg=starts_atg, ends_stop=ends_stop,
                                 divisible_by_3=div3, internal_stops=internal,
                                 pct_n=round(pct_n, 2))


def main():
    os.makedirs(DIR_QC, exist_ok=True)
    log("STEP 25: complete CDS QC record over the carried-forward set")

    # ---- reference CDS lengths, and the prior statuses we must preserve -------------
    old_path = os.path.join(DIR_QC, "07_cds_qc_log.tsv")
    if not os.path.exists(old_path):
        log(f"ABORT: {old_path} not found")
        return 2
    old = pd.read_csv(old_path, sep="\t")
    old["v"] = old.virus.map(lambda x: NORM.get(x, x))
    old["a"] = old.accession.map(strip_ver)
    old_by_key = {(r.v, r.gene, r.a): r for r in old.itertuples()}
    ref_len = {(r.v, r.gene): int(r.ref_length_nt)
               for r in old.itertuples() if pd.notna(r.ref_length_nt)}
    log(f"prior QC log: {len(old)} rows, reference lengths for {len(ref_len)} cohorts")

    # extraction log gives reference lengths for cohorts the QC log never covered
    ext_path = os.path.join(DIR_QC, "06_cds_extraction_log.tsv")
    ext_by_key = {}
    if os.path.exists(ext_path):
        ext = pd.read_csv(ext_path, sep="\t")
        ext["v"] = ext.virus.map(lambda x: NORM.get(x, x))
        ext["a"] = ext.accession.map(strip_ver)
        for r in ext.itertuples():
            ext_by_key[(r.v, r.gene, r.a)] = r
            if pd.notna(getattr(r, "ref_length", None)):
                ref_len.setdefault((r.v, r.gene), int(r.ref_length))
        log(f"extraction log: {len(ext)} rows")

    # ---- tier membership, from the genome FASTAs (authoritative) --------------------
    tier = {}
    for fname, label in (("HSV1_filtered.fasta", "Tier_1"),
                         ("HSV2_filtered.fasta", "Tier_1"),
                         ("HSV2_filtered_tier2.fasta", "Tier_2")):
        p = os.path.join(DIR_GENOMES, fname)
        if not os.path.exists(p):
            continue
        for h, _ in read_fasta(p):
            tier[strip_ver(h.split()[0])] = label
    log(f"tier assignments from 02_filtered_genomes: {len(tier)} accessions")

    # ---- the set that actually reaches the selection analysis -----------------------
    final = set()
    for p in glob.glob(os.path.join(DIR_READY, "*", "*_idmap.tsv")):
        v = os.path.basename(os.path.dirname(p))
        g = os.path.basename(p).replace("_idmap.tsv", "")
        d = pd.read_csv(p, sep="\t")
        for a in d.accession:
            final.add((v, g, strip_ver(str(a).split("|")[0])))
    log(f"sequences in the final analysis cohorts: {len(final)}")

    # ---- assess every carried-forward sequence --------------------------------------
    rows = []
    for path in sorted(glob.glob(os.path.join(DIR_PASSED, "*", "*.fasta"))):
        base = os.path.basename(path)
        if base.endswith("_dedup.fasta") or "resistance" in base:
            continue
        virus = os.path.basename(os.path.dirname(path))
        gene = base[:-len(".fasta")]
        for header, seq in read_fasta(path):
            acc = strip_ver(header.split("|")[0])
            key = (virus, gene, acc)
            status, reasons, m = assess(seq, ref_len.get((virus, gene)))

            prior = old_by_key.get(key)
            prior_status = prior.qc_status if prior is not None else "NOT_IN_PRIOR_LOG"
            # The resistance-variant exemption was a curation decision; keep it.
            if prior_status == "PASS_RESISTANCE_VARIANT":
                status = "PASS_RESISTANCE_VARIANT"

            mt = re.search(r"cohort_tier=([^|]+)", header)
            rows.append(dict(
                virus=virus, gene=gene, accession=acc,
                cohort_tier=tier.get(acc, "UNKNOWN"),
                cohort_tier_in_header=(mt.group(1).strip().replace(" ", "_")
                                       if mt else "ABSENT"),
                in_final_analysis=key in final,
                qc_status=status,
                prior_qc_status=prior_status,
                status_changed=(status != prior_status),
                failure_reasons="; ".join(reasons),
                fasta_header=header,
                **m))

    df = pd.DataFrame(rows)
    out = os.path.join(DIR_QC, "07d_cds_qc_log_complete.tsv")
    df.to_csv(out, sep="\t", index=False)
    log(f"wrote 07d_cds_qc_log_complete.tsv ({len(df)} rows)")

    # ---- summary --------------------------------------------------------------------
    log("")
    log("status over the carried-forward set:")
    for k, v in df.qc_status.value_counts().items():
        log(f"    {k:28} {v}")
    log("")
    log("coverage by the prior log:")
    for k, v in df.prior_qc_status.value_counts().items():
        log(f"    {k:28} {v}")

    fails = df[df.qc_status == "FAIL"]
    lp = df[df.qc_status == "PASS_LENGTH_POLYMORPHISM"]
    log("")
    log(f"PASS_LENGTH_POLYMORPHISM: {len(lp)} sequences, "
        f"{int(lp.in_final_analysis.sum())} of them in the final analysis")
    if len(lp):
        for (v, g), grp in lp.groupby(["virus", "gene"]):
            log(f"    {v} {g}: {len(grp)} "
                f"(ratios {grp.length_ratio.min():.3f}-{grp.length_ratio.max():.3f}, "
                f"{int(grp.in_final_analysis.sum())} analysed)")
    log("")
    log(f"FAIL on a correctness criterion: {len(fails)} sequences, "
        f"{int(fails.in_final_analysis.sum())} of them in the final analysis")
    if int(fails.in_final_analysis.sum()):
        log("    *** these are in the analysed cohort and must be reported as a defect ***")
        for r in fails[fails.in_final_analysis].itertuples():
            log(f"    {r.virus} {r.gene} {r.accession}: {r.failure_reasons}")

    log("")
    log("tier composition of the final analysis cohorts:")
    fa = df[df.in_final_analysis]
    for v, grp in fa.groupby("virus"):
        counts = grp.cohort_tier.value_counts().to_dict()
        n_acc = grp.accession.nunique()
        t2 = grp[grp.cohort_tier == "Tier_2"].accession.nunique()
        log(f"    {v}: {len(grp)} sequences, {n_acc} accessions, "
            f"{t2} of them Tier_2 ({100.0 * t2 / max(1, n_acc):.1f}%)  {counts}")

    hdr_bad = int((fa.cohort_tier_in_header == "Tier").sum()) if len(fa) else 0
    log("")
    log(f"headers whose cohort_tier field is truncated to 'Tier': {hdr_bad} "
        f"(the field is written with a space, so ID-based readers cut it; the "
        f"cohort_tier column of this file is derived from 02_filtered_genomes instead)")

    rec = pd.DataFrame([
        dict(quantity="CDS rows in prior QC log (07_cds_qc_log.tsv)", n=len(old)),
        dict(quantity="CDS carried forward (05_cds_qc_passed/)", n=len(df)),
        dict(quantity="carried forward with no row in prior log",
             n=int((df.prior_qc_status == "NOT_IN_PRIOR_LOG").sum())),
        dict(quantity="carried forward marked FAIL by prior log",
             n=int((df.prior_qc_status == "FAIL").sum())),
        dict(quantity="PASS_GOOD (this log)",
             n=int((df.qc_status == "PASS_GOOD").sum())),
        dict(quantity="PASS_LENGTH_POLYMORPHISM (this log)", n=len(lp)),
        dict(quantity="PASS_RESISTANCE_VARIANT (this log)",
             n=int((df.qc_status == "PASS_RESISTANCE_VARIANT").sum())),
        dict(quantity="FAIL on a correctness criterion (this log)", n=len(fails)),
        dict(quantity="FAIL and in the final analysis",
             n=int(fails.in_final_analysis.sum())),
        dict(quantity="sequences in the final analysis cohorts",
             n=int(df.in_final_analysis.sum())),
    ])
    rec_out = os.path.join(DIR_QC, "07e_qc_reconciliation.tsv")
    rec.to_csv(rec_out, sep="\t", index=False)
    log(f"wrote 07e_qc_reconciliation.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
