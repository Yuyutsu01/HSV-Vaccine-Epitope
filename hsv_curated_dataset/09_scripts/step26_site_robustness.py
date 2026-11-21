"""
Step 26: robustness assessment of every positively selected site (verification pass, 2026-11-21).

Why this exists
---------------
The site-level screen reports 92 sites as positively selected by at least one method and 2
by two methods. Nothing in steps 17-18 asks the question a reviewer asks first: **how many
sequences actually carry the variant that produced the signal?**

On alignments this conserved, a single non-synonymous singleton at an otherwise invariant
codon drives beta/alpha towards infinity, and MEME in particular will fire with an
essentially unbounded beta. Such a call is a property of the likelihood surface at a
near-invariant site, not evidence of adaptive evolution. The two concordant sites turn out
to rest on 4 and 2 isolates respectively, with FUBAR posteriors of 0.956 and 0.952 against
a 0.95 threshold, and FEL supporting neither.

This step therefore annotates every reported site with:
  - how many sequences and how many weighted isolates carry a non-majority codon
  - how many distinct non-majority codons there are
  - the minimum number of nucleotide substitutions separating the majority codon from each
    variant (a 2- or 3-nt change in a single isolate is more consistent with a sequencing
    or assembly artifact than with a point mutation)
  - whether any supporting sequence fails CDS QC in 07d_cds_qc_log_complete.tsv
  - whether the supporting accessions are adjacent in accession order, which flags
    same-submission batch artifacts
  - the resulting robustness verdict

Verdict rule, declared here rather than tuned
---------------------------------------------
  NOT_ROBUST   fewer than MIN_ISOLATES isolates carry the variant, or every variant
               observation comes from a QC-failing sequence
  WEAK         MIN_ISOLATES or more isolates, but either only one distinct variant codon
               supported by a single sequence, or any variant is a multi-nucleotide change
               seen in only one isolate
  ROBUST       otherwise

MIN_ISOLATES is 5, i.e. a site must show the variant in at least 5 of the thousands of
isolates represented before the signal is treated as biological. This is a deliberately
blunt threshold and it is applied uniformly; it is not chosen per site.

This step is read-only with respect to the analysis. It annotates, it does not re-filter.
"""

import os
import re
import sys
import glob
from datetime import datetime
from collections import Counter, defaultdict

import pandas as pd

DIR_BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DIR_SEL = os.path.join(DIR_BASE, "14_selection")
DIR_READY = os.path.join(DIR_BASE, "06_alignments_v2", "hyphy_ready")
DIR_PARTS = os.path.join(DIR_SEL, "partitions")
DIR_QC = os.path.join(DIR_BASE, "07_qc_reports")
LOG_PATH = os.path.join(DIR_SEL, "site_robustness.log")

MIN_ISOLATES = 5
STOPS = ("TAA", "TAG", "TGA")


def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def read_fasta(path):
    seqs, key, buf = {}, None, []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if key is not None:
                    seqs[key] = "".join(buf)
                key, buf = line[1:].split()[0], []
            elif key is not None:
                buf.append(line)
    if key is not None:
        seqs[key] = "".join(buf)
    return seqs


def nt_diff(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def acc_number(acc):
    m = re.search(r"(\d+)", str(acc))
    return int(m.group(1)) if m else None


def main():
    log("STEP 26: robustness assessment of reported positively selected sites")

    conc_path = os.path.join(DIR_SEL, "selection_concordance.tsv")
    if not os.path.exists(conc_path):
        log(f"ABORT: {conc_path} not found; run step18 first")
        return 2
    conc = pd.read_csv(conc_path, sep="\t")
    site = pd.read_csv(os.path.join(DIR_SEL, "site_level_results.tsv"), sep="\t")
    log(f"sites to assess: {len(conc)}")

    qc = None
    qc_path = os.path.join(DIR_QC, "07d_cds_qc_log_complete.tsv")
    if os.path.exists(qc_path):
        qc = pd.read_csv(qc_path, sep="\t")
        qc["acc0"] = qc.accession.map(lambda a: re.sub(r"\.\d+$", "", str(a)))
        qc_bad = {(r.virus, r.gene, r.acc0) for r in qc.itertuples()
                  if r.qc_status == "FAIL"}
        log(f"QC log: {len(qc)} rows, {len(qc_bad)} FAIL keys")
    else:
        qc_bad = set()
        log("WARNING: 07d_cds_qc_log_complete.tsv missing; run step25 for the QC column")

    # cache alignments, partition alignments and idmaps
    aln_cache, idmap_cache = {}, {}

    def get_aln(cohort, gene, partition):
        key = (cohort, gene, partition)
        if key in aln_cache:
            return aln_cache[key]
        if partition and partition > 0:
            p = os.path.join(DIR_PARTS, cohort, f"{gene}_part{partition}.fasta")
        else:
            p = os.path.join(DIR_READY, cohort, f"{gene}.fasta")
        aln_cache[key] = read_fasta(p) if os.path.exists(p) else {}
        return aln_cache[key]

    def get_idmap(cohort, gene):
        key = (cohort, gene)
        if key in idmap_cache:
            return idmap_cache[key]
        p = os.path.join(DIR_READY, cohort, f"{gene}_idmap.tsv")
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t")
            idmap_cache[key] = (dict(zip(d.safe_id, d.accession)),
                                dict(zip(d.safe_id, d.cluster_size)))
        else:
            idmap_cache[key] = ({}, {})
        return idmap_cache[key]

    rows = []
    for r in conc.itertuples():
        gene, cohort, ref_res = r.gene, r.cohort, r.ref_residue

        # the per-method rows give the partition and the column within that unit
        sr = site[(site.gene == gene) & (site.cohort == cohort)
                  & (site.ref_residue == ref_res) & site.positively_selected.astype(bool)]
        if sr.empty:
            continue
        partition = int(sr.partition.iloc[0])
        col = int(sr.site_in_unit.iloc[0])

        seqs = get_aln(cohort, gene, partition)
        acc_of, size_of = get_idmap(cohort, gene)
        if not seqs:
            rows.append(dict(gene=gene, cohort=cohort, ref_residue=ref_res,
                             verdict="UNASSESSED", note="alignment not found"))
            continue

        i = (col - 1) * 3
        codons = {sid: s[i:i + 3] for sid, s in seqs.items() if len(s) >= i + 3}
        counts = Counter(c for c in codons.values() if "-" not in c and "N" not in c)
        if not counts:
            rows.append(dict(gene=gene, cohort=cohort, ref_residue=ref_res,
                             verdict="UNASSESSED", note="no resolved codon"))
            continue
        major = counts.most_common(1)[0][0]

        variants = defaultdict(list)
        for sid, c in codons.items():
            if c != major and "-" not in c and "N" not in c:
                variants[c].append(sid)

        var_seqs = sum(len(v) for v in variants.values())
        var_isolates = sum(size_of.get(sid, 1) for v in variants.values() for sid in v)
        total_isolates = sum(size_of.get(sid, 1) for sid in codons)
        n_variant_codons = len(variants)

        max_ntdiff = max((nt_diff(major, c) for c in variants), default=0)
        multi_nt_singleton = any(
            nt_diff(major, c) >= 2 and sum(size_of.get(s, 1) for s in v) <= 1
            for c, v in variants.items())

        supporting_acc = []
        any_qc_fail = False
        all_qc_fail = bool(variants)
        for c, v in variants.items():
            for sid in v:
                acc = re.sub(r"\.\d+$", "", str(acc_of.get(sid, "")).split("|")[0])
                supporting_acc.append(acc)
                if (cohort, gene, acc) in qc_bad:
                    any_qc_fail = True
                else:
                    all_qc_fail = False

        nums = sorted(n for n in (acc_number(a) for a in supporting_acc) if n)
        adjacent_batch = len(nums) >= 2 and any(
            b - a <= 5 for a, b in zip(nums, nums[1:]))

        if var_isolates < MIN_ISOLATES or (variants and all_qc_fail):
            verdict = "NOT_ROBUST"
        elif (n_variant_codons == 1 and var_seqs == 1) or multi_nt_singleton:
            verdict = "WEAK"
        else:
            verdict = "ROBUST"

        rows.append(dict(
            gene=gene, cohort=cohort, ref_residue=ref_res,
            methods=getattr(r, "methods", ""),
            n_methods=getattr(r, "n_methods", ""),
            partition=partition, codon_col_in_unit=col,
            majority_codon=major, majority_n_seqs=counts[major],
            n_variant_codons=n_variant_codons,
            variant_codons=";".join(f"{c}x{len(v)}" for c, v in sorted(variants.items())),
            variant_seqs=var_seqs, variant_isolates=var_isolates,
            total_isolates=total_isolates,
            variant_isolate_pct=round(100.0 * var_isolates / max(1, total_isolates), 3),
            max_nt_changes=max_ntdiff,
            multi_nt_singleton=multi_nt_singleton,
            any_supporting_seq_qc_fail=any_qc_fail,
            all_supporting_seqs_qc_fail=all_qc_fail if variants else False,
            supporting_accessions=";".join(sorted(set(supporting_acc))),
            adjacent_accession_batch=adjacent_batch,
            verdict=verdict, note=""))

    df = pd.DataFrame(rows)
    out = os.path.join(DIR_SEL, "site_robustness.tsv")
    df.to_csv(out, sep="\t", index=False)
    log(f"wrote site_robustness.tsv ({len(df)} rows)")

    log("")
    log(f"robustness verdicts (MIN_ISOLATES = {MIN_ISOLATES}):")
    for k, v in df.verdict.value_counts().items():
        log(f"    {k:14} {v}")

    log("")
    log("the two sites that satisfy the >=2-method concordance rule:")
    two = df[df.n_methods.astype(str) == "2"] if "n_methods" in df.columns else df.iloc[:0]
    for r in two.itertuples():
        log(f"    {r.cohort} {r.gene} residue {r.ref_residue:.0f}: {r.verdict}")
        log(f"        majority {r.majority_codon} x{r.majority_n_seqs}, "
            f"variants {r.variant_codons}")
        log(f"        {r.variant_isolates} of {r.total_isolates} isolates "
            f"({r.variant_isolate_pct}%) carry a variant; "
            f"max {r.max_nt_changes} nt change(s)")
        log(f"        supporting accessions: {r.supporting_accessions}"
            f"{'  [adjacent batch]' if r.adjacent_accession_batch else ''}"
            f"{'  [includes QC-FAIL sequence]' if r.any_supporting_seq_qc_fail else ''}")

    robust = df[df.verdict == "ROBUST"]
    log("")
    if robust.empty:
        log("NO reported site is ROBUST under this rule.")
        log("The site-level screen therefore provides no site that survives a minimum-")
        log("variant-support check, and the headline must say so.")
    else:
        log(f"{len(robust)} site(s) ROBUST:")
        for r in robust.itertuples():
            log(f"    {r.cohort} {r.gene} residue {r.ref_residue:.0f} "
                f"({r.methods}): {r.variant_isolates}/{r.total_isolates} isolates")

    log("")
    log("distribution of variant isolate support across all reported sites:")
    if "variant_isolates" in df.columns and df.variant_isolates.notna().any():
        vi = df.variant_isolates.dropna()
        log(f"    min {int(vi.min())}  median {int(vi.median())}  max {int(vi.max())}")
        log(f"    sites with 1 isolate: {int((vi == 1).sum())}; "
            f"<5 isolates: {int((vi < MIN_ISOLATES).sum())} of {len(vi)}")
    log("")
    log("Reported alongside the primary tables, not instead of them: "
        "selected_sites_final.tsv is unchanged and this file qualifies it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
