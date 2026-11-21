"""
Step 18 (CORRECTION A4): Assemble HyPhy site results, apply FDR, require method concordance.

The previous pipeline published raw p-values and even a p<0.10 column across 198 + 62 runs
testing 91-1,318 codons each, with no multiple-testing correction anywhere in the repository
(CORRECTIONS_AND_PENDING_WORK.md, A4). This step fixes that:

  - FEL and MEME site p-values: Benjamini-Hochberg FDR WITHIN each analysis unit (one family
    per gene x cohort x partition x method), reported as q-values.
  - BUSTED gene-wide p-values: BH across the whole gene x cohort family.
  - aBSREL: HyPhy's own Holm-Bonferroni corrected branch p-values are used as reported.
  - FUBAR: posterior probability, reported at >= 0.95 (the >= 0.90 column the old summaries
    headlined is retained only as a secondary column).
  - A site is called positively selected only when it is significant at q < 0.05 in at least
    two independent methods, as the proposal specified.

Every site is reported in reference residue numbering, by composing:
  partition offset (step16) -> masked alignment codon column -> reference residue (step09c map).

Outputs:
  - 14_selection/site_level_results.tsv        every site x method, with q-values
  - 14_selection/gene_level_results.tsv        BUSTED, BH-corrected
  - 14_selection/branch_level_results.tsv      aBSREL selected branches
  - 14_selection/selection_concordance.tsv     per site: methods supporting selection
  - 14_selection/selected_sites_final.tsv      q<0.05 in >= 2 methods (the usable call set)
  - 14_selection/fdr.log
"""

import os
import re
import sys
import json
import glob
import argparse
from datetime import datetime

import pandas as pd
from statsmodels.stats.multitest import multipletests

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_SEL = os.path.join(DIR_BASE, "14_selection")
DIR_HYPHY = os.path.join(DIR_SEL, "hyphy")
DIR_READY = os.path.join(DIR_BASE, "06_alignments_v2", "hyphy_ready")
DIR_MAPS = os.path.join(DIR_BASE, "06_alignments_v2", "coordinate_maps")
DIR_CONS = os.path.join(DIR_BASE, "11_conservation")

Q_CUT = 0.05
FUBAR_PP = 0.95
FUBAR_PP_SECONDARY = 0.90
MIN_METHODS = 2

LOG_PATH = os.path.join(DIR_SEL, "fdr.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_SEL, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def load_json(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception as exc:
        log(f"  unreadable {os.path.basename(path)}: {type(exc).__name__}")
        return None


def mle_table(js):
    """HyPhy MLE block -> (list of header names, list of row lists)."""
    mle = js.get("MLE") or {}
    headers = [h[0] if isinstance(h, list) else str(h) for h in (mle.get("headers") or [])]
    content = mle.get("content") or {}
    rows = []
    for key in sorted(content.keys(), key=lambda k: int(k) if str(k).isdigit() else 0):
        block = content[key]
        if isinstance(block, list):
            rows.extend(block)
    return headers, rows


def col_index(headers, *candidates):
    low = [h.lower() for h in headers]
    for cand in candidates:
        c = cand.lower()
        for i, h in enumerate(low):
            if h == c:
                return i
    for cand in candidates:
        c = cand.lower()
        for i, h in enumerate(low):
            if c in h:
                return i
    return None


def partition_offsets():
    """{(gene, cohort, partition): start codon offset (0-based)} from step16 ranges."""
    path = os.path.join(DIR_SEL, "gard_summary_per_type.tsv")
    out = {}
    if not os.path.exists(path):
        return out
    df = pd.read_csv(path, sep="\t")
    for r in df.itertuples():
        ranges = str(getattr(r, "partition_ranges", "None") or "None")
        if ranges in ("None", "nan"):
            continue
        for n, rng in enumerate(ranges.split(";"), start=1):
            if "-" not in rng:
                continue
            start_nt = int(rng.split("-")[0])
            out[(r.gene, r.cohort, n)] = (start_nt - 1) // 3
    return out


def colmap_index():
    """{(gene, cohort): {masked_codon_col: ref_residue}}"""
    out = {}
    for cohort in ("HSV1", "HSV2"):
        for path in glob.glob(os.path.join(DIR_MAPS, cohort, "*_colmap.tsv")):
            gene = os.path.basename(path).replace("_colmap.tsv", "")
            df = pd.read_csv(path, sep="\t")
            sub = df[df.retained_after_masking.astype(str).str.lower() == "true"]
            out[(gene, cohort)] = {
                int(m): (int(r) if pd.notna(r) else None)
                for m, r in zip(sub.masked_codon_col, sub.ref_residue) if pd.notna(m)}
    return out


def conservation_index():
    path = os.path.join(DIR_CONS, "per_residue_conservation.tsv")
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    idx = {}
    for vt, rescol, idcol in (("HSV1", "hsv1_ref_residue", "hsv1_pct_identity_isolate"),
                             ("HSV2", "hsv2_ref_residue", "hsv2_pct_identity_isolate")):
        sub = df[df[rescol].notna()]
        for gene, res, ident, both in zip(sub.gene, sub[rescol], sub[idcol],
                                          sub.conserved_both_types):
            idx[(gene, vt, int(res))] = (float(ident), bool(both))
    return idx


def parse_unit(path):
    """One HyPhy output -> (kind, rows). kind in {site, gene, branch}."""
    name = os.path.basename(path)
    m = re.match(r"(.+)_(BUSTED|FEL|MEME|FUBAR|ABSREL)\.json$", name, re.I)
    if not m:
        return None, []
    unit, method = m.group(1), m.group(2).upper()
    cohort = os.path.basename(os.path.dirname(path))
    pm = re.match(r"(.+)_part(\d+)$", unit)
    gene, partition = (pm.group(1), int(pm.group(2))) if pm else (unit, 0)

    js = load_json(path)
    if js is None:
        return None, []

    base = dict(unit=unit, gene=gene, cohort=cohort, partition=partition, method=method)

    if method == "BUSTED":
        tr = js.get("test results") or {}
        return "gene", [dict(**base, p_value=tr.get("p-value"),
                             lrt=tr.get("LRT"),
                             evidence_ratio_tested=None)]

    if method == "ABSREL":
        rows = []
        ba = (js.get("branch attributes") or {}).get("0") or {}
        for branch, attrs in ba.items():
            if not isinstance(attrs, dict):
                continue
            p = attrs.get("Corrected P-value", attrs.get("p-value"))
            if p is None:
                continue
            rows.append(dict(**base, branch=branch, corrected_p_value=p,
                             uncorrected_p_value=attrs.get("Uncorrected P-value"),
                             rate_classes=attrs.get("Rate classes"),
                             lrt=attrs.get("LRT")))
        return "branch", rows

    headers, table = mle_table(js)
    if not table:
        return None, []

    rows = []
    if method == "FUBAR":
        i_a = col_index(headers, "alpha")
        i_b = col_index(headers, "beta")
        i_pos = col_index(headers, "Prob[alpha<beta]", "prob[alpha<beta]")
        i_neg = col_index(headers, "Prob[alpha>beta]", "prob[alpha>beta]")
        for n, row in enumerate(table, start=1):
            rows.append(dict(**base, site_in_unit=n,
                             alpha=row[i_a] if i_a is not None else None,
                             beta=row[i_b] if i_b is not None else None,
                             posterior_positive=row[i_pos] if i_pos is not None else None,
                             posterior_negative=row[i_neg] if i_neg is not None else None,
                             p_value=None))
        return "site", rows

    # FEL and MEME
    i_p = col_index(headers, "p-value", "p_value")
    i_a = col_index(headers, "alpha")
    i_b = col_index(headers, "beta", "&beta;<sup>+</sup>")
    for n, row in enumerate(table, start=1):
        p = row[i_p] if i_p is not None and i_p < len(row) else None
        a = row[i_a] if i_a is not None and i_a < len(row) else None
        b = row[i_b] if i_b is not None and i_b < len(row) else None
        rows.append(dict(**base, site_in_unit=n, alpha=a, beta=b, p_value=p,
                         posterior_positive=None, posterior_negative=None))
    return "site", rows


def main():
    ap = argparse.ArgumentParser()
    # These tables are the end of the selection arm and feed the candidate shortlist.
    # Built from a partial selection run they would silently omit whole genes, with
    # nothing in the output saying which, so completeness is required by default.
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="assemble from whatever selection output exists (labels the "
                         "result provisional; for validating the parsers only)")
    args = ap.parse_args()

    log("=" * 78)
    log("STEP 18: FDR correction and cross-method concordance")

    files = sorted(glob.glob(os.path.join(DIR_HYPHY, "*", "*.json")))
    if not files:
        log(f"No HyPhy outputs under {DIR_HYPHY}. Run step16 then step17 first.")
        return 1
    log(f"Found {len(files)} HyPhy output file(s)")

    # Every gene x cohort with a HyPhy-ready alignment must be represented.
    expected = {(os.path.basename(f).replace(".fasta", ""),
                 os.path.basename(os.path.dirname(f)))
                for f in glob.glob(os.path.join(DIR_READY, "*", "*.fasta"))
                if "best-gard" not in f}
    seen = {(re.sub(r"_part\d+$", "", os.path.basename(f).rsplit("_", 1)[0]),
             os.path.basename(os.path.dirname(f))) for f in files}
    absent = sorted(expected - seen)
    if absent and not args.allow_incomplete:
        log(f"ABORT: selection output is missing {len(absent)} of {len(expected)} "
            f"gene x cohort unit(s): " + ", ".join(f"{g}/{c}" for g, c in absent[:12])
            + (" ..." if len(absent) > 12 else ""))
        log("Finish step17 for those units, or pass --allow-incomplete to assemble a "
            "provisional set for checking the parsers.")
        return 2
    # The output filenames do not change between a provisional and a final assembly, so
    # leave a marker on disk saying which this is; a complete run removes it.
    flag = os.path.join(DIR_SEL, "PROVISIONAL_RESULTS.txt")
    if absent:
        log(f"WARNING: --allow-incomplete set. PROVISIONAL RESULT -- {len(absent)} unit(s) "
            f"absent and NOT represented in any output table: "
            + ", ".join(f"{g}/{c}" for g, c in absent))
        with open(flag, "w", encoding="utf-8") as fh:
            fh.write("The *_results.tsv, selection_concordance.tsv and "
                     "selected_sites_final.tsv files in this directory are PROVISIONAL.\n"
                     f"Assembled {datetime.now():%Y-%m-%d %H:%M} with these gene x cohort "
                     "unit(s) absent from every table:\n")
            for g, c in absent:
                fh.write(f"  {g}/{c}\n")
            fh.write("\nRe-run step18 with no flags once step17 covers them; that run "
                     "deletes this file.\n")
    elif os.path.exists(flag):
        os.remove(flag)
        log("Selection output is complete; removed the PROVISIONAL_RESULTS marker.")

    site_rows, gene_rows, branch_rows = [], [], []
    for path in files:
        kind, rows = parse_unit(path)
        if kind == "site":
            site_rows.extend(rows)
        elif kind == "gene":
            gene_rows.extend(rows)
        elif kind == "branch":
            branch_rows.extend(rows)

    offsets = partition_offsets()
    colmaps = colmap_index()
    cons = conservation_index()

    sites = pd.DataFrame(site_rows)
    if len(sites):
        # Partition-local site index -> masked alignment codon column -> reference residue.
        off = [offsets.get((g, c, p), 0) if p else 0
               for g, c, p in zip(sites.gene, sites.cohort, sites.partition)]
        sites["alignment_codon_col"] = [o + n for o, n in zip(off, sites.site_in_unit)]
        sites["ref_residue"] = [
            colmaps.get((g, c), {}).get(int(col))
            for g, c, col in zip(sites.gene, sites.cohort, sites.alignment_codon_col)]

        ident, both = [], []
        for g, c, r in zip(sites.gene, sites.cohort, sites.ref_residue):
            hit = cons.get((g, c, int(r))) if pd.notna(r) else None
            ident.append(hit[0] if hit else float("nan"))
            both.append(bool(hit[1]) if hit else False)
        sites["within_type_identity"] = ident
        sites["conserved_both_types"] = both

        # Benjamini-Hochberg within each gene x cohort x partition x method family.
        sites["q_value"] = float("nan")
        for key, grp in sites.groupby(["gene", "cohort", "partition", "method"]):
            if key[3] not in ("FEL", "MEME"):
                continue
            mask = grp.p_value.notna()
            if not mask.any():
                continue
            pv = grp.loc[mask, "p_value"].astype(float).clip(0, 1)
            _rej, q, _a, _b = multipletests(pv, alpha=Q_CUT, method="fdr_bh")
            sites.loc[grp.index[mask], "q_value"] = q
            log(f"  BH {key[0]} {key[1]} part{key[2]} {key[3]}: "
                f"{len(pv)} sites, {(q < Q_CUT).sum()} at q<{Q_CUT} "
                f"(raw p<0.05 would have been {(pv < 0.05).sum()})")

        sites["significant"] = False
        fel_meme = sites.method.isin(["FEL", "MEME"])
        sites.loc[fel_meme, "significant"] = sites.loc[fel_meme, "q_value"] < Q_CUT
        fub = sites.method == "FUBAR"
        sites.loc[fub, "significant"] = (
            sites.loc[fub, "posterior_positive"].astype(float) >= FUBAR_PP)
        sites["significant_secondary"] = False
        sites.loc[fub, "significant_secondary"] = (
            sites.loc[fub, "posterior_positive"].astype(float) >= FUBAR_PP_SECONDARY)

        # FEL distinguishes purifying from diversifying by the sign of beta - alpha.
        def direction(r):
            if r.method == "FUBAR":
                return "diversifying" if (r.posterior_positive or 0) >= FUBAR_PP else (
                    "purifying" if (r.posterior_negative or 0) >= FUBAR_PP else "neutral")
            try:
                if r.beta is None or r.alpha is None:
                    return "unknown"
                return "diversifying" if float(r.beta) > float(r.alpha) else "purifying"
            except (TypeError, ValueError):
                return "unknown"
        sites["direction"] = [direction(r) for r in sites.itertuples()]
        sites["positively_selected"] = sites.significant & (sites.direction == "diversifying")

        sites.to_csv(os.path.join(DIR_SEL, "site_level_results.tsv"), sep="\t", index=False)
        log(f"Wrote site_level_results.tsv ({len(sites)} site x method rows)")

    genes = pd.DataFrame(gene_rows)
    if len(genes):
        mask = genes.p_value.notna()
        genes["q_value"] = float("nan")
        if mask.any():
            pv = genes.loc[mask, "p_value"].astype(float).clip(0, 1)
            _rej, q, _a, _b = multipletests(pv, alpha=Q_CUT, method="fdr_bh")
            genes.loc[mask, "q_value"] = q
        genes["significant"] = genes.q_value < Q_CUT
        genes.to_csv(os.path.join(DIR_SEL, "gene_level_results.tsv"), sep="\t", index=False)
        log(f"Wrote gene_level_results.tsv ({len(genes)} BUSTED runs, "
            f"{int(genes.significant.sum())} at q<{Q_CUT}; "
            f"raw p<0.05 would have been {int((genes.p_value.astype(float) < 0.05).sum())})")

    branches = pd.DataFrame(branch_rows)
    if len(branches):
        branches["significant"] = branches.corrected_p_value.astype(float) < Q_CUT
        branches.to_csv(os.path.join(DIR_SEL, "branch_level_results.tsv"),
                        sep="\t", index=False)
        sel = branches[branches.significant]
        log(f"Wrote branch_level_results.tsv ({len(branches)} branches, "
            f"{len(sel)} selected at Holm-corrected p<{Q_CUT})")

    if len(sites):
        pos = sites[sites.positively_selected & sites.ref_residue.notna()]
        conc = (pos.groupby(["gene", "cohort", "ref_residue"])
                   .agg(methods=("method", lambda s: ";".join(sorted(set(s)))),
                        n_methods=("method", lambda s: len(set(s))),
                        min_q=("q_value", "min"),
                        within_type_identity=("within_type_identity", "first"),
                        conserved_both_types=("conserved_both_types", "first"))
                   .reset_index()
                   .sort_values(["gene", "cohort", "ref_residue"]))
        conc.to_csv(os.path.join(DIR_SEL, "selection_concordance.tsv"), sep="\t", index=False)

        final = conc[conc.n_methods >= MIN_METHODS]
        final.to_csv(os.path.join(DIR_SEL, "selected_sites_final.tsv"), sep="\t", index=False)
        log(f"Wrote selection_concordance.tsv ({len(conc)} sites positively selected in "
            f">=1 method) and selected_sites_final.tsv ({len(final)} sites supported by "
            f">={MIN_METHODS} methods)")
        for r in final.itertuples():
            log(f"    {r.cohort} {r.gene} residue {int(r.ref_residue)}: "
                f"{r.methods} (min q={r.min_q:.4g})")
        if final.empty:
            log("    No site met the concordance requirement. For a virus as conserved as "
                "HSV this is an expected outcome, not a pipeline failure.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
