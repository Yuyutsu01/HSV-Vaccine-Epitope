"""
Step 19 (PENDING WORK B5 + B7): Benchmark epitope prediction against curated IEDB epitopes.

Proposal ML section ("benchmarking the predictors") and pipeline step 9 ("validation against
known biology"). This is a small evaluation study, not a trained model.

Honest labelling, as the proposal's curation step 10 demands:
  - POSITIVES: peptides from the reference proteins that exactly match, or fall inside, a
    curated human-host, positively assayed IEDB epitope for that gene and type (step14).
  - BACKGROUND: every other enumerated peptide from the same proteins. These are NOT true
    negatives -- most have simply never been tested. Metrics are therefore reported as
    "positives versus untested background", and the headline number is the gene-grouped
    mean, so one heavily studied antigen (gD, gB) cannot carry the score.

Evaluated scores:
  - MHCflurry presentation score and affinity (step15)
  - conservation alone (within-type isolate identity), as a null-ish baseline
  - a combined rank, to show whether conservation adds anything to binding

Metrics: ROC AUC, average precision, precision@k, and enrichment of IEDB positives in the
top decile, computed per gene x type and then averaged across genes (grouped evaluation).

Outputs:
  - 15_benchmark/benchmark_per_gene.tsv
  - 15_benchmark/benchmark_summary.tsv
  - 15_benchmark/labelled_peptides.tsv
  - 15_benchmark/benchmark.log
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_EPI = os.path.join(DIR_BASE, "13_epitopes")
DIR_IEDB = os.path.join(DIR_BASE, "12_iedb")
DIR_OUT = os.path.join(DIR_BASE, "15_benchmark")

MIN_POSITIVES = 5          # a gene needs this many positives to be scored
MIN_EPITOPE_LEN = 8        # class I ligand length range; longer IEDB records are
MAX_EPITOPE_LEN = 14       # peptide-scan constructs, not single class I epitopes
TOP_K = (20, 50, 100)
LOG_PATH = os.path.join(DIR_OUT, "benchmark.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def label_peptides(pred, iedb):
    """
    A peptide is positive when it lies inside a curated epitope's reference range for the
    same gene and type, or when its sequence is a substring of that epitope.
    """
    pred = pred.copy()
    pred["is_iedb_positive"] = False
    pred["iedb_support"] = ""

    for (vt, gene), grp in iedb.groupby(["virus_type", "gene"]):
        mask = (pred.virus_type == vt) & (pred.gene == gene)
        if not mask.any():
            continue
        sub = pred[mask]
        hit = np.zeros(len(sub), dtype=bool)
        kinds = [set() for _ in range(len(sub))]

        starts = sub.ref_start.to_numpy()
        ends = sub.ref_end.to_numpy()
        peps = sub.peptide.to_numpy()

        for row in grp.itertuples():
            inside = (starts >= row.ref_start) & (ends <= row.ref_end)
            substr = np.array([p in row.peptide for p in peps])
            found = inside | substr
            hit |= found
            for i in np.nonzero(found)[0]:
                kinds[i].add(row.epitope_class)

        pred.loc[sub.index, "is_iedb_positive"] = hit
        pred.loc[sub.index, "iedb_support"] = [";".join(sorted(k)) for k in kinds]
    return pred


def precision_at_k(y, score, k):
    if len(y) <= k:
        return float("nan")
    order = np.argsort(-score)
    return float(np.mean(y[order[:k]]))


def evaluate(df, score_col, higher_is_better=True):
    y = df.is_iedb_positive.to_numpy().astype(int)
    s = df[score_col].to_numpy(dtype=float)
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    if y.sum() < MIN_POSITIVES or y.sum() == len(y):
        return None
    if not higher_is_better:
        s = -s
    out = dict(n_peptides=len(y), n_positives=int(y.sum()),
               positive_rate=round(float(y.mean()), 4),
               auc=round(float(roc_auc_score(y, s)), 4),
               average_precision=round(float(average_precision_score(y, s)), 4))
    for k in TOP_K:
        out[f"precision_at_{k}"] = round(precision_at_k(y, s, k), 4)
    cut = np.quantile(s, 0.9)
    top = y[s >= cut]
    out["top_decile_enrichment"] = (round(float(top.mean() / y.mean()), 3)
                                    if len(top) and y.mean() > 0 else float("nan"))
    return out


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 19: predictor benchmarking against curated IEDB epitopes")

    pred_path = os.path.join(DIR_EPI, "class1_peptide_predictions.tsv")
    iedb_path = os.path.join(DIR_IEDB, "iedb_curated_epitopes.tsv")
    for p in (pred_path, iedb_path):
        if not os.path.exists(p):
            log(f"ABORT: missing {p}")
            return 1

    pred = pd.read_csv(pred_path, sep="\t")
    iedb = pd.read_csv(iedb_path, sep="\t")
    log(f"Loaded {len(pred)} predicted peptides and {len(iedb)} curated IEDB mappings")

    # T-cell relevant evidence only, for a class I predictor, and only epitopes of a length
    # a class I molecule actually presents. Without the length ceiling, the long overlapping
    # peptide scans IEDB holds for gD and gB (20-35 aa) blanket those proteins, so nearly
    # every enumerated 9-mer falls inside a "positive" range and the AUC collapses toward
    # chance for exactly the best-studied antigens.
    iedb_t = iedb[iedb.epitope_class.isin(["T-cell", "MHC-ligand"])
                  & iedb.peptide_length.between(MIN_EPITOPE_LEN, MAX_EPITOPE_LEN)]
    log(f"Using {len(iedb_t)} T-cell / MHC-ligand mappings "
        f"({iedb_t.peptide.nunique()} unique peptides, "
        f"{MIN_EPITOPE_LEN}-{MAX_EPITOPE_LEN} aa) as positives")
    dropped = iedb[iedb.epitope_class.isin(["T-cell", "MHC-ligand"])
                   & ~iedb.peptide_length.between(MIN_EPITOPE_LEN, MAX_EPITOPE_LEN)]
    log(f"  excluded {dropped.peptide.nunique()} longer T-cell/MHC-ligand peptides "
        f"(peptide-scan constructs, not class I ligands)")

    pred = label_peptides(pred, iedb_t)
    n_pos = int(pred.is_iedb_positive.sum())
    log(f"Labelled {n_pos} of {len(pred)} peptides as IEDB-supported "
        f"({100 * n_pos / len(pred):.2f}%)")

    pred["combined_rank_score"] = (
        pred.presentation_score.rank(pct=True).fillna(0) * 0.5
        + pred.min_within_type_identity.rank(pct=True).fillna(0) * 0.5)

    pred.to_csv(os.path.join(DIR_OUT, "labelled_peptides.tsv"), sep="\t", index=False)

    scores = [("presentation_score", True), ("affinity", False),
              ("min_within_type_identity", True), ("combined_rank_score", True)]

    rows = []
    for (vt, gene), grp in pred.groupby(["virus_type", "gene"]):
        if grp.is_iedb_positive.sum() < MIN_POSITIVES:
            continue
        for col, hib in scores:
            if col not in grp.columns:
                continue
            res = evaluate(grp, col, hib)
            if res:
                rows.append(dict(virus_type=vt, gene=gene, score=col, **res))

    per_gene = pd.DataFrame(rows)
    if per_gene.empty:
        log("No gene had enough curated positives to benchmark.")
        return 1
    per_gene.to_csv(os.path.join(DIR_OUT, "benchmark_per_gene.tsv"), sep="\t", index=False)
    log(f"Wrote benchmark_per_gene.tsv ({per_gene.gene.nunique()} genes scored)")

    # Gene-grouped summary: mean across genes, so one well-studied antigen cannot dominate.
    agg = (per_gene.groupby("score")
                   .agg(genes_scored=("gene", "nunique"),
                        mean_auc=("auc", "mean"), median_auc=("auc", "median"),
                        min_auc=("auc", "min"), max_auc=("auc", "max"),
                        mean_average_precision=("average_precision", "mean"),
                        mean_precision_at_50=("precision_at_50", "mean"),
                        mean_top_decile_enrichment=("top_decile_enrichment", "mean"),
                        total_positives=("n_positives", "sum"))
                   .round(4).reset_index().sort_values("mean_auc", ascending=False))
    agg.to_csv(os.path.join(DIR_OUT, "benchmark_summary.tsv"), sep="\t", index=False)

    log("Gene-grouped benchmark (positives vs UNTESTED background -- not true negatives):")
    for r in agg.to_dict("records"):
        log(f"    {r['score']:<28} genes={r['genes_scored']:>2} "
            f"mean AUC={r['mean_auc']:.3f} (range {r['min_auc']:.3f}-{r['max_auc']:.3f}) "
            f"mean AP={r['mean_average_precision']:.3f} "
            f"P@50={r['mean_precision_at_50']:.3f} "
            f"top-decile enrichment={r['mean_top_decile_enrichment']:.2f}x")

    log("Per-gene AUC for the MHCflurry presentation score:")
    ps = per_gene[per_gene.score == "presentation_score"].sort_values("auc", ascending=False)
    for r in ps.to_dict("records"):
        log(f"    {r['virus_type']} {r['gene']:<6} AUC={r['auc']:.3f} "
            f"n_pos={r['n_positives']:>4} of {r['n_peptides']}")

    log("Caveat: IEDB coverage is uneven and absence of a record is not a negative, so "
        "these AUCs measure recovery of known epitopes, not true discrimination.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
