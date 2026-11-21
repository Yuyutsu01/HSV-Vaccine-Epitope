"""
Step 21 (PENDING WORK B6): Integration and ranking of candidate epitopes.

Proposal pipeline step 8. For every candidate epitope, record predicted binding,
conservation, and whether it overlaps a positively selected site; keep epitopes that are
predicted, highly conserved and NOT positively selected; highlight those shared by HSV-1 and
HSV-2; note N-glycosylation sequons that could shield an epitope.

Inputs (each optional -- the table records which evidence was available):
  - 13_epitopes/class1_peptide_predictions.tsv     class I binding + conservation (step15)
  - 11_conservation/per_residue_conservation.tsv   per-residue conservation (step13)
  - 12_iedb/iedb_residue_flags.tsv                 known-epitope coverage (step14)
  - 14_selection/selected_sites_final.tsv          FDR-corrected, concordant positively
                                                   selected sites (step18). Absent until the
                                                   HyPhy stages finish; the selection filter
                                                   is then reported as NOT YET APPLIED rather
                                                   than silently passing.
  - optional class II / B-cell predictions, if the licensed tools have been run:
      13_epitopes/classii_predictions.tsv  (columns: virus_type, gene, peptide, ref_start,
                                            ref_end, allele, percentile_rank)
      13_epitopes/bcell_linear_scores.tsv  (columns: virus_type, gene, ref_residue, score)

Outputs:
  - 16_shortlist/epitope_candidates_ranked.tsv   all binders, scored and annotated
  - 16_shortlist/shortlist_primary.tsv           conserved, predicted, not positively selected
  - 16_shortlist/shortlist_cross_protective.tsv  the above and identical in both types
  - 16_shortlist/shortlist_by_gene.tsv
  - 16_shortlist/integration.log
"""

import os
import re
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from Bio import SeqIO

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")
DIR_EPI = os.path.join(DIR_BASE, "13_epitopes")
DIR_CONS = os.path.join(DIR_BASE, "11_conservation")
DIR_IEDB = os.path.join(DIR_BASE, "12_iedb")
DIR_SEL = os.path.join(DIR_BASE, "14_selection")
DIR_OUT = os.path.join(DIR_BASE, "16_shortlist")

REF_ACC = {"HSV1": "NC_001806.2", "HSV2": "NC_001798.2"}
PRIORITY_20 = ["US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53",
               "UL10", "UL49A", "RL2", "RS1", "UL54", "UL48", "UL46", "UL47",
               "UL49", "UL23", "UL30"]

CONSERVED_MIN_IDENTITY = 99.0     # every residue of the peptide, isolate-weighted
LOG_PATH = os.path.join(DIR_OUT, "integration.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def load_reference_proteins():
    out = {}
    for vt, acc in REF_ACC.items():
        rec = SeqIO.read(os.path.join(DIR_REF, f"{acc}.gbk"), "genbank")
        seen = set()
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            gene = (feat.qualifiers.get("gene", [""])[0]
                    or feat.qualifiers.get("locus_tag", [""])[0])
            if gene in PRIORITY_20 and gene not in seen:
                prot = feat.qualifiers.get("translation", [""])[0]
                if prot:
                    out[(vt, gene)] = prot
                    seen.add(gene)
    return out


def glycosylation_sequons(refs):
    """N-X-S/T with X != P. Returns {(type, gene): set(residues covered by a sequon)}."""
    out = {}
    for key, prot in refs.items():
        covered = set()
        for m in re.finditer(r"N[^P][ST]", prot):
            covered.update(range(m.start() + 1, m.start() + 4))
        out[key] = covered
    return out


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 21: integration and ranking of candidate epitopes")

    pred_path = os.path.join(DIR_EPI, "class1_peptide_predictions.tsv")
    if not os.path.exists(pred_path):
        log(f"ABORT: {pred_path} missing -- run step15 first")
        return 1
    pred = pd.read_csv(pred_path, sep="\t")
    log(f"Loaded {len(pred)} predicted peptides")

    cand = pred[pred.is_binder].copy()
    log(f"{len(cand)} peptides pass the class I binder cut")

    refs = load_reference_proteins()
    sequons = glycosylation_sequons(refs)

    # Known IEDB epitope coverage, per residue.
    iedb_path = os.path.join(DIR_IEDB, "iedb_residue_flags.tsv")
    known = {}
    if os.path.exists(iedb_path):
        fl = pd.read_csv(iedb_path, sep="\t")
        known = {(r.virus_type, r.gene, int(r.ref_residue)): bool(r.known_epitope)
                 for r in fl.itertuples()}
        log(f"Loaded IEDB residue flags ({sum(known.values())} residues inside a known epitope)")
    else:
        log("WARNING: IEDB residue flags missing -- known-epitope column will be blank")

    # FDR-corrected, concordant positively selected sites.
    sel_path = os.path.join(DIR_SEL, "selected_sites_final.tsv")
    selected, selection_available = set(), False
    if os.path.exists(sel_path):
        sel = pd.read_csv(sel_path, sep="\t")
        selected = {(r.cohort, r.gene, int(r.ref_residue)) for r in sel.itertuples()}
        selection_available = True
        log(f"Loaded {len(selected)} positively selected site(s) "
            f"(q<0.05, >=2 methods) from step18")
    else:
        log("NOTE: selected_sites_final.tsv not present yet (HyPhy stages still running). "
            "The positive-selection filter is recorded as NOT YET APPLIED, not as passed.")

    # Optional licensed-tool outputs.
    # step29 writes the full 610k-row table gzipped; accept either form.
    classii_path = next(
        (c for c in (os.path.join(DIR_EPI, "classii_predictions.tsv"),
                     os.path.join(DIR_EPI, "classii_predictions.tsv.gz"))
         if os.path.exists(c)),
        os.path.join(DIR_EPI, "classii_predictions.tsv"))
    bcell_path = os.path.join(DIR_EPI, "bcell_linear_scores.tsv")
    classii = pd.read_csv(classii_path, sep="\t") if os.path.exists(classii_path) else None
    bcell = pd.read_csv(bcell_path, sep="\t") if os.path.exists(bcell_path) else None
    log(f"Class II predictions: {'loaded' if classii is not None else 'NOT AVAILABLE (licence-gated)'}")
    # BepiPred-3.0 is NOT licence-gated (free for academic use); it is produced by step28.
    # Only NetMHCIIpan-4.3, the class II predictor, sits behind the DTU academic form.
    log(f"B-cell linear scores: {'loaded' if bcell is not None else 'NOT AVAILABLE (run step28_bcell_predictors.py --bepipred)'}")

    bcell_idx = {}
    if bcell is not None:
        bcell_idx = {(r.virus_type, r.gene, int(r.ref_residue)): float(r.score)
                     for r in bcell.itertuples()}

    classii_idx = {}
    if classii is not None:
        for r in classii.itertuples():
            key = (r.virus_type, r.gene)
            classii_idx.setdefault(key, []).append(
                (int(r.ref_start), int(r.ref_end), float(r.percentile_rank)))

    # Annotate each candidate.
    glyc, known_col, pos_sel, overlap_sel = [], [], [], []
    bc, c2 = [], []
    for r in cand.itertuples():
        key = (r.virus_type, r.gene)
        rng = range(int(r.ref_start), int(r.ref_end) + 1)

        glyc.append(bool(sequons.get(key, set()) & set(rng)))
        known_col.append(any(known.get((r.virus_type, r.gene, i), False) for i in rng))

        if selection_available:
            hits = [i for i in rng if (r.virus_type, r.gene, i) in selected]
            pos_sel.append(bool(hits))
            overlap_sel.append(";".join(str(i) for i in hits))
        else:
            pos_sel.append(None)
            overlap_sel.append("NOT_YET_APPLIED")

        if bcell_idx:
            vals = [bcell_idx.get((r.virus_type, r.gene, i)) for i in rng]
            vals = [v for v in vals if v is not None]
            bc.append(round(float(np.mean(vals)), 4) if vals else float("nan"))
        else:
            bc.append(float("nan"))

        if classii_idx:
            best = [p for s, e, p in classii_idx.get(key, [])
                    if not (e < r.ref_start or s > r.ref_end)]
            c2.append(round(min(best), 4) if best else float("nan"))
        else:
            c2.append(float("nan"))

    cand["overlaps_glycosylation_sequon"] = glyc
    cand["inside_known_iedb_epitope"] = known_col
    cand["overlaps_positively_selected_site"] = pos_sel
    cand["positively_selected_residues"] = overlap_sel
    cand["bcell_linear_mean_score"] = bc
    cand["classii_best_percentile_rank"] = c2
    cand["selection_filter_applied"] = selection_available

    cand["highly_conserved"] = cand.min_within_type_identity >= CONSERVED_MIN_IDENTITY

    # Low-complexity flag (added 2026-11-21, verification pass finding V11).
    #
    # The HSV proteome carries serine-, glycine-, alanine- and proline-rich tandem repeats,
    # heavily in RL2/ICP0 and RS1/ICP4. Enumerating every 9-mer therefore produces peptides
    # such as SSSSSSSSA and GAGGGGGGR, which are perfectly conserved, score as binders, and
    # are not credible epitopes: class I predictors are known to behave poorly on
    # homopolymers, and a repeat peptide cannot be assigned to a unique position in the
    # protein (SAASSSSSA occurs at RL2 595, 602 and 609 in HSV-2).
    #
    # Flagged rather than dropped, so the shortlist's size is not changed silently by a
    # filter introduced after the fact. To exclude them instead, add
    # `& ~cand.low_complexity` to the `primary` mask below.
    def _shannon(pep):
        counts = {}
        for ch in pep:
            counts[ch] = counts.get(ch, 0) + 1
        n = len(pep)
        return -sum((c / n) * np.log2(c / n) for c in counts.values()) if n else 0.0

    def _max_frac(pep):
        if not len(pep):
            return 1.0
        counts = {}
        for ch in pep:
            counts[ch] = counts.get(ch, 0) + 1
        return max(counts.values()) / len(pep)

    cand["aa_entropy_bits"] = cand.peptide.map(_shannon).round(3)
    cand["max_residue_fraction"] = cand.peptide.map(_max_frac).round(3)
    cand["n_distinct_residues"] = cand.peptide.map(lambda s: len(set(s)))
    cand["low_complexity"] = ((cand.aa_entropy_bits < 2.2)
                              | (cand.max_residue_fraction >= 0.5)
                              | (cand.n_distinct_residues <= 4))

    # Composite rank: binding strength and conservation, both as percentiles, with a penalty
    # for a glycosylation sequon and a bonus for cross-type identity. Deliberately simple and
    # inspectable -- a ranking aid, not a trained model (proposal, ML section).
    cand["binding_percentile"] = (-cand.affinity).rank(pct=True)
    cand["conservation_percentile"] = cand.min_within_type_identity.rank(pct=True)
    cand["priority_score"] = (
        0.45 * cand.binding_percentile.fillna(0)
        + 0.35 * cand.conservation_percentile.fillna(0)
        + 0.15 * cand.fully_conserved_both_types.astype(float)
        - 0.05 * cand.overlaps_glycosylation_sequon.astype(float)).round(4)

    cand = cand.sort_values("priority_score", ascending=False)
    cand.to_csv(os.path.join(DIR_OUT, "epitope_candidates_ranked.tsv"),
                sep="\t", index=False)
    log(f"Wrote epitope_candidates_ranked.tsv ({len(cand)} candidates)")

    not_selected = (cand.overlaps_positively_selected_site != True) \
        if selection_available else pd.Series(True, index=cand.index)
    primary = cand[cand.highly_conserved & not_selected]
    primary.to_csv(os.path.join(DIR_OUT, "shortlist_primary.tsv"), sep="\t", index=False)

    cross = primary[primary.fully_conserved_both_types].copy()

    # Independent check on cross-protective status, not routed through the alignment.
    # "Identical between types at every aligned column" is not the same claim as "this
    # peptide occurs in the other type's protein": an indel inside the window makes the
    # residues correspond one by one while the contiguous peptide exists in only one
    # type. The synthesised peptide is what matters, so verify it directly.
    refs = load_reference_proteins()
    other = {"HSV1": "HSV2", "HSV2": "HSV1"}
    present = [
        bool(refs.get((other[vt], g), "").find(pep) >= 0)
        for vt, g, pep in zip(cross.virus_type, cross.gene, cross.peptide)]
    cross["present_in_other_type_reference"] = present
    # This check is the DEFINITION of cross-protective, not an invariant, so peptides
    # that fail it are excluded -- but to their own file, with the reason, so the
    # exclusion is auditable and the count is reportable. Anything still failing after
    # the exclusion would be a bug, and the assertion below says so.
    missing = cross[~cross.present_in_other_type_reference]
    if len(missing):
        missing.to_csv(os.path.join(DIR_OUT, "excluded_not_in_other_type.tsv"),
                       sep="\t", index=False)
        log(f"  excluded {len(missing)} peptide(s) flagged conserved in both types that "
            f"do not occur in the other type's reference protein (indel inside the "
            f"window; written to excluded_not_in_other_type.tsv): "
            f"{sorted(set(missing.peptide))}")
        cross = cross[cross.present_in_other_type_reference]
    assert bool(cross.present_in_other_type_reference.all()), \
        "cross-protective shortlist still contains a peptide absent from the other type"

    cross.to_csv(os.path.join(DIR_OUT, "shortlist_cross_protective.tsv"),
                 sep="\t", index=False)

    by_gene = (primary.groupby(["virus_type", "gene"])
                      .agg(shortlisted=("peptide", "nunique"),
                           cross_protective=("fully_conserved_both_types", "sum"),
                           inside_known_epitope=("inside_known_iedb_epitope", "sum"),
                           glycosylation_overlap=("overlaps_glycosylation_sequon", "sum"),
                           best_priority=("priority_score", "max"))
                      .reset_index().sort_values("shortlisted", ascending=False))
    by_gene.to_csv(os.path.join(DIR_OUT, "shortlist_by_gene.tsv"), sep="\t", index=False)

    log(f"Primary shortlist: {len(primary)} peptides "
        f"(binder, every residue >= {CONSERVED_MIN_IDENTITY}% conserved within type"
        + (", not positively selected)" if selection_available
           else ", selection filter PENDING)"))
    log(f"Cross-protective subset: {len(cross)} peptides identical across both types")
    log(f"    distinct peptide sequences: {primary.peptide.nunique()} primary, "
        f"{cross.peptide.nunique()} cross-protective "
        f"(a row is one peptide in one virus type, so a cross-type peptide appears twice)")
    log(f"    low-complexity rows (repeat-derived, flagged NOT dropped): "
        f"{int(primary.low_complexity.sum())} of {len(primary)} primary "
        f"({100.0 * primary.low_complexity.mean():.1f}%), "
        f"{int(cross.low_complexity.sum())} of {len(cross)} cross-protective; "
        f"{int(primary.nlargest(50, 'priority_score').low_complexity.sum())} in the top 50")
    log("Shortlist by gene:")
    for r in by_gene.to_dict("records"):
        log(f"    {r['virus_type']} {r['gene']:<6} {r['shortlisted']:>4} shortlisted, "
            f"{int(r['cross_protective']):>3} cross-protective, "
            f"{int(r['inside_known_epitope']):>3} inside a known IEDB epitope, "
            f"{int(r['glycosylation_overlap']):>3} overlap a glycosylation sequon")

    log("Top 15 cross-protective candidates:")
    cols = ["virus_type", "gene", "peptide", "ref_start", "ref_end", "best_allele",
            "affinity", "min_within_type_identity", "inside_known_iedb_epitope",
            "priority_score"]
    for r in cross.head(15)[cols].to_dict("records"):
        log(f"    {r['virus_type']} {r['gene']:<6} {r['peptide']:<11} "
            f"{int(r['ref_start']):>4}-{int(r['ref_end']):<4} {str(r['best_allele']):<12} "
            f"aff={r['affinity']:>8.1f}nM cons={r['min_within_type_identity']:.1f}% "
            f"known={'Y' if r['inside_known_iedb_epitope'] else 'n'} "
            f"score={r['priority_score']}")

    if not selection_available:
        log("REMINDER: re-run this step after step17/step18 so the positive-selection filter "
            "is actually applied. The shortlist is provisional until then.")
    log("Class I population coverage is computed separately by step23 over the per-allele "
        "binder table; see 18_population_coverage/coverage_summary.tsv.")
    if classii is not None:
        n_strong = int((cand.classii_best_percentile_rank <= 1).sum())
        log(f"Class II: classii_best_percentile_rank is the BEST rank over all "
            f"{classii.allele.nunique()} alleles and every 15-mer overlapping the peptide, "
            f"so it is a minimum over roughly a hundred values and is NOT comparable to a "
            f"single-allele rank. It is below 1% for {n_strong} of {len(cand)} candidates "
            f"({100.0 * n_strong / len(cand):.0f}%), which makes it useless as a filter and "
            f"is why it is carried as an annotation only. Per-allele class II binding is in "
            f"13_epitopes/classii_binders_by_allele.tsv.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
