"""
Step 15 (PENDING WORK B4, CD8 arm): HLA class I epitope prediction with conservation joined.

Proposal pipeline step 7 (CD8 T-cell arm). NetMHCpan-4.1 requires an academic licence and is
not installable in this environment, so MHCflurry 2.3 (Class1PresentationPredictor, an
already-trained neural predictor) is used as the free equivalent the proposal allows.
Nothing is trained here; the predictor is used as an instrument.

Peptides are enumerated from the type reference proteins so that every prediction carries a
reference residue range, and each peptide is immediately annotated with the per-residue
conservation computed in step13 (isolate-weighted, within type and between types). This is
the join the integration step (B6) consumes.

Allele panel: the 27-allele HLA class I reference panel (Weiskopf/Sette), chosen for broad
global population coverage.

Outputs:
  - 13_epitopes/class1_peptide_predictions.tsv   (per peptide: best allele, affinity, conservation)
  - 13_epitopes/class1_binders_by_allele.tsv     (every peptide x allele passing the binder cut)
  - 13_epitopes/class1_prediction_summary.tsv
  - 13_epitopes/class1.log
"""

import os
import sys
from datetime import datetime

import pandas as pd
from Bio import SeqIO

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")
DIR_CONS = os.path.join(DIR_BASE, "11_conservation")
DIR_OUT = os.path.join(DIR_BASE, "13_epitopes")

REF_ACC = {"HSV1": "NC_001806.2", "HSV2": "NC_001798.2"}
PRIORITY_20 = ["US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53",
               "UL10", "UL49A", "RL2", "RS1", "UL54", "UL48", "UL46", "UL47",
               "UL49", "UL23", "UL30"]

# 27-allele HLA class I reference panel (broad global coverage).
HLA_PANEL = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01",
    "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-A*30:01",
    "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01",
    "HLA-A*68:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01", "HLA-B*40:01",
    "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01", "HLA-B*53:01", "HLA-B*57:01",
    "HLA-B*58:01",
]

PEPTIDE_LENGTHS = (9, 10)
AFFINITY_CUT = 500.0          # nM, conventional binder threshold
PRESENTATION_CUT = 0.50       # MHCflurry presentation score

LOG_PATH = os.path.join(DIR_OUT, "class1.log")


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


def conservation_index():
    """{(gene, type, ref_residue): (within_identity, between_identical, conserved_both)}"""
    path = os.path.join(DIR_CONS, "per_residue_conservation.tsv")
    if not os.path.exists(path):
        log(f"WARNING: {path} missing -- run step13 first; conservation columns will be empty")
        return {}
    df = pd.read_csv(path, sep="\t")
    idx = {}
    for vt, rescol, idcol in (("HSV1", "hsv1_ref_residue", "hsv1_pct_identity_isolate"),
                             ("HSV2", "hsv2_ref_residue", "hsv2_pct_identity_isolate")):
        sub = df[df[rescol].notna()]
        # Between-type identity is taken from the REFERENCE residues, because the peptide
        # being annotated was cut out of a reference protein (see step13).
        for gene, res, ident, between, both in zip(
                sub.gene, sub[rescol], sub[idcol],
                sub.between_type_reference_identical, sub.conserved_both_types):
            idx[(gene, vt, int(res))] = (float(ident), bool(between), bool(both))
    return idx


def enumerate_peptides(refs):
    rows = []
    for (vt, gene), prot in sorted(refs.items()):
        for L in PEPTIDE_LENGTHS:
            for i in range(len(prot) - L + 1):
                rows.append(dict(virus_type=vt, gene=gene, peptide=prot[i:i + L],
                                 length=L, ref_start=i + 1, ref_end=i + L))
    return pd.DataFrame(rows)


def annotate_conservation(df, cons):
    mins, means, allbetween, allboth, covered = [], [], [], [], []
    for vt, gene, s, e in zip(df.virus_type, df.gene, df.ref_start, df.ref_end):
        vals, betw, both, n = [], [], [], 0
        for res in range(s, e + 1):
            hit = cons.get((gene, vt, res))
            if hit:
                n += 1
                vals.append(hit[0])
                betw.append(hit[1])
                both.append(hit[2])
        if vals:
            mins.append(round(min(vals), 2))
            means.append(round(sum(vals) / len(vals), 2))
            # A between-type or both-types claim holds only if EVERY residue of the
            # peptide was checked. Masked columns carry no conservation row, so a
            # partially covered peptide would otherwise inherit all() over a subset --
            # all([]) is True and all() over 8 of 9 residues says nothing about the 9th.
            full = (n == (e - s + 1))
            allbetween.append(bool(full and all(betw)))
            allboth.append(bool(full and all(both)))
        else:
            mins.append(float("nan"))
            means.append(float("nan"))
            allbetween.append(False)
            allboth.append(False)
        covered.append(n)
    df["min_within_type_identity"] = mins
    df["mean_within_type_identity"] = means
    df["all_residues_identical_between_types"] = allbetween
    df["fully_conserved_both_types"] = allboth
    df["residues_with_conservation_data"] = covered
    # Both flags above are now False wherever coverage is partial, so this must hold.
    partial = df[(df.residues_with_conservation_data < df.length)
                 & (df.all_residues_identical_between_types
                    | df.fully_conserved_both_types)]
    if len(partial):
        raise RuntimeError(
            f"{len(partial)} peptide(s) carry a between-type flag on partial residue "
            f"coverage; first: {partial.iloc[0].to_dict()}")
    return df


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 15: HLA class I prediction (MHCflurry 2.3) with conservation annotation")

    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    panel = [a for a in HLA_PANEL if a in set(predictor.supported_alleles)]
    missing = sorted(set(HLA_PANEL) - set(panel))
    log(f"Allele panel: {len(panel)}/{len(HLA_PANEL)} supported"
        + (f"; unsupported: {missing}" if missing else ""))

    refs = load_reference_proteins()
    log(f"Reference proteins: {len(refs)}")

    peptides = enumerate_peptides(refs)
    log(f"Enumerated {len(peptides)} peptides "
        f"({'/'.join(str(x) for x in PEPTIDE_LENGTHS)}-mers) across {len(refs)} proteins")

    cons = conservation_index()
    peptides = annotate_conservation(peptides, cons)

    unique_pep = sorted(peptides.peptide.unique())
    log(f"Predicting {len(unique_pep)} unique peptides x {len(panel)} alleles "
        f"= {len(unique_pep) * len(panel)} predictions")

    # Each allele is passed as its own sample, so the result has one row per
    # peptide x allele. Passing the panel as a single sample would collapse to the best
    # allele only, which leaves no per-allele data for the population coverage step.
    res = predictor.predict(peptides=unique_pep,
                            alleles={a: [a] for a in panel},
                            include_affinity_percentile=True,
                            verbose=0)
    res = res.rename(columns={"sample_name": "allele"})
    log(f"MHCflurry returned {len(res)} peptide x allele rows; "
        f"columns: {list(res.columns)}")

    best = (res.sort_values("presentation_score", ascending=False)
               .drop_duplicates("peptide")
               .set_index("peptide"))
    if "best_allele" not in best.columns:
        best["best_allele"] = best["allele"]

    keep = [c for c in ("best_allele", "affinity", "affinity_percentile",
                        "presentation_score", "processing_score") if c in best.columns]
    merged = peptides.join(best[keep], on="peptide")

    merged["is_binder_affinity"] = merged.affinity <= AFFINITY_CUT
    merged["is_binder_presentation"] = merged.presentation_score >= PRESENTATION_CUT
    merged["is_binder"] = merged.is_binder_affinity | merged.is_binder_presentation

    out1 = os.path.join(DIR_OUT, "class1_peptide_predictions.tsv")
    merged.to_csv(out1, sep="\t", index=False)
    log(f"Wrote {out1} ({len(merged)} rows)")

    # Per-allele detail, restricted to binders so the file stays usable.
    per_allele = res[(res.affinity <= AFFINITY_CUT)
                     | (res.presentation_score >= PRESENTATION_CUT)].copy()
    pep_meta = peptides.drop_duplicates("peptide").set_index("peptide")[
        ["virus_type", "gene", "ref_start", "ref_end", "min_within_type_identity",
         "fully_conserved_both_types"]]
    per_allele = per_allele.join(pep_meta, on="peptide")
    out2 = os.path.join(DIR_OUT, "class1_binders_by_allele.tsv")
    per_allele.to_csv(out2, sep="\t", index=False)
    log(f"Wrote {out2} ({len(per_allele)} peptide-allele binder rows)")

    rows = []
    for (vt, gene), grp in merged.groupby(["virus_type", "gene"]):
        b = grp[grp.is_binder]
        cb = b[b.fully_conserved_both_types]
        rows.append(dict(virus_type=vt, gene=gene, peptides=len(grp), binders=len(b),
                         pct_binders=round(100 * len(b) / len(grp), 2) if len(grp) else 0,
                         binders_fully_conserved_within_type=int(
                             (b.min_within_type_identity >= 99.0).sum()),
                         binders_conserved_both_types=len(cb)))
    summary = pd.DataFrame(rows).sort_values(["virus_type", "gene"])
    out3 = os.path.join(DIR_OUT, "class1_prediction_summary.tsv")
    summary.to_csv(out3, sep="\t", index=False)
    log(f"Wrote {out3}")
    for r in summary.to_dict("records"):
        log(f"    {r['virus_type']} {r['gene']}: {r['binders']}/{r['peptides']} binders "
            f"({r['pct_binders']}%), {r['binders_conserved_both_types']} conserved in both types")
    return 0


if __name__ == "__main__":
    sys.exit(main())
