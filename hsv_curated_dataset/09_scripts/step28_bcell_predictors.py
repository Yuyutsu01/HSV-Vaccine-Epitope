"""
Step 28: B-cell epitope layer -- BepiPred-3.0 (linear) and DiscoTope-3.0 (conformational).

Why this exists
---------------
Until this step the pipeline carried a CD8 (HLA class I) arm only. Methods section 16
listed the B-cell arm as pending because both predictors are published by DTU Health Tech
and the earlier audit assumed they were gated behind the same academic-licence form as
NetMHCIIpan-4.3. They are not: both are distributed as source with a licence that permits
academic non-commercial use, so both can be run locally and reproducibly.

  BepiPred-3.0    ESM-2 (650M) protein-language-model encodings -> 5-fold FFNN ensemble,
                  per-residue antibody-contact probability. Sequence input only.
                  Licence: free for academic / non-commercial use (DTU Health Tech).
                  Cite: Clifford et al. 2022, Protein Science, doi:10.1002/pro.4497

  DiscoTope-3.0   ESM-IF1 inverse-folding representations + structural features
                  -> 100-model XGBoost ensemble, per-residue epitope propensity on a
                  3D structure. Structure input, single-chain mode.
                  Licence: CC BY-NC 4.0.
                  Cite: Hoie et al. 2024, Front. Immunol., doi:10.3389/fimmu.2024.1322712

This script does not run the predictors. It is the adapter: the predictors are third-party
code in their own environments (see PROVENANCE below), and this step converts their raw
output into the pipeline's coordinate system, with the checks that conversion needs.

What it checks rather than assumes
----------------------------------
  * Every BepiPred residue letter must equal the reference protein residue at the same
    1-based position. A single mismatch means the FASTA and the output are out of register
    and the whole file is rejected -- a silent off-by-one here would put every B-cell score
    on the wrong residue.
  * Every DiscoTope row is joined to a reference residue through the existing
    17_structures/structure_residue_annotation.tsv map (pdb_id, chain, pdb_resnum), and the
    residue letter is checked against that table's `aa`. Rows outside the solved,
    mappable region are reported as unmapped, not dropped silently.
  * Counts per gene are logged on both sides so the join can be audited.

PROVENANCE (recorded here because it is not reconstructible from the outputs)
----------------------------------------------------------------------------
  BepiPred-3.0  github.com/UberClifford/BepiPred3.0-Predictor, bp3 0.0.12.7,
                fair-esm 1.0.3, esm2_t33_650M_UR50D, torch 2.5.1+cu121 on CUDA.
                Run: bepipred3_CLI.py -pred vt_pred -rolling_window_size 9
                (default variable threshold 0.1512, default top 20%).
  DiscoTope-3.0 github.com/Magnushhoie/discotope3_web, bundled 100-model XGBoost
                ensemble from models.zip, esm_if1_gvp4_t16_142M_UR50,
                torch 2.5.1+cu121, torch_geometric 2.7.0, biotite 1.6.0, xgboost 3.1.2.
                Run: --struc_type solved, single-chain mode (the default and the mode the
                model was trained in), so each chain is scored in isolation. For 2C36 and
                3M1C that means the partner chain's buried interface is scored as exposed;
                those interfaces are flagged in the output rather than corrected.

Inputs (paths given on the command line, since the predictors live outside the project):
  --bepipred   raw_output.csv from BepiPred-3.0
  --discotope  output directory from DiscoTope-3.0 (searched for *_discotope3.csv)

Outputs:
  13_epitopes/bcell_linear_scores.tsv          virus_type, gene, ref_residue, score, ...
  13_epitopes/bcell_conformational_scores.tsv  pdb_id, chain, pdb_resnum, virus_type,
                                               gene, ref_residue, discotope_score, rsa
  13_epitopes/bcell_predictor_provenance.tsv   tool, version, model, licence, citation
  13_epitopes/bcell.log
"""

import os
import re
import sys
import glob
import argparse
from datetime import datetime

import pandas as pd
from Bio import SeqIO

DIR_BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")
DIR_EPI = os.path.join(DIR_BASE, "13_epitopes")
DIR_STRUC = os.path.join(DIR_BASE, "17_structures")
LOG_PATH = os.path.join(DIR_EPI, "bcell.log")

REF_ACC = {"HSV1": "NC_001806.2", "HSV2": "NC_001798.2"}
PRIORITY_20 = [
    "RL2", "RS1", "UL1", "UL10", "UL22", "UL23", "UL27", "UL30", "UL44", "UL46",
    "UL47", "UL48", "UL49", "UL49A", "UL53", "UL54", "US4", "US6", "US7", "US8",
]

# BepiPred-3.0 defaults, recorded so the thresholds in the output are not folklore.
BP3_VAR_THRESHOLD = 0.1512
BP3_ROLLING_WINDOW = 9
# DiscoTope-3.0 does not ship a single recommended cutoff; the paper reports performance
# across thresholds. 0.90 is the value the DTU web server highlights as high-confidence.
DT3_HIGH_CONF = 0.90

PROVENANCE = [
    dict(tool="BepiPred-3.0", version="bp3 0.0.12.7",
         model="esm2_t33_650M_UR50D + 5-fold FFNN ensemble (BP3C50IDFFNN)",
         mode="vt_pred, threshold 0.1512, rolling window 9",
         licence="free for academic / non-commercial use (DTU Health Tech)",
         citation="Clifford et al. 2022 Protein Science doi:10.1002/pro.4497",
         source="github.com/UberClifford/BepiPred3.0-Predictor"),
    dict(tool="DiscoTope-3.0", version="discotope3_web (100-model XGBoost ensemble)",
         model="esm_if1_gvp4_t16_142M_UR50 inverse folding + structural features",
         mode="struc_type=solved, single-chain (training mode)",
         licence="CC BY-NC 4.0",
         citation="Hoie et al. 2024 Front. Immunol. doi:10.3389/fimmu.2024.1322712",
         source="github.com/Magnushhoie/discotope3_web"),
]


def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_EPI, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def load_reference_proteins():
    """Same extraction as step15, so B-cell and class I coordinates are the same."""
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


# ------------------------------------------------------------------ BepiPred-3.0

def adapt_bepipred(path, refs):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    need = ["Accession", "Residue", "BepiPred-3.0 score",
            "BepiPred-3.0 linear epitope score"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        log(f"ABORT: BepiPred output missing columns {missing}")
        return None
    log(f"BepiPred rows: {len(df)}, accessions: {df.Accession.nunique()}")

    rows, mismatches, unknown = [], [], set()
    for acc, grp in df.groupby("Accession", sort=False):
        m = re.match(r"^(HSV1|HSV2)_(.+)$", str(acc))
        if not m:
            unknown.add(acc)
            continue
        vt, gene = m.group(1), m.group(2)
        prot = refs.get((vt, gene))
        if prot is None:
            unknown.add(acc)
            continue
        if len(grp) != len(prot):
            mismatches.append(f"{acc}: {len(grp)} scored residues vs "
                              f"{len(prot)} reference residues")
            continue
        bad = 0
        for i, (res, s_raw, s_lin) in enumerate(zip(
                grp.Residue, grp["BepiPred-3.0 score"],
                grp["BepiPred-3.0 linear epitope score"]), start=1):
            if str(res).strip().upper() != prot[i - 1]:
                bad += 1
                continue
            rows.append(dict(virus_type=vt, gene=gene, ref_residue=i,
                             aa=prot[i - 1],
                             score=round(float(s_lin), 5),
                             bepipred_raw_score=round(float(s_raw), 5),
                             bepipred_linear_score=round(float(s_lin), 5),
                             predicted_epitope_residue=bool(
                                 float(s_raw) >= BP3_VAR_THRESHOLD)))
        if bad:
            mismatches.append(f"{acc}: {bad} residue-letter mismatches")

    if unknown:
        log(f"WARNING: {len(unknown)} BepiPred accessions not in the 20-gene panel: "
            f"{sorted(unknown)[:5]}")
    if mismatches:
        log("ABORT: BepiPred output is out of register with the reference proteins. "
            "Every per-residue score would be assigned to the wrong residue.")
        for mm in mismatches[:10]:
            log(f"    {mm}")
        return None

    out = pd.DataFrame(rows)
    log(f"BepiPred adapted: {len(out)} residues over "
        f"{out.groupby(['virus_type', 'gene']).ngroups} proteins, all in register")
    return out


# ----------------------------------------------------------------- DiscoTope-3.0

def adapt_discotope(out_dir):
    files = sorted(glob.glob(os.path.join(out_dir, "**", "*_discotope3.csv"),
                             recursive=True))
    if not files:
        log(f"WARNING: no *_discotope3.csv under {out_dir}; "
            "conformational layer not produced")
        return None
    log(f"DiscoTope CSV files: {len(files)}")

    ann_path = os.path.join(DIR_STRUC, "structure_residue_annotation.tsv")
    if not os.path.exists(ann_path):
        log(f"ABORT: {ann_path} missing; run step22 first (it owns the "
            "pdb_resnum -> ref_residue map)")
        return None
    ann = pd.read_csv(ann_path, sep="\t")
    amap = {(str(r.pdb_id).upper(), str(r.chain), int(r.pdb_resnum)):
            (r.virus_type, r.gene, int(r.ref_residue), str(r.aa))
            for r in ann.itertuples()}
    log(f"coordinate map: {len(amap)} (pdb, chain, resnum) keys from step22")

    rows, unmapped, aa_bad = [], 0, 0
    for f in files:
        d = pd.read_csv(f)
        # 'DiscoTope-3.0_score' is not a valid Python identifier, so itertuples would
        # rename it positionally. Rename it explicitly instead of trusting the position.
        score_cols = [c for c in d.columns if "DiscoTope" in c and "score" in c.lower()]
        if len(score_cols) != 1:
            log(f"ABORT: expected exactly one DiscoTope score column in {f}, "
                f"found {score_cols}")
            return None
        d = d.rename(columns={score_cols[0]: "dt3_score"})
        for r in d.itertuples():
            tag = str(r.pdb)
            m = re.match(r"^(.+?)_([A-Za-z0-9])$", tag)
            if not m:
                continue
            pdb_id, chain = m.group(1).upper(), m.group(2)
            key = (pdb_id, chain, int(r.res_id))
            hit = amap.get(key)
            if hit is None:
                unmapped += 1
                continue
            vt, gene, ref_res, aa = hit
            if str(r.residue).strip().upper() != aa.strip().upper():
                aa_bad += 1
                continue
            score = float(r.dt3_score)
            rows.append(dict(
                pdb_id=pdb_id, chain=chain, pdb_resnum=int(r.res_id),
                virus_type=vt, gene=gene, ref_residue=ref_res, aa=aa,
                discotope_score=round(score, 5),
                discotope_high_confidence=bool(score >= DT3_HIGH_CONF),
                rsa_discotope=round(float(r.rsa), 5),
                chain_length_scored=int(r.length)))

    log(f"DiscoTope rows mapped: {len(rows)}; unmapped (outside the step22 map, i.e. "
        f"chains or residues with no reference assignment): {unmapped}; "
        f"residue-letter mismatches dropped: {aa_bad}")
    if aa_bad:
        log("WARNING: residue-letter mismatches mean the PDB numbering in the DiscoTope "
            "output and in structure_residue_annotation.tsv disagree for those rows. "
            "They are dropped, not guessed.")
    if not rows:
        return None
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bepipred", help="raw_output.csv from BepiPred-3.0")
    ap.add_argument("--discotope", help="output directory from DiscoTope-3.0")
    args = ap.parse_args()

    log("STEP 28: B-cell epitope layer (BepiPred-3.0 linear, DiscoTope-3.0 structural)")
    os.makedirs(DIR_EPI, exist_ok=True)

    refs = load_reference_proteins()
    log(f"reference proteins: {len(refs)} "
        f"({sum(len(p) for p in refs.values())} residues)")

    wrote = []

    if args.bepipred:
        if not os.path.exists(args.bepipred):
            log(f"ABORT: {args.bepipred} not found")
            return 2
        lin = adapt_bepipred(args.bepipred, refs)
        if lin is None:
            return 2
        p = os.path.join(DIR_EPI, "bcell_linear_scores.tsv")
        lin.to_csv(p, sep="\t", index=False)
        wrote.append(p)
        log(f"wrote bcell_linear_scores.tsv ({len(lin)} rows)")
        n_epi = int(lin.predicted_epitope_residue.sum())
        log(f"    residues at or above the BepiPred threshold {BP3_VAR_THRESHOLD}: "
            f"{n_epi} of {len(lin)} ({100.0 * n_epi / len(lin):.1f}%)")
        log("    per-gene mean linear epitope score, highest 8:")
        top = (lin.groupby(["virus_type", "gene"]).score.mean()
               .sort_values(ascending=False).head(8))
        for (vt, gene), v in top.items():
            log(f"        {vt} {gene:6} {v:.4f}")

    if args.discotope:
        if not os.path.isdir(args.discotope):
            log(f"ABORT: {args.discotope} is not a directory")
            return 2
        conf = adapt_discotope(args.discotope)
        if conf is not None:
            p = os.path.join(DIR_EPI, "bcell_conformational_scores.tsv")
            conf.to_csv(p, sep="\t", index=False)
            wrote.append(p)
            log(f"wrote bcell_conformational_scores.tsv ({len(conf)} rows)")
            log("    per chain: residues scored, mean score, high-confidence count")
            for (pdb_id, chain), g in conf.groupby(["pdb_id", "chain"]):
                log(f"        {pdb_id} {chain}  {g.virus_type.iloc[0]} "
                    f"{g.gene.iloc[0]:6} n={len(g):4}  "
                    f"mean={g.discotope_score.mean():.4f}  "
                    f"high-conf(>={DT3_HIGH_CONF})={int(g.discotope_high_confidence.sum())}")

    prov = pd.DataFrame(PROVENANCE)
    p = os.path.join(DIR_EPI, "bcell_predictor_provenance.tsv")
    prov.to_csv(p, sep="\t", index=False)
    wrote.append(p)

    log("")
    log("wrote:")
    for p in wrote:
        log(f"    {os.path.relpath(p, DIR_BASE)}")
    log("")
    log("Scope note: both predictors are run on the TYPE REFERENCE protein and on solved "
        "structures of reference-like strains, not per isolate. A B-cell score is "
        "therefore a property of the reference antigen, and the conservation layer is "
        "what carries the per-isolate variability. Nothing here is validated against "
        "experimentally mapped HSV B-cell epitopes in this study.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
