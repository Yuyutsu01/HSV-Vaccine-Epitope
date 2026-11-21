"""
Step 29: CD4 / HLA class II binding predictions (NetMHCIIpan-4.3) -> pipeline coordinates.

Why this is an adapter and not a predictor
------------------------------------------
NetMHCIIpan-4.3 is the one tool in this design that is genuinely licence-gated: it is
distributed only through DTU Health Tech's academic registration form, and it has no public
source repository. The predictions are therefore produced either from a local install or
from the DTU web server, and this step converts whatever came back into the pipeline's
coordinate system. It never calls the predictor.

Panel
-----
The 27-allele IEDB class II reference set (Greenbaum et al.), the direct analogue of the
27-allele class I panel used in step15. Each of the 27 names was checked against
NetMHCIIpan-4.3's own MHC_allele_names.txt before the runs were submitted; all 27 are
present. Note that the server's free-text allele spelling (DRB1_0101,
HLA-DQA10501-DQB10201) differs from its own dropdown spelling (DRB1*0101,
HLA-DQA1*0501-DQB1*0201) and from IEDB's (HLA-DRB1*01:01, HLA-DQA1*05:01/DQB1*02:01).
This step normalises all three onto IEDB nomenclature, because that is what the population
coverage step (step23) joins on.

One panel difference is recorded rather than reconciled: IEDB states it has no predictor
for HLA-DPA1*02:01/DPB1*14:01, while NetMHCIIpan-4.3 does. Class II coverage computed here
therefore rests on 27 alleles where an IEDB-computed figure would rest on 26.

Parsing
-------
The NetMHCIIpan `-xls` output is a tab-separated file with a TWO-ROW header: the first row
carries the allele names, each above its own block of columns, mostly padded with empty
cells; the second row carries the column names within each block (Rank_EL, Score_EL,
Rank_BA, Affinity, ...). The leading columns (Pos, Peptide, ID, Target) belong to no
allele. This layout is not documented in a stable form, so the parser DETECTS it rather
than assuming fixed positions:

  * the header row pair is found by locating the row whose cells contain 'Peptide' and
    'Pos';
  * allele blocks are delimited by the non-empty cells of the row above it;
  * within each block, columns are matched by name, case-insensitively.

Run with --inspect to print what the parser found without writing anything. Do that first
on any new output file.

Checks
------
Every peptide must occur in the reference protein at the position the predictor reports.
A mismatch means the coordinates are wrong, and the step aborts rather than writing a file
that silently puts class II binders on the wrong residues. This is the same check that
step28 applies to the B-cell scores.

Inputs:
  --xls   one or more NetMHCIIpan output files (repeatable, or a glob)

Outputs:
  13_epitopes/classii_predictions.tsv.gz    virus_type, gene, peptide, ref_start, ref_end,
                                            allele, percentile_rank, (+ BA columns if run)
  13_epitopes/classii_binders_by_allele.tsv binders only, for population coverage
  13_epitopes/classii.log
"""

import os
import re
import sys
import glob
import gzip
import argparse
from datetime import datetime

import pandas as pd
from Bio import SeqIO

DIR_BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")
DIR_EPI = os.path.join(DIR_BASE, "13_epitopes")
LOG_PATH = os.path.join(DIR_EPI, "classii.log")

REF_ACC = {"HSV1": "NC_001806.2", "HSV2": "NC_001798.2"}
PRIORITY_20 = [
    "RL2", "RS1", "UL1", "UL10", "UL22", "UL23", "UL27", "UL30", "UL44", "UL46",
    "UL47", "UL48", "UL49", "UL49A", "UL53", "UL54", "US4", "US6", "US7", "US8",
]

# NetMHCIIpan / IEDB class II convention. Looser than class I by design.
STRONG_RANK = 2.0
WEAK_RANK = 10.0


def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_EPI, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def load_reference_proteins():
    """Identical extraction to step15 and step28, so all three arms share coordinates."""
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


# ------------------------------------------------------------------ allele normalisation

def normalise_allele(raw):
    """Map any of the three spellings onto IEDB nomenclature.

    DRB1_0101              -> HLA-DRB1*01:01
    DRB1*0101              -> HLA-DRB1*01:01
    HLA-DQA10501-DQB10201  -> HLA-DQA1*05:01/DQB1*02:01
    HLA-DQA1*0501-DQB1*0201-> HLA-DQA1*05:01/DQB1*02:01
    Anything unrecognised is returned unchanged and reported, never silently dropped.
    """
    a = str(raw).strip()
    if not a:
        return a

    def split_digits(d):
        d = d.replace(":", "")
        return f"{d[:2]}:{d[2:]}" if len(d) >= 4 else d

    # DP / DQ alpha-beta pair
    m = re.match(r"^HLA-(DPA1|DQA1)\*?(\d+:?\d*)-(DPB1|DQB1)\*?(\d+:?\d*)$", a)
    if m:
        return (f"HLA-{m.group(1)}*{split_digits(m.group(2))}"
                f"/{m.group(3)}*{split_digits(m.group(4))}")
    # already IEDB form
    if re.match(r"^HLA-D[PQ]A1\*\d+:\d+/D[PQ]B1\*\d+:\d+$", a):
        return a
    # DR single chain
    m = re.match(r"^(?:HLA-)?(DRB[1345])[_*]?(\d+:?\d*)$", a)
    if m:
        return f"HLA-{m.group(1)}*{split_digits(m.group(2))}"
    if re.match(r"^HLA-DRB[1345]\*\d+:\d+$", a):
        return a
    return a


# ---------------------------------------------------------------------------- parsing

def find_header(rows):
    """Return (allele_row_index, column_row_index) or (None, None)."""
    for i, r in enumerate(rows):
        low = [c.strip().lower() for c in r]
        if "peptide" in low and ("pos" in low or "position" in low):
            return (i - 1 if i > 0 else None), i
    return None, None


def parse_netmhciipan(path):
    """Parse one NetMHCIIpan -xls file into long format. Returns (DataFrame, report)."""
    # The archived raw output is gzipped: each NetMHCIIpan run is ~20 MB of xls and cannot
    # be regenerated (the server deletes results after 24 hours), so it is kept compressed
    # in the repository rather than left out of it. Plain .xls still works.
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if not rows:
        return None, {"error": "file is empty"}

    ai, ci = find_header(rows)
    if ci is None:
        return None, {"error": "no header row containing both 'Peptide' and 'Pos' was found"}

    cols = [c.strip() for c in rows[ci]]
    allele_row = [c.strip() for c in rows[ai]] if ai is not None else [""] * len(cols)
    allele_row += [""] * (len(cols) - len(allele_row))

    # Leading columns are exactly those BEFORE the first non-empty cell of the allele row.
    # Deciding this by column name instead would be wrong: the real output repeats 'Core'
    # inside every allele block (Core, Inverted, Score_EL, Rank_EL, Score_BA, nM, Rank_BA),
    # so a name-based rule strips each block's Core out of its own block.
    first_allele_col = next((i for i, c in enumerate(allele_row) if c), len(cols))
    lead_idx = set(range(first_allele_col))

    # forward-fill the allele row across its block
    cur, allele_of = "", {}
    for i, c in enumerate(cols):
        if allele_row[i]:
            cur = allele_row[i]
        if i not in lead_idx and cur:
            allele_of[i] = cur

    blocks = {}
    for i, al in allele_of.items():
        blocks.setdefault(al, {})[cols[i].lower()] = i

    def lead(name):
        for i in sorted(lead_idx):
            if cols[i].lower() == name:
                return i
        return None

    i_pos, i_pep, i_id = lead("pos"), lead("peptide"), lead("id")
    # The first line is a '#' comment holding the exact command the server ran. Read the
    # %Rank thresholds from it rather than hard-coding a convention: NetMHCIIpan-4.3
    # defaults to -rankS 1 -rankW 5, while 4.0 and the IEDB convention use 2 and 10. Taking
    # them from the file means the reported thresholds are always the ones actually used.
    cmdline = rows[0][0] if rows and rows[0] and rows[0][0].startswith("#") else ""
    def opt(flag, default=None):
        m = re.search(rf"{flag}\s+([0-9.]+)", cmdline)
        return float(m.group(1)) if m else default

    report = {
        "file": os.path.basename(path),
        "command_line": cmdline.lstrip("#").strip(),
        "rank_strong": opt("-rankS"),
        "rank_weak": opt("-rankW"),
        "data_rows": len(rows) - ci - 1,
        "n_alleles": len(blocks),
        "alleles_raw": sorted(blocks),
        "block_columns": sorted({c for b in blocks.values() for c in b}),
        "lead_columns": [cols[i] for i in sorted(lead_idx)],
        "has_pos": i_pos is not None, "has_peptide": i_pep is not None,
        "has_id": i_id is not None,
    }
    if i_pos is None or i_pep is None or i_id is None:
        report["error"] = "missing one of the Pos / Peptide / ID columns"
        return None, report

    def pick(block, *names):
        for n in names:
            if n in block:
                return block[n]
        return None

    out = []
    for r in rows[ci + 1:]:
        if len(r) <= i_pep:
            continue
        pep = r[i_pep].strip()
        if not pep or not pep.isalpha():
            continue
        try:
            pos = int(float(r[i_pos]))
        except (ValueError, IndexError):
            continue
        sid = r[i_id].strip()
        for al, block in blocks.items():
            def val(*names):
                j = pick(block, *names)
                if j is None or j >= len(r) or r[j].strip() == "":
                    return None
                try:
                    return float(r[j])
                except ValueError:
                    return None
            rank_el = val("rank_el", "%rank_el", "rank", "%rank")
            if rank_el is None:
                continue
            out.append(dict(
                seq_id=sid, pos=pos, peptide=pep, allele_raw=al,
                percentile_rank=rank_el,
                score_el=val("score_el", "score"),
                rank_ba=val("rank_ba", "%rank_ba"),
                affinity_nm=val("nm", "affinity", "affinity(nm)", "aff(nm)", "ic50"),
                core=(r[block["core"]].strip()
                      if "core" in block and block["core"] < len(r) else None),
                # NetMHCIIpan-4.3 models INVERTED binding, where the peptide sits in the
                # groove back-to-front. When this flag is set the reported Core is the
                # REVERSE of a substring of the peptide, not a substring of it. In this
                # dataset it fires on 4.9% of DP rows and never on DR or DQ. Carrying the
                # flag is what makes the Core column interpretable downstream.
                inverted=val("inverted")))
    return pd.DataFrame(out), report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xls", action="append", required=True,
                    help="NetMHCIIpan output file or glob (repeatable)")
    ap.add_argument("--inspect", action="store_true",
                    help="report the detected layout and exit without writing")
    args = ap.parse_args()

    paths = []
    for p in args.xls:
        hits = sorted(glob.glob(p))
        paths.extend(hits if hits else [p])
    paths = [p for p in paths if os.path.exists(p)]
    if not paths:
        log(f"ABORT: none of {args.xls} exist")
        return 2

    log("STEP 29: NetMHCIIpan-4.3 class II predictions -> pipeline coordinates")
    log(f"input files: {len(paths)}")

    frames, thresholds = [], set()
    for p in paths:
        df, rep = parse_netmhciipan(p)
        log("")
        log(f"--- {rep.get('file', p)}")
        if rep.get("command_line"):
            log(f"    command: {rep['command_line']}")
        if rep.get("rank_strong") is not None:
            thresholds.add((rep["rank_strong"], rep["rank_weak"]))
            log(f"    thresholds from the run: strong <= {rep['rank_strong']}%, "
                f"weak <= {rep['rank_weak']}%")
        for k in ("data_rows", "n_alleles", "lead_columns", "block_columns"):
            if k in rep:
                log(f"    {k}: {rep[k]}")
        if "alleles_raw" in rep:
            log(f"    alleles: {rep['alleles_raw']}")
            log("    normalised:")
            for a in rep["alleles_raw"]:
                n = normalise_allele(a)
                flag = "" if n != a else "   <-- UNRECOGNISED SPELLING, left as-is"
                log(f"        {a:28} -> {n}{flag}")
        if "error" in rep:
            log(f"    ERROR: {rep['error']}")
            if not args.inspect:
                return 2
            continue
        log(f"    parsed rows: {len(df)}")
        frames.append(df)

    if args.inspect:
        log("")
        log("--inspect: nothing written. Re-run without --inspect to produce the tables.")
        return 0

    if not frames:
        log("ABORT: nothing parsed")
        return 2
    df = pd.concat(frames, ignore_index=True)
    log("")
    log(f"combined rows: {len(df)}")

    # ---- map sequence ids onto (virus_type, gene) and verify coordinates
    refs = load_reference_proteins()
    vt, gn, bad_id, bad_coord = [], [], set(), []
    for r in df.itertuples():
        m = re.match(r"^(HSV1|HSV2)[_|](.+?)$", str(r.seq_id))
        if not m or (m.group(1), m.group(2)) not in refs:
            bad_id.add(r.seq_id)
            vt.append(None); gn.append(None)
            continue
        v, g = m.group(1), m.group(2)
        prot = refs[(v, g)]
        if prot[r.pos - 1:r.pos - 1 + len(r.peptide)] != r.peptide:
            bad_coord.append((v, g, r.pos, r.peptide))
        vt.append(v); gn.append(g)
    df["virus_type"], df["gene"] = vt, gn

    if bad_id:
        log(f"ABORT: {len(bad_id)} sequence id(s) do not match the 'HSV1_GENE' convention "
            f"or are not in the 20-gene panel: {sorted(bad_id)[:5]}")
        return 2
    if bad_coord:
        log(f"ABORT: {len(bad_coord)} peptide(s) do not occur in the reference protein at "
            "the reported position. The coordinates are wrong and nothing was written.")
        for v, g, p, pep in bad_coord[:10]:
            log(f"    {v} {g} pos {p}: {pep} != "
                f"{refs[(v, g)][p - 1:p - 1 + len(pep)]}")
        return 2
    log(f"coordinate check: all {len(df)} peptides occur at the reported position")

    df["allele"] = df.allele_raw.map(normalise_allele)
    unrec = sorted(set(df.loc[df.allele == df.allele_raw, "allele_raw"]))
    if unrec:
        log(f"WARNING: {len(unrec)} allele spelling(s) not recognised and left unchanged; "
            f"step23 will not join these: {unrec}")

    df["ref_start"] = df.pos
    df["ref_end"] = df.pos + df.peptide.str.len() - 1

    # Use the thresholds the runs actually used, not a hard-coded convention. If the input
    # files were produced with different settings, say so and fall back to the defaults
    # rather than applying one file's thresholds to another's scores.
    strong, weak = STRONG_RANK, WEAK_RANK
    if len(thresholds) == 1:
        strong, weak = thresholds.pop()
        log(f"binder thresholds taken from the run: strong <= {strong}%, weak <= {weak}%")
    elif len(thresholds) > 1:
        log(f"WARNING: input files were run with DIFFERENT %Rank thresholds {thresholds}; "
            f"falling back to the script defaults {STRONG_RANK}/{WEAK_RANK} and recording "
            "the raw rank, which is unaffected.")
    df["binder_level"] = pd.cut(
        df.percentile_rank, [-1, strong, weak, 1e9],
        labels=["strong", "weak", "non_binder"]).astype(str)
    df["rank_threshold_strong"], df["rank_threshold_weak"] = strong, weak

    cols = ["virus_type", "gene", "peptide", "ref_start", "ref_end", "allele",
            "percentile_rank", "score_el", "rank_ba", "affinity_nm", "core",
            "binder_level", "allele_raw"]
    out = df[[c for c in cols if c in df.columns]].sort_values(
        ["virus_type", "gene", "ref_start", "allele"])
    # 610,848 rows is ~76 MB uncompressed, which is too large to carry in git sensibly.
    # Gzipped it is ~17 MB, and pandas reads it transparently from the .gz path.
    p = os.path.join(DIR_EPI, "classii_predictions.tsv.gz")
    out.to_csv(p, sep="\t", index=False, compression="gzip")
    log(f"wrote classii_predictions.tsv.gz ({len(out)} rows, "
        f"{out.allele.nunique()} alleles, "
        f"{out.groupby(['virus_type', 'gene']).ngroups} proteins)")

    binders = out[out.percentile_rank <= weak]
    p = os.path.join(DIR_EPI, "classii_binders_by_allele.tsv")
    binders.to_csv(p, sep="\t", index=False)
    log(f"wrote classii_binders_by_allele.tsv ({len(binders)} rows at "
        f"%Rank_EL <= {weak})")

    log("")
    log(f"binder levels (strong <= {strong}%, weak <= {weak}%):")
    for k, v in out.binder_level.value_counts().items():
        log(f"    {k:12} {v:>8}  ({100.0 * v / len(out):.1f}%)")
    log("")
    log("strong binders per allele:")
    s = out[out.binder_level == "strong"].allele.value_counts()
    for a in sorted(out.allele.unique()):
        log(f"    {a:34} {int(s.get(a, 0)):>6}")
    log("")
    log("Scope note: these are predicted class II binders on the reference proteins. No "
        "class II predictor is benchmarked in this study, and %Rank thresholds for class II "
        "are looser than class I by convention, so the continuous rank is what should be "
        "used downstream.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
