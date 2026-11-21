"""
Step 10b (CORRECTION): Prepare HyPhy-ready alignments and rebuild per-type phylogenies.

Supersedes step11_prereq2_build_trees.py. Changes:
  - Trees are built from the rebuilt, masked MSAs (step09b/09c), not the old star alignments.
  - PER TYPE ONLY. No combined HSV-1 + HSV-2 trees are built, because a tree spanning two
    species cannot support a single codon model (CORRECTIONS_AND_PENDING_WORK.md, A2).
  - Sequence names are rewritten to short safe IDs up front, with a mapping table, so
    IQ-TREE never mangles them and HyPhy taxon names match exactly. This removes the need
    to parse IQ-TREE rename warnings after the fact.
  - --keep-ident retains identical sequences so every alignment taxon has a tree tip.
  - The reference anchor row (REF_*) is dropped: these cohorts are clinical isolates.
    Reference numbering is recovered from the step09c coordinate maps instead.

Outputs:
  - 06_alignments_v2/hyphy_ready/<cohort>/<gene>.fasta
  - 06_alignments_v2/hyphy_ready/<cohort>/<gene>_idmap.tsv
  - 06_alignments_v2/trees_v2/<cohort>/<gene>.treefile (+ IQ-TREE side files)
  - 06_alignments_v2/trees_v2/tree_summary_v2.tsv
  - 06_alignments_v2/trees_v2/trees.log

Usage:
  python step10b_prep_and_trees.py            # prep + all trees
  python step10b_prep_and_trees.py --prep     # prep only
  python step10b_prep_and_trees.py --genes US6,UL27
"""

import os
import re
import sys
import glob
import time
import argparse
import subprocess
from datetime import datetime

import pandas as pd

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_V2 = os.path.join(DIR_BASE, "06_alignments_v2")
DIR_MASKED = os.path.join(DIR_V2, "masked")
DIR_READY = os.path.join(DIR_V2, "hyphy_ready")
DIR_TREES = os.path.join(DIR_V2, "trees_v2")

IQTREE = os.path.join(DIR_BASE, "bin", "iqtree", "iqtree-2.4.0-Windows", "bin", "iqtree2.exe")
COHORTS = ("HSV1", "HSV2")
THREADS = "4"
BOOTSTRAP = "1000"

LOG_PATH = os.path.join(DIR_TREES, "trees.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_TREES, exist_ok=True)
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


def accession_of(header):
    """First field of the pipe-delimited header is the accession."""
    return header.split("|")[0]


def prep(gene, cohort):
    src = os.path.join(DIR_MASKED, cohort, f"{gene}_masked.fasta")
    if not os.path.exists(src):
        return None
    ids, seqs = read_fasta(src)

    rows, keep_ids, keep_seqs = [], [], []
    n = 0
    for header, seq in zip(ids, seqs):
        if header.startswith("REF_"):
            continue
        n += 1
        safe = f"S{n:04d}"
        keep_ids.append(safe)
        keep_seqs.append(seq)
        m = re.search(r"cluster_size=(\d+)", header)
        rows.append(dict(safe_id=safe, accession=accession_of(header),
                         cluster_size=int(m.group(1)) if m else 1,
                         original_header=header))

    os.makedirs(os.path.join(DIR_READY, cohort), exist_ok=True)
    out_fa = os.path.join(DIR_READY, cohort, f"{gene}.fasta")
    with open(out_fa, "w", encoding="utf-8") as fh:
        for i, s in zip(keep_ids, keep_seqs):
            fh.write(f">{i}\n{s}\n")
    pd.DataFrame(rows).to_csv(
        os.path.join(DIR_READY, cohort, f"{gene}_idmap.tsv"), sep="\t", index=False)

    return dict(gene=gene, cohort=cohort, n_taxa=len(keep_ids),
                length_nt=len(keep_seqs[0]) if keep_seqs else 0, fasta=out_fa)


def build_tree(info):
    gene, cohort = info["gene"], info["cohort"]
    out_dir = os.path.join(DIR_TREES, cohort)
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, gene)

    # A .treefile alone is NOT evidence of a finished run: when IQ-TREE dies during
    # the bootstrap stage it leaves the ML tree behind with no support values, and
    # caching on existence alone silently keeps an unsupported tree. The .iqtree
    # report is written only on completion, and across all 40 cohorts its presence
    # matches bootstrap support exactly, so require it as well.
    if os.path.exists(prefix + ".treefile") and os.path.exists(prefix + ".iqtree"):
        log(f"{gene} {cohort}: tree cached, skipped")
        return dict(**{k: info[k] for k in ("gene", "cohort", "n_taxa", "length_nt")},
                    model="cached", runtime_s=0.0, status="CACHED", safe_kernel=False)
    if os.path.exists(prefix + ".treefile"):
        log(f"{gene} {cohort}: previous tree has no .iqtree report "
            f"(died before bootstrapping finished); re-inferring")

    def invoke(safe):
        cmd = [IQTREE, "-s", info["fasta"], "-m", "MFP", "-B", BOOTSTRAP,
               "--keep-ident", "-T", THREADS, "--prefix", prefix, "--quiet", "-redo"]
        if safe:
            cmd.append("-safe")
        return subprocess.run(cmd, capture_output=True, text=True)

    t0 = time.time()
    res = invoke(safe=False)
    # These cohorts are near-identical sequences, which drives the likelihood
    # derivative to underflow on some genes. `-safe` is IQ-TREE's own documented
    # remedy (a slower, higher-precision likelihood kernel); it is only used on
    # the genes that need it, and the kernel does not change the inferred model.
    used_safe = False
    if res.returncode != 0 and "Numerical underflow" in (res.stderr or "") + (res.stdout or ""):
        log(f"{gene} {cohort}: numerical underflow, retrying with -safe")
        res = invoke(safe=True)
        used_safe = True
    rt = round(time.time() - t0, 1)

    model = ""
    if os.path.exists(prefix + ".iqtree"):
        with open(prefix + ".iqtree", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("Best-fit model according to BIC:"):
                    model = line.split(":", 1)[1].strip()
                    break

    ok = res.returncode == 0 and os.path.exists(prefix + ".treefile")
    log(f"{gene} {cohort}: {'OK' if ok else 'FAILED'} n={info['n_taxa']} "
        f"len={info['length_nt']} model={model or '?'} {rt}s")
    if not ok:
        log(f"  stderr: {(res.stderr or '')[:300]}")

    return dict(**{k: info[k] for k in ("gene", "cohort", "n_taxa", "length_nt")},
                model=model, runtime_s=rt, status="OK" if ok else "FAILED",
                safe_kernel=used_safe)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prep", action="store_true", help="prepare alignments only")
    ap.add_argument("--genes", default="", help="comma-separated gene subset")
    args = ap.parse_args()

    wanted = [g for g in args.genes.split(",") if g] or None

    log("=" * 78)
    log("STEP 10b: HyPhy-ready alignments and per-type trees (combined cohorts excluded)")
    if not os.path.exists(IQTREE):
        log(f"ERROR: IQ-TREE not found at {IQTREE}")
        return 1

    prepped = []
    for cohort in COHORTS:
        for path in sorted(glob.glob(os.path.join(DIR_MASKED, cohort, "*_masked.fasta"))):
            gene = os.path.basename(path).replace("_masked.fasta", "")
            if wanted and gene not in wanted:
                continue
            info = prep(gene, cohort)
            if info:
                prepped.append(info)
    log(f"Prepared {len(prepped)} HyPhy-ready alignments (reference rows dropped)")

    if args.prep:
        return 0

    rows = [build_tree(info) for info in prepped]
    df = pd.DataFrame(rows)
    out = os.path.join(DIR_TREES, "tree_summary_v2.tsv")
    df.to_csv(out, sep="\t", index=False)
    log(f"Wrote {out}: {(df.status != 'FAILED').sum()}/{len(df)} trees available")
    return 0 if (df.status != "FAILED").all() else 1


if __name__ == "__main__":
    sys.exit(main())
