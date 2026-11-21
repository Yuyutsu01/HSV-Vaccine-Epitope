"""
Step 17 (CORRECTION A2 + A3 + A4): Selection analysis per type, breakpoint-aware, FDR-corrected.

Replaces step11_stage1_pipeline.py and step12_stage2_pipeline.py. Changes:
  - PER TYPE ONLY. Pooled HSV-1 + HSV-2 cohorts are not analysed: between-type amino-acid
    identity runs from 37% (gG) to 91% (UL30), and a single codon model across that scores
    type divergence as positive selection (CORRECTIONS_AND_PENDING_WORK.md, A2).
  - BREAKPOINT-AWARE. Where step16 found GARD breakpoints, each partition is analysed
    separately with its own tree. Previously nothing was ever split (A3).
  - ONE PASS, not a trigger-gated two-stage screen. The old design ran FEL/MEME only on
    genes whose uncorrected Stage 1 p-values looked interesting, which is itself a
    selection-on-significance step. All five methods now run on every unit.
  - FDR CORRECTION. Site-level p-values from FEL and MEME get Benjamini-Hochberg within
    each analysis unit; aBSREL's Holm-Bonferroni corrected branch p-values are used;
    BUSTED p-values get BH across the gene x cohort family. FUBAR is reported at posterior
    >= 0.95. Nothing in the old pipeline was corrected for multiple testing (A4).

Methods: BUSTED (gene-wide), FEL (pervasive), MEME (episodic), FUBAR (fast Bayesian),
aBSREL (branch-site).

Inputs:
  - 06_alignments_v2/hyphy_ready/<cohort>/<gene>.fasta        (step10b)
  - 06_alignments_v2/trees_v2/<cohort>/<gene>.treefile        (step10b)
  - 14_selection/partitions/<cohort>/<gene>_part<N>.fasta     (step16, when breakpoints exist)
  - 14_selection/gard_summary_per_type.tsv                    (step16)

Outputs:
  - 14_selection/hyphy/<cohort>/<unit>_<METHOD>.json
  - 14_selection/selection_runs.tsv          (one row per run, with provenance)
  - 14_selection/selection.log

Site-level results are assembled and FDR-corrected by step18.

Usage:
  python step17_selection_per_type.py
  python step17_selection_per_type.py --genes US6,UL27 --methods FEL,MEME
"""

import os
import re
import sys
import json
import glob
import argparse
import threading
import hashlib
import subprocess
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hyphy_runner import (DIR_BASE, HYPHY_IMAGE, HYPHY_CPUS, docker_available,
                          wait_for_docker, ensure_image, hyphy_provenance, rel)

DIR_READY = os.path.join(DIR_BASE, "06_alignments_v2", "hyphy_ready")
DIR_TREES = os.path.join(DIR_BASE, "06_alignments_v2", "trees_v2")
DIR_SEL = os.path.join(DIR_BASE, "14_selection")
DIR_PARTS = os.path.join(DIR_SEL, "partitions")
DIR_HYPHY = os.path.join(DIR_SEL, "hyphy")

COHORTS = ("HSV1", "HSV2")
METHODS = ("BUSTED", "FEL", "MEME", "FUBAR", "ABSREL")
GENE_ORDER = ["US6", "UL27", "UL22", "UL1", "UL44", "UL10", "UL53", "US4", "US7", "US8",
              "UL49A", "UL48", "UL54", "RL2", "RS1", "UL46", "UL47", "UL49", "UL23", "UL30"]

RUN_TIMEOUT = int(os.environ.get("HSV_SEL_TIMEOUT", "43200"))
LOG_PATH = os.path.join(DIR_SEL, "selection.log")


_LOG_LOCK = threading.Lock()


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    with _LOG_LOCK:                   # workers log concurrently
        print(line, flush=True)
        os.makedirs(DIR_SEL, exist_ok=True)
        with open(LOG_PATH, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def run_method(method, alignment, tree, output):
    """
    HyPhy call. aBSREL is spelled 'absrel' on the command line. FUBAR does not take a
    --tree argument in the same way as the likelihood methods but accepts one; passing the
    fixed tree keeps every method on the identical topology.
    """
    cmd = ["docker", "run", "--rm",
           "--mount", f"type=bind,src={DIR_BASE},dst=/data",
           "-w", "/data", HYPHY_IMAGE,
           "hyphy", "ENV=TOLERATE_NUMERICAL_ERRORS=1", f"CPU={HYPHY_CPUS}",
           method.lower(),
           "--alignment", rel(alignment),
           "--tree", rel(tree),
           "--output", rel(output)]
    t0 = datetime.now()
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=RUN_TIMEOUT)
    except subprocess.TimeoutExpired:
        return False, (datetime.now() - t0).total_seconds(), f"timeout after {RUN_TIMEOUT}s"
    rt = round((datetime.now() - t0).total_seconds(), 1)

    if res.returncode != 0:
        tail = (res.stdout or "")[-600:] + (res.stderr or "")[-300:]
        msg = next((l for l in reversed(tail.splitlines()) if l.strip()), "unknown error")
        return False, rt, msg[:300]
    if not os.path.exists(output):
        return False, rt, "exited 0 but wrote no output"
    try:
        with open(output, encoding="utf-8") as fh:
            json.load(fh)
    except Exception:
        return False, rt, "output is not valid JSON"
    return True, rt, "ok"


def gard_breakpoints():
    """{(gene, cohort): n_breakpoints} from step16."""
    path = os.path.join(DIR_SEL, "gard_summary_per_type.tsv")
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    return {(r.gene, r.cohort): int(r.num_breakpoints or 0)
            for r in df.itertuples() if str(getattr(r, "status", "")) != "FAILED"}


def partition_tree(part_fasta, prefix):
    """Build a tree for one GARD partition (each segment needs its own topology)."""
    iqtree = os.path.join(DIR_BASE, "bin", "iqtree",
                          "iqtree-2.4.0-Windows", "bin", "iqtree2.exe")
    tf = prefix + ".treefile"
    # Reuse a tree only if it is COMPLETE (the .iqtree report exists) and describes THIS
    # alignment. Partition spans change when the breakpoints change and `prefix` is
    # reused across runs, so an existence check alone silently pairs a new partition
    # with the old one's tree. Keyed on the alignment's content hash rather than its
    # mtime, because re-exporting identical partitions touches every file.
    with open(part_fasta, "rb") as fh:
        aln_md5 = hashlib.md5(fh.read()).hexdigest()
    stamp = prefix + ".alnmd5"
    if os.path.exists(tf) and os.path.exists(prefix + ".iqtree"):
        prev = None
        if os.path.exists(stamp):
            with open(stamp, encoding="utf-8") as fh:
                prev = fh.read().strip()
        if prev == aln_md5:
            return tf
        log(f"{os.path.basename(prefix)}: tree does not match this partition alignment "
            f"({'no hash recorded' if prev is None else 'alignment changed'}); rebuilding")
    cmd = [iqtree, "-s", part_fasta, "-m", "GTR+F+G4", "--keep-ident",
           "-T", "4", "--prefix", prefix, "--quiet", "-redo"]
    subprocess.run(cmd, capture_output=True, text=True)
    ok = os.path.exists(tf) and os.path.exists(prefix + ".iqtree")
    if ok:
        with open(stamp, "w", encoding="utf-8") as fh:
            fh.write(aln_md5 + "\n")
    return tf if ok else None


def units(genes, cohorts, bps, skip=frozenset()):
    """
    One analysis unit per gene x cohort, split into GARD partitions where breakpoints exist.
    Returns dicts with the alignment, the tree and a label.
    """
    out = []
    for gene in genes:
        for cohort in cohorts:
            # `skip` holds (gene, cohort) PAIRS. Filtering by separate gene and cohort
            # lists does not work: dropping RS1/HSV1 while RS1 survives via HSV2 and
            # HSV1 survives via other genes lets this cross-product rebuild the very
            # pair that was excluded, and it would then be analysed unpartitioned.
            if (gene, cohort) in skip:
                continue
            aln = os.path.join(DIR_READY, cohort, f"{gene}.fasta")
            tree = os.path.join(DIR_TREES, cohort, f"{gene}.treefile")
            if not os.path.exists(aln):
                continue

            nbp = bps.get((gene, cohort), 0)
            parts = sorted(glob.glob(os.path.join(DIR_PARTS, cohort, f"{gene}_part*.fasta")))

            if nbp and parts:
                for p in parts:
                    n = re.search(r"_part(\d+)\.fasta$", p).group(1)
                    prefix = os.path.join(DIR_PARTS, cohort, f"{gene}_part{n}")
                    pt = partition_tree(p, prefix)
                    if not pt:
                        log(f"{gene} {cohort} part{n}: no tree could be built, skipped")
                        continue
                    out.append(dict(gene=gene, cohort=cohort, unit=f"{gene}_part{n}",
                                    partition=int(n), n_breakpoints=nbp,
                                    alignment=p, tree=pt))
            else:
                # The .iqtree report is written only when IQ-TREE completes. A bare
                # .treefile is the ML tree from a run that died during bootstrapping,
                # carrying no support values, and must not be analysed as if it were
                # a finished tree. main() refuses to start when any are missing.
                if not (os.path.exists(tree)
                        and os.path.exists(tree[:-len(".treefile")] + ".iqtree")):
                    log(f"{gene} {cohort}: no completed tree, skipped "
                        f"(re-run after step10b finishes)")
                    continue
                out.append(dict(gene=gene, cohort=cohort, unit=gene, partition=0,
                                n_breakpoints=nbp, alignment=aln, tree=tree))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="")
    ap.add_argument("--cohorts", default=",".join(COHORTS))
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--wait-docker", type=int, default=300)
    # Off by default: a partial GARD screen silently mixes partition-aware and
    # unpartitioned results in one output set.
    ap.add_argument("--allow-partial-gard", action="store_true")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--allow-missing-trees", action="store_true")
    args = ap.parse_args()

    genes = [g for g in args.genes.split(",") if g] or GENE_ORDER
    cohorts = [c for c in args.cohorts.split(",") if c]
    methods = [m.upper() for m in args.methods.split(",") if m]

    os.makedirs(DIR_SEL, exist_ok=True)
    log("=" * 78)
    log("STEP 17: selection analysis per type, breakpoint-aware "
        f"({', '.join(methods)})")

    ok, msg = docker_available()
    if not ok:
        log(f"Docker not ready ({msg}); waiting up to {args.wait_docker}s")
        ok, msg = wait_for_docker(max_wait=args.wait_docker)
    if not ok:
        log(f"ABORT: HyPhy needs the Docker daemon. {msg}")
        return 2
    log(msg)
    ok, msg = ensure_image()
    if not ok:
        log(f"ABORT: {msg}")
        return 2
    prov = hyphy_provenance()
    log(f"HyPhy provenance: {prov['version']} | CPU={prov['cpus']}")

    bps = gard_breakpoints()
    if not bps:
        log("ABORT: no GARD summary found. Selection must follow recombination screening "
            "(step16); running unpartitioned is the error correction A6/A7 fixes.")
        return 2

    # Selection run on a partial GARD screen mixes partition-aware results for the genes
    # GARD finished with unpartitioned results for the genes it had not reached, in one
    # output set, with nothing in the files to tell them apart. Require a complete screen.
    expected = {(g, c) for g in genes for c in cohorts
                if os.path.exists(os.path.join(DIR_READY, c, f"{g}.fasta"))}
    missing = sorted(expected - set(bps))
    if missing and not args.allow_partial_gard:
        log(f"ABORT: GARD screen is incomplete -- {len(missing)} of {len(expected)} "
            f"unit(s) have no result yet: {', '.join(f'{g}/{c}' for g, c in missing[:12])}"
            + (" ..." if len(missing) > 12 else ""))
        log("Wait for step16 to finish, or pass --allow-partial-gard to analyse only "
            "the units GARD has completed (the summary then records which).")
        return 2
    skip = set()
    if missing:
        log(f"WARNING: --allow-partial-gard set; skipping {len(missing)} unit(s) with no "
            f"GARD result rather than running them unpartitioned: "
            + ", ".join(f"{g}/{c}" for g, c in missing))
        skip = set(missing)

    # A gene whose tree is missing or unfinished was previously skipped with only a log
    # line, so the result set would quietly omit it and nothing in the output would say
    # so. Treat it the same as an incomplete GARD screen: refuse to start.
    no_tree = []
    for gene, cohort in sorted(expected - skip):
        if bps.get((gene, cohort), 0) and glob.glob(
                os.path.join(DIR_PARTS, cohort, f"{gene}_part*.fasta")):
            continue                      # partitioned: gets its own partition trees
        tf = os.path.join(DIR_TREES, cohort, f"{gene}.treefile")
        if not (os.path.exists(tf)
                and os.path.exists(tf[:-len(".treefile")] + ".iqtree")):
            no_tree.append((gene, cohort))
    if no_tree and not args.allow_missing_trees:
        log(f"ABORT: {len(no_tree)} unit(s) have no completed tree: "
            + ", ".join(f"{g}/{c}" for g, c in no_tree[:12])
            + (" ..." if len(no_tree) > 12 else ""))
        log("A .treefile without its .iqtree report is an unbootstrapped tree from an "
            "interrupted run. Wait for step10b, or pass --allow-missing-trees to "
            "analyse only the units that have one.")
        return 2
    if no_tree:
        log(f"WARNING: --allow-missing-trees set; {len(no_tree)} unit(s) will be absent "
            f"from this result set: {', '.join(f'{g}/{c}' for g, c in no_tree)}")

    work = units(genes, cohorts, bps, skip=skip)
    # Belt and braces: assert no skipped pair survived into the work list.
    leaked = sorted({(u["gene"], u["cohort"]) for u in work} & skip)
    if leaked:
        log(f"ABORT: internal error -- excluded unit(s) present in the work list: "
            + ", ".join(f"{g}/{c}" for g, c in leaked))
        return 3
    log(f"{len(work)} analysis unit(s) x {len(methods)} method(s) = "
        f"{len(work) * len(methods)} runs")

    rows = []
    rows_lock = threading.Lock()

    def record(row):
        with rows_lock:
            rows.append(row)
            pd.DataFrame(rows).to_csv(
                os.path.join(DIR_SEL, "selection_runs.tsv"), sep="\t", index=False)

    def run_one(u, method):
        out_dir = os.path.join(DIR_HYPHY, u["cohort"])
        os.makedirs(out_dir, exist_ok=True)
        out_json = os.path.join(out_dir, f"{u['unit']}_{method}.json")
        # Same completion-marker rule as step16: a JSON that merely parses may be a
        # partial write from a run that was killed, and caching on that silently turns
        # an interrupted analysis into a result.
        done_marker = out_json + ".done"

        cached = False
        if os.path.exists(out_json):
            try:
                with open(out_json, encoding="utf-8") as fh:
                    json.load(fh)
            except Exception:
                os.remove(out_json)
            else:
                if os.path.exists(done_marker):
                    cached = True
                else:
                    log(f"{u['unit']} {u['cohort']} {method}: previous output has no "
                        f"completion marker; discarding and re-running")
                    os.remove(out_json)

        if cached:
            status, rt, m = "CACHED", 0.0, "ok"
        else:
            # Defensive: drop any HyPhy sidecar cache older than the alignment it sits
            # beside. FUBAR keys its grid cache on the filename alone, so a cache built
            # for a different partition span is silently reloaded and the run dies on a
            # dimension mismatch. step16 clears these when spans change; this catches
            # anything that slipped through.
            for side in glob.glob(u["alignment"] + ".*cache"):
                try:
                    if os.path.getmtime(side) < os.path.getmtime(u["alignment"]):
                        os.remove(side)
                        log(f"{u['unit']} {u['cohort']}: removed stale "
                            f"{os.path.basename(side)} (older than its alignment)")
                except OSError:
                    pass
            good, rt, m = run_method(method, u["alignment"], u["tree"], out_json)
            status = "OK" if good else "FAILED"
            if good:
                with open(done_marker, "w", encoding="utf-8") as fh:
                    fh.write(f"completed in {rt}s\n")
            else:
                log(f"{u['unit']} {u['cohort']} {method}: FAILED {rt}s -- {m}")

        if status != "FAILED":
            log(f"{u['unit']} {u['cohort']} {method}: {status} {rt}s")

        record(dict(unit=u["unit"], gene=u["gene"], cohort=u["cohort"],
                    partition=u["partition"], n_breakpoints=u["n_breakpoints"],
                    method=method, status=status, runtime_s=rt, message=m,
                    alignment=os.path.relpath(u["alignment"], DIR_BASE),
                    tree=os.path.relpath(u["tree"], DIR_BASE),
                    output=os.path.relpath(out_json, DIR_BASE),
                    hyphy_version=prov["version"]))

    def guarded(u, method):
        """One run's failure must not take the other workers down with it."""
        try:
            run_one(u, method)
        except Exception as exc:
            log(f"{u['unit']} {u['cohort']} {method}: ERROR {type(exc).__name__}: {exc}")
            record(dict(unit=u["unit"], gene=u["gene"], cohort=u["cohort"],
                        partition=u["partition"], n_breakpoints=u["n_breakpoints"],
                        method=method, status="FAILED", runtime_s=0.0,
                        message=f"{type(exc).__name__}: {exc}",
                        hyphy_version=prov["version"]))

    # Every (unit, method) pair is an independent HyPhy analysis writing to its own
    # file, so they parallelise freely. HyPhy's own threading plateaus near 1.5 cores,
    # so throughput comes from running several at once with a low CPU= each.
    jobs = [(u, m) for u in work for m in methods]
    if args.workers > 1:
        log(f"Running {len(jobs)} analyses on {args.workers} workers "
            f"(HyPhy CPU={HYPHY_CPUS} each)")
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(lambda j: guarded(*j), jobs))
    else:
        for u, m in jobs:
            guarded(u, m)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DIR_SEL, "selection_runs.tsv"), sep="\t", index=False)
    log(f"Selection complete: {(df.status != 'FAILED').sum()}/{len(df)} runs available")
    log("Run step18 to assemble site-level results with Benjamini-Hochberg FDR.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
