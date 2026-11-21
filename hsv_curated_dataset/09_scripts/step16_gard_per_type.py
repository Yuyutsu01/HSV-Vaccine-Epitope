"""
Step 16 (CORRECTION A2 + A3): Recombination screening PER TYPE, before selection, with
breakpoint partitions exported for the selection stage.

What changed from step12_gard_pipeline.py:
  - GARD is run on HSV-1 and HSV-2 cohorts separately. The previous runs used pooled
    HSV-1 + HSV-2 alignments, where GARD detects the species split rather than
    recombination; that produced all four "RECOMBINANT" verdicts in the old README
    (ΔAICc up to 6,619 for UL30) (CORRECTIONS_AND_PENDING_WORK.md, A2).
  - This stage now runs BEFORE selection, and writes each breakpoint partition out as its
    own in-frame FASTA so the selection stage can analyse segments separately, as the
    proposal requires. Previously GARD ran last and no alignment was ever split (A3).
  - All 20 genes are screened, not a top-10 subset, because every gene that goes into a
    selection test needs its own recombination screen.

Inputs:  06_alignments_v2/hyphy_ready/<cohort>/<gene>.fasta   (step10b)
Outputs:
  - 14_selection/gard/<cohort>/<gene>_GARD.json
  - 14_selection/partitions/<cohort>/<gene>_part<N>.fasta     (only where breakpoints found)
  - 14_selection/gard_summary_per_type.tsv
  - 14_selection/gard_provenance.json
  - 14_selection/gard.log

Usage:
  python step16_gard_per_type.py                       # all genes, glycoproteins first
  python step16_gard_per_type.py --genes US6,UL27
  python step16_gard_per_type.py --cohorts HSV1
"""

import os
import re
import sys
import json
import glob
import argparse
import threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hyphy_runner import (DIR_BASE, docker_available, wait_for_docker, ensure_image,
                          hyphy_provenance, run_hyphy)

DIR_READY = os.path.join(DIR_BASE, "06_alignments_v2", "hyphy_ready")
DIR_OUT = os.path.join(DIR_BASE, "14_selection")
DIR_GARD = os.path.join(DIR_OUT, "gard")
DIR_PARTS = os.path.join(DIR_OUT, "partitions")

COHORTS = ("HSV1", "HSV2")

# Glycoproteins and the primary vaccine antigens first, so the headline genes finish early.
GENE_ORDER = ["US6", "UL27", "UL22", "UL1", "UL44", "UL10", "UL53", "US4", "US7", "US8",
              "UL49A", "UL48", "UL54", "RL2", "RS1", "UL46", "UL47", "UL49", "UL23", "UL30"]

# No per-run ceiling by default. GARD's cost grows with each breakpoint stage it adds,
# so the genes that take longest are precisely the ones with the most recombination to
# find -- exactly the results worth waiting for. A timeout here would discard hours of
# computation and leave that gene with no result at all, which is worse than waiting.
# Set HSV_GARD_TIMEOUT to a number of seconds to impose one.
RUN_TIMEOUT = int(os.environ.get("HSV_GARD_TIMEOUT", "0")) or None

LOG_PATH = os.path.join(DIR_OUT, "gard.log")


_LOG_LOCK = threading.Lock()


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    with _LOG_LOCK:                   # workers log concurrently
        print(line, flush=True)
        os.makedirs(DIR_OUT, exist_ok=True)
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


def parse_gard(path):
    """Breakpoints and model fit from a GARD JSON."""
    with open(path, encoding="utf-8") as fh:
        js = json.load(fh)

    base = js.get("baselineScore")
    best = js.get("bestModelAICc")

    # breakpointData maps SEGMENT index -> {"bps": [[start, end]]} -- one entry per
    # segment of the selected partitioning, each carrying only its own span. Reading
    # just the highest-index entry therefore yields the LAST segment and mistakes its
    # start for the only breakpoint: RS1/HSV2's five segments came out as one
    # breakpoint. Collect every segment and take the interior boundaries.
    sites = js.get("input", {}).get("number of sites")
    segments = []
    bd = js.get("breakpointData") or {}
    if isinstance(bd, dict):
        for key in sorted(bd, key=lambda k: int(k) if str(k).isdigit() else -1):
            entry = bd[key]
            raw = entry.get("bps") if isinstance(entry, dict) else entry
            for span in (raw or []):
                if isinstance(span, list) and len(span) >= 2:
                    segments.append((int(span[0]), int(span[1])))
                elif isinstance(span, (int, float)):
                    segments.append((int(span), int(span)))

    segments.sort()
    # Each segment after the first begins at a breakpoint. Sentinel edges (site 1 and
    # the final site) are segment bounds, not breakpoints, so they fall out naturally.
    bps = sorted({s for s, _e in segments[1:]
                  if 0 < s < (sites or 10 ** 9)})

    delta = (base - best) if isinstance(base, (int, float)) and isinstance(best, (int, float)) else None
    return dict(baseline_aicc=round(base, 2) if isinstance(base, (int, float)) else None,
                best_aicc=round(best, 2) if isinstance(best, (int, float)) else None,
                single_tree_aicc=js.get("singleTreeAICc"),
                delta_aicc=round(delta, 2) if isinstance(delta, (int, float)) else None,
                breakpoints=bps, sites=sites,
                sequences=js.get("input", {}).get("number of sequences"),
                models_evaluated=js.get("totalModelCount"))


def clear_stale_partitions(cohort, gene):
    """
    Remove any partition FASTA this gene x cohort left behind earlier.

    step17 picks partitions up with a glob, so a stale file is silently analysed as
    though the current GARD run had produced it. That happens whenever a re-run finds
    FEWER breakpoints than before (orphaned high-numbered parts) or none at all (the
    export is skipped entirely, so every old part survives). Clearing first makes the
    partition directory a statement about the current run only.
    """
    stale = glob.glob(os.path.join(DIR_PARTS, cohort, f"{gene}_part*.fasta"))
    for p in stale:
        os.remove(p)
    # HyPhy keeps its own sidecars next to the alignment, keyed on the FILENAME only --
    # FUBAR's grid cache above all. A cache left from a different span is reloaded for
    # the new partition and FUBAR dies on a dimension mismatch ("Had 400x159 and
    # 400x676"), so the sidecars must go whenever the spans change.
    sidecars = [p for pat in ("_part*.fasta.*cache", "_part*.fasta.*-gard*",
                              "_part*.fasta.*.bf")
                for p in glob.glob(os.path.join(DIR_PARTS, cohort, f"{gene}{pat}"))]
    for p in sidecars:
        os.remove(p)
    if stale or sidecars:
        log(f"{gene} {cohort}: cleared {len(stale)} partition file(s) and "
            f"{len(sidecars)} HyPhy sidecar(s) from a previous run")
    return len(stale)


def export_partitions(aln_path, breakpoints, cohort, gene):
    """
    Write one FASTA per GARD segment. Breakpoints are nucleotide positions; they are snapped
    outward to codon boundaries so every partition stays in frame for codon models.
    """
    ids, seqs = read_fasta(aln_path)
    if not seqs:
        return []
    L = len(seqs[0])

    cuts = [0]
    for b in breakpoints:
        c = max(3, min(L - 3, (int(b) // 3) * 3))
        if c not in cuts:
            cuts.append(c)
    cuts.append(L)
    cuts = sorted(set(cuts))

    out_dir = os.path.join(DIR_PARTS, cohort)
    os.makedirs(out_dir, exist_ok=True)

    written = []
    for n, (s, e) in enumerate(zip(cuts[:-1], cuts[1:]), start=1):
        if e - s < 30:                       # too short for a codon model
            continue
        path = os.path.join(out_dir, f"{gene}_part{n}.fasta")
        with open(path, "w", encoding="utf-8") as fh:
            for i, seq in zip(ids, seqs):
                fh.write(f">{i}\n{seq[s:e]}\n")
        written.append(dict(partition=n, start_nt=s + 1, end_nt=e,
                            length_nt=e - s, path=path))
    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="")
    ap.add_argument("--cohorts", default=",".join(COHORTS))
    ap.add_argument("--wait-docker", type=int, default=300)
    # HyPhy's OpenMP build does not scale past roughly 1.5 cores on alignments this
    # small, so raising CPU= alone leaves most of the allocation idle. The genes are
    # independent analyses writing to separate files, so running several at once is
    # pure throughput with no effect on any result. Set HSV_HYPHY_CPUS low (2) and
    # --workers to fill the Docker CPU allocation.
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    genes = [g for g in args.genes.split(",") if g] or GENE_ORDER
    cohorts = [c for c in args.cohorts.split(",") if c]

    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 16: GARD recombination screening, PER TYPE, before selection")

    ok, msg = docker_available()
    if not ok:
        log(f"Docker not ready ({msg}); waiting up to {args.wait_docker}s")
        ok, msg = wait_for_docker(max_wait=args.wait_docker)
    if not ok:
        log(f"ABORT: HyPhy needs the Docker daemon. {msg}")
        log("Start Docker Desktop and re-run this script; completed runs are cached.")
        return 2
    log(msg)

    ok, msg = ensure_image()
    if not ok:
        log(f"ABORT: HyPhy image unavailable: {msg}")
        return 2
    prov = hyphy_provenance()
    log(f"HyPhy provenance: {prov['version']} | image {prov['image']} | CPU={prov['cpus']}")
    with open(os.path.join(DIR_OUT, "gard_provenance.json"), "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=2)

    rows = []
    rows_lock = threading.Lock()

    def record(row):
        with rows_lock:
            rows.append(row)
            pd.DataFrame(rows).to_csv(
                os.path.join(DIR_OUT, "gard_summary_per_type.tsv"), sep="\t", index=False)

    def run_unit(gene, cohort):
            aln = os.path.join(DIR_READY, cohort, f"{gene}.fasta")
            if not os.path.exists(aln):
                log(f"{gene} {cohort}: no HyPhy-ready alignment, skipped")
                return

            ids, seqs = read_fasta(aln)
            out_dir = os.path.join(DIR_GARD, cohort)
            os.makedirs(out_dir, exist_ok=True)
            out_json = os.path.join(out_dir, f"{gene}_GARD.json")

            # GARD rewrites its JSON after every breakpoint it adds, so an interrupted
            # run leaves a file that parses and carries every expected key but reports
            # only the breakpoints found so far. Accepting that as a cache hit silently
            # turns a killed run into a truncated result, so a run counts as complete
            # only if run_hyphy returned ok and wrote the marker beside the JSON.
            done_marker = out_json + ".done"
            cached = False
            if os.path.exists(out_json):
                try:
                    with open(out_json, encoding="utf-8") as fh:
                        json.load(fh)
                except Exception:
                    log(f"{gene} {cohort}: discarding unreadable previous output")
                    os.remove(out_json)
                else:
                    if os.path.exists(done_marker):
                        cached = True
                    else:
                        log(f"{gene} {cohort}: previous output has no completion marker "
                            f"(interrupted mid-search); discarding and re-running")
                        os.remove(out_json)

            # A marker saying TRUNCATED means the search was stopped before it reached
            # its own stopping criterion, so the partitioning is GARD's best model at
            # that point but is not established as exhaustive. Surface it as a column
            # rather than leaving it only in a sidecar file.
            truncated = False
            if os.path.exists(done_marker):
                try:
                    with open(done_marker, encoding="utf-8") as fh:
                        truncated = "TRUNCATED" in fh.read(200).upper()
                except Exception:
                    pass

            if cached:
                log(f"{gene} {cohort}: cached"
                    + (" (TRUNCATED SEARCH)" if truncated else ""))
                status, rt = ("CACHED_TRUNCATED_SEARCH" if truncated else "CACHED"), 0.0
            else:
                log(f"{gene} {cohort}: GARD start (n={len(ids)}, {len(seqs[0])} nt)")
                good, rt, m = run_hyphy("gard", aln, out_json, timeout=RUN_TIMEOUT)
                status = "OK" if good else "FAILED"
                if good:
                    with open(done_marker, "w", encoding="utf-8") as fh:
                        fh.write(f"completed in {rt}s\n")
                if not good:
                    log(f"{gene} {cohort}: FAILED after {rt}s -- {m}")
                    record(dict(gene=gene, cohort=cohort, status=status,
                                n_sequences=len(ids), sites_nt=len(seqs[0]),
                                runtime_s=rt, error=m))
                    return

            info = parse_gard(out_json)
            # Clear first, unconditionally: a gene that previously found breakpoints and
            # now finds none would otherwise keep its old partitions, and step17 would
            # analyse a gene as partitioned that this run says is not.
            clear_stale_partitions(cohort, gene)
            parts = (export_partitions(aln, info["breakpoints"], cohort, gene)
                     if info["breakpoints"] else [])

            verdict = ("NO BREAKPOINT DETECTED" if not info["breakpoints"]
                       else f"{len(info['breakpoints'])} BREAKPOINT(S)")
            record(dict(
                gene=gene, cohort=cohort, status=status, n_sequences=len(ids),
                sites_nt=len(seqs[0]), runtime_s=rt,
                baseline_aicc=info["baseline_aicc"], best_aicc=info["best_aicc"],
                single_tree_aicc=info["single_tree_aicc"],
                delta_aicc=info["delta_aicc"], models_evaluated=info["models_evaluated"],
                num_breakpoints=len(info["breakpoints"]),
                breakpoint_positions=";".join(str(b) for b in info["breakpoints"]) or "None",
                n_partitions_exported=len(parts),
                partition_ranges=";".join(f"{p['start_nt']}-{p['end_nt']}" for p in parts) or "None",
                search_truncated=truncated,
                verdict=verdict, error=""))
            log(f"{gene} {cohort}: {verdict} (dAICc={info['delta_aicc']}), "
                f"{len(parts)} partition(s) exported, {rt}s")

    def guarded(gene, cohort):
        """One unit's failure must not take the other workers down with it."""
        try:
            run_unit(gene, cohort)
        except Exception as exc:
            log(f"{gene} {cohort}: ERROR {type(exc).__name__}: {exc}")
            record(dict(gene=gene, cohort=cohort, status="FAILED",
                        runtime_s=0.0, error=f"{type(exc).__name__}: {exc}"))

    units = [(g, c) for g in genes for c in cohorts]
    if args.workers > 1:
        log(f"Running {len(units)} units on {args.workers} workers "
            f"(HyPhy CPU={os.environ.get('HSV_HYPHY_CPUS', '4')} each)")
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(lambda u: guarded(*u), units))
    else:
        for g, c in units:
            guarded(g, c)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DIR_OUT, "gard_summary_per_type.tsv"), sep="\t", index=False)
    done = (df.status != "FAILED").sum() if len(df) else 0
    log(f"GARD complete: {done}/{len(df)} runs available")
    log("Note: no breakpoint detected is not evidence of clonality -- report n per cohort.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
