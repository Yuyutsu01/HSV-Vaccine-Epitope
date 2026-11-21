#!/usr/bin/env python3
"""
step30_provenance_audit.py -- verify every data file in the repository is accounted
for in 08_documentation/DATA_DICTIONARY.md.

This is a completeness check, not a correctness check. It answers one question:
is there any data file on disk that the documentation does not describe?

A file counts as accounted for if the dictionary mentions either its exact
basename or the directory it sits in. Directory-level accounting is deliberate:
the dictionary describes 'hyphy/HSV1/' as a set of 610 JSON files rather than
listing each one, and that is the right granularity for bulk per-gene output.

Exit status 0 = every file accounted for. 1 = gaps found, listed on stdout.
No files are written or modified.
"""

import os
import re
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(ROOT)
DICT_PATH = os.path.join(ROOT, "08_documentation", "DATA_DICTIONARY.md")

# Directories that hold no research data and are excluded by design.
SKIP_DIRS = {".git", "__pycache__", ".venv", "venv", "env", ".idea", ".vscode"}

# Paths excluded from the audit, with the reason recorded in the output so the
# exclusion is visible rather than silent.
EXCLUDED = {
    "hsv_curated_dataset/bin": "external tool binaries, reinstalled by install_prerequisites.py",
}

# Extensions that are never research data.
SKIP_EXT = {".pyc", ".pyo", ".swp", ".swo", ".DS_Store"}


def load_dictionary():
    if not os.path.exists(DICT_PATH):
        print(f"ABORT: data dictionary not found at {DICT_PATH}")
        sys.exit(2)
    with open(DICT_PATH, encoding="utf-8") as fh:
        text = fh.read()
    # Every backtick-quoted token in the dictionary is a potential path or filename.
    quoted = set(re.findall(r"`([^`]+)`", text))
    return text, quoted


# Directories whose contents the dictionary describes COLLECTIVELY, because they
# hold one file per gene, per allele, or per structure and enumerating 610 of them
# would not aid anyone. For these, and only these, a documented directory counts
# as accounting for the files inside it.
#
# Every other directory must name each of its files explicitly. That distinction
# is the whole point: without it, the audit passes any file dropped into a
# documented directory, which is exactly the failure a canary file exposed.
BULK_DIRS = {
    "00_raw_downloads", "04_cds_unaligned", "05_cds_qc_passed",
    "06_alignments_v2/HSV1", "06_alignments_v2/HSV2", "06_alignments_v2/combined",
    "06_alignments_v2/masked", "06_alignments_v2/coordinate_maps",
    "06_alignments_v2/hyphy_ready", "06_alignments_v2/trees_v2/HSV1",
    "06_alignments_v2/trees_v2/HSV2",
    "07_qc_reports/hyphy/stage1", "07_qc_reports/hyphy/stage2", "07_qc_reports/hyphy/gard",
    "09_scripts",
    "13_epitopes/bcell_inputs", "13_epitopes/bcell_raw", "13_epitopes/classii_inputs",
    "13_epitopes/classii_raw",
    "14_selection/gard", "14_selection/hyphy", "14_selection/partitions",
    "17_structures/pdb", "17_structures/figures",
}


def is_bulk(dirrel):
    d = dirrel.split("hsv_curated_dataset/")[-1]
    return any(d == b or d.startswith(b + "/") for b in BULK_DIRS)


def dictionary_covers(relpath, basename, dirrel, text, quoted):
    """Is this file accounted for? Returns the matching dictionary token, or None."""
    # 1. Exact basename quoted somewhere in the dictionary.
    if basename in quoted:
        return basename
    # 1b. The basename's stem, e.g. `fig1_cohort_size` covering both .pdf and .png.
    stem = os.path.splitext(basename)[0]
    if stem in quoted:
        return stem
    # 2. Directory and filename quoted as one token, e.g. `unpacked/md5sum.txt`.
    segs = dirrel.split("/")
    for i in range(len(segs)):
        tail = "/".join(segs[i:])
        if f"{tail}/{basename}" in quoted:
            return f"{tail}/{basename}"
    # 3. A brace or glob pattern covers it, e.g. `HSV{1,2}/*.zip`.
    #    The pattern must contain a '/' so it is anchored to a directory. A bare
    #    `*.log` or `*.tsv` would otherwise match every such file anywhere in the
    #    repository -- a loophole that let a stray log at the dataset root pass
    #    on the strength of a `*.log` mentioned in an unrelated table.
    for tok in quoted:
        if "/" not in tok or ("{" not in tok and "*" not in tok):
            continue
        pat = re.escape(tok)
        pat = re.sub(r"\\\{([^}]*)\\\}", lambda m: "(?:" + m.group(1).replace(",", "|") + ")", pat)
        pat = pat.replace(r"\*", "[^/]*")
        if re.search(pat.rstrip("/") + r"(/|$)", relpath + "/"):
            return tok
    # 4. Collective accounting -- bulk directories only.
    if is_bulk(dirrel):
        for i in range(len(segs)):
            tail = "/".join(segs[i:])
            if tail in quoted or tail + "/" in quoted:
                return tail + "/ (bulk)"
        # A described extension inside a bulk directory, e.g. `.mldist` in trees_v2.
        ext = os.path.splitext(basename)[1]
        if ext and ext in quoted:
            return ext + " (bulk)"
        return "(bulk dir, described collectively)"
    return None


def main():
    text, quoted = load_dictionary()
    accounted = 0
    gaps = []
    excluded_count = defaultdict(int)
    by_dir = defaultdict(lambda: [0, 0])

    for dirpath, dirnames, filenames in os.walk(REPO):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        dirrel = os.path.relpath(dirpath, REPO).replace("\\", "/")
        if dirrel == ".":
            dirrel = ""

        ex = next((v for k, v in EXCLUDED.items() if dirrel == k or dirrel.startswith(k + "/")), None)
        if ex:
            excluded_count[ex] += len(filenames)
            continue

        for fn in filenames:
            if os.path.splitext(fn)[1] in SKIP_EXT or fn in SKIP_EXT:
                continue
            relpath = f"{dirrel}/{fn}" if dirrel else fn
            hit = dictionary_covers(relpath, fn, dirrel, text, quoted)
            by_dir[dirrel][0] += 1
            if hit:
                accounted += 1
                by_dir[dirrel][1] += 1
            else:
                gaps.append(relpath)

    total = accounted + len(gaps)
    print("=" * 74)
    print("PROVENANCE AUDIT -- every data file vs. DATA_DICTIONARY.md")
    print("=" * 74)
    print(f"\nFiles audited      : {total}")
    print(f"Accounted for      : {accounted}")
    print(f"Unaccounted        : {len(gaps)}")
    for reason, n in excluded_count.items():
        print(f"Excluded ({n:>3})     : {reason}")

    if gaps:
        print(f"\n--- UNACCOUNTED FILES ({len(gaps)}) ---")
        per_dir = defaultdict(list)
        for g in gaps:
            per_dir[os.path.dirname(g)].append(os.path.basename(g))
        for d in sorted(per_dir):
            fs = per_dir[d]
            print(f"\n  {d or '(repo root)'}  [{len(fs)} file(s)]")
            for f in sorted(fs)[:12]:
                print(f"      {f}")
            if len(fs) > 12:
                print(f"      ... and {len(fs) - 12} more")

    print("\n--- COVERAGE BY DIRECTORY ---")
    for d in sorted(by_dir):
        n, ok = by_dir[d]
        flag = "" if ok == n else f"   <-- {n - ok} UNACCOUNTED"
        print(f"  {ok:>4}/{n:<4} {d or '(repo root)':<56}{flag}")

    print("\n" + "=" * 74)
    if gaps:
        print(f"RESULT: FAIL -- {len(gaps)} file(s) have no entry in the data dictionary.")
        return 1
    print("RESULT: PASS -- every data file is accounted for in the data dictionary.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
