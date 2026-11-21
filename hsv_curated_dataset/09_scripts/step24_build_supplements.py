"""
Step 24 (CORRECTION A10): Build the publication supplements.

Curation step 0 asked for a count of sequences surviving each filter, and curation step 11 for
the final accession list and a table of removed sequences with the reason for each removal.
Neither existed. This step assembles them from the QC logs that the pipeline already writes,
so the numbers in the methods section are generated rather than transcribed.

Outputs:
  - 19_supplements/S1_data_flow.tsv              sequences surviving each filter
  - 19_supplements/S2_final_accessions.tsv       accessions in the analysed cohorts
  - 19_supplements/S3_removed_sequences.tsv      every exclusion, with its reason
  - 19_supplements/S4_per_gene_cohorts.tsv       n per gene x type at each stage
  - 19_supplements/S5_software_versions.tsv      tool provenance
  - 19_supplements/supplements.log
"""

import os
import sys
import glob
import json
import subprocess
from datetime import datetime

import pandas as pd

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_QC = os.path.join(DIR_BASE, "07_qc_reports")
DIR_META = os.path.join(DIR_BASE, "03_metadata")
DIR_V2 = os.path.join(DIR_BASE, "06_alignments_v2")
DIR_READY = os.path.join(DIR_V2, "hyphy_ready")
DIR_OUT = os.path.join(DIR_BASE, "19_supplements")

LOG_PATH = os.path.join(DIR_OUT, "supplements.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def read_tsv(path):
    return pd.read_csv(path, sep="\t") if os.path.exists(path) else None


def data_flow():
    rows = []

    meta = read_tsv(os.path.join(DIR_META, "raw_metadata_all_genomes.tsv"))
    curated = read_tsv(os.path.join(DIR_META, "curated_metadata.tsv"))
    filt = read_tsv(os.path.join(DIR_QC, "04_genome_filtering_log.tsv"))

    if meta is not None:
        rows.append(dict(stage="1. Records retrieved from NCBI Virus",
                         level="genome", n=len(meta), detail="HSV-1 taxid 10298 + HSV-2 taxid 10310"))
    if curated is not None:
        rows.append(dict(stage="2. Metadata curated", level="genome", n=len(curated),
                         detail="one row per genome record"))
    if filt is not None:
        status = filt.filter_status.astype(str)
        passed = int((status == "PASS").sum())
        flagged = int((status == "FLAGGED").sum())
        excluded = int((status == "EXCLUDED").sum())
        rows.append(dict(stage="3. Passed genome-level filters (PASS)", level="genome",
                         n=passed,
                         detail="length, %N<=1%, not engineered/synthetic, not duplicate"))
        rows.append(dict(stage="3a. Flagged but retained (FLAGGED)", level="genome",
                         n=flagged,
                         detail="kept in the analysis; reasons in S3_removed_sequences.tsv"))
        rows.append(dict(stage="3b. Excluded at genome level (EXCLUDED)", level="genome",
                         n=excluded, detail="see S3_removed_sequences.tsv"))

    # These file counts are what actually entered per-gene extraction, and they do not sum
    # to stage 3: they include the flagged-but-kept records, one reference anchor per type,
    # and the tier-2 set, which is drawn from the stage-3b excluded pool rather than from
    # the PASS records. 02_filtered_genomes/ is the authority for what was extracted.
    ready_total = 0
    for name, path, note in (
            ("HSV-1", "HSV1_filtered.fasta", "tier 1 + anchor"),
            ("HSV-2", "HSV2_filtered.fasta", "tier 1 + anchor"),
            ("HSV-2 tier 2", "HSV2_filtered_tier2.fasta",
             "readmitted from stage 3b (1-2% ambiguity); POOLED with tier 1 downstream, "
             "not analysed separately")):
        p = os.path.join(DIR_BASE, "02_filtered_genomes", path)
        if os.path.exists(p):
            with open(p, encoding="utf-8", errors="ignore") as fh:
                n = sum(1 for l in fh if l.startswith(">"))
            ready_total += n
            rows.append(dict(stage=f"4. Analysis-ready genomes ({name})", level="genome",
                             n=n, detail=f"{path} -- {note}"))
    if ready_total:
        rows.append(dict(stage="4a. Analysis-ready genomes (total)", level="genome",
                         n=ready_total,
                         detail="exceeds stage 3 because it includes FLAGGED records, "
                                "reference anchors and the readmitted tier-2 set"))

    qc = read_tsv(os.path.join(DIR_QC, "07_cds_qc_log.tsv"))
    if qc is not None:
        rows.append(dict(stage="5. CDS rows in the per-sequence QC log", level="CDS",
                         n=len(qc),
                         detail="07_cds_qc_log.tsv; does NOT cover the full carried-forward "
                                "set (see stage 6b)"))
        ok = int(qc.qc_status.astype(str).str.startswith("PASS").sum())
        rows.append(dict(stage="6. Of those, passing QC", level="CDS", n=ok,
                         detail="ATG start, stop codon, in frame, no internal stop, %N limit, "
                                "length ratio 0.95-1.05"))
        rows.append(dict(stage="6a. Of those, failing QC", level="CDS", n=len(qc) - ok,
                         detail="see S3_removed_sequences.tsv"))

    # The count that matters: what deduplication and every later stage actually read.
    carried = 0
    for f in glob.glob(os.path.join(DIR_BASE, "05_cds_qc_passed", "*", "*.fasta")):
        b = os.path.basename(f)
        if b.endswith("_dedup.fasta") or "resistance" in b:
            continue
        with open(f, encoding="utf-8", errors="ignore") as fh:
            carried += sum(1 for l in fh if l.startswith(">"))
    if carried:
        rows.append(dict(stage="6b. CDS carried forward to deduplication", level="CDS",
                         n=carried,
                         detail="05_cds_qc_passed/; exceeds stage 6 because a later "
                                "extraction top-up updated the FASTAs and 06_cds_extraction_log "
                                "but not 07_cds_qc_log -- this is the number to cite"))

    dedup = read_tsv(os.path.join(DIR_QC, "08_dedup_log.tsv"))
    if dedup is not None:
        rows.append(dict(stage="7. Unique sequences after exact-duplicate collapse",
                         level="CDS", n=int(dedup.unique_sequences.sum()),
                         detail=f"from {int(dedup.total_clean_sequences.sum())} clean CDS; "
                                "cluster sizes retained for population frequencies"))

    aln = read_tsv(os.path.join(DIR_V2, "alignment_qc_v2.tsv"))
    if aln is not None:
        per_type = aln[aln.cohort != "combined"]
        rows.append(dict(stage="8. Sequences in rebuilt per-type codon alignments",
                         level="CDS", n=int(per_type.n_sequences.sum()),
                         detail=f"{len(per_type)} alignments, includes one reference anchor each"))

    ready = glob.glob(os.path.join(DIR_READY, "*", "*_idmap.tsv"))
    if ready:
        total = sum(len(pd.read_csv(p, sep="\t")) for p in ready)
        rows.append(dict(stage="9. Sequences entering selection analysis", level="CDS",
                         n=total,
                         detail=f"{len(ready)} gene x type cohorts, reference anchor removed"))

    return pd.DataFrame(rows)


def removed_sequences():
    rows = []

    filt = read_tsv(os.path.join(DIR_QC, "04_genome_filtering_log.tsv"))
    if filt is not None:
        bad = filt[filt.filter_status.astype(str) != "PASS"]
        for r in bad.itertuples():
            rows.append(dict(level="genome", identifier=r.accession, gene="",
                             virus_type=getattr(r, "database_type", ""),
                             stage="genome filtering",
                             reason=str(getattr(r, "decision_reasons", "") or
                                        getattr(r, "filter_status", ""))))

    disagree = read_tsv(os.path.join(DIR_QC, "04e_species_disagreement_full_list.tsv"))
    if disagree is not None:
        for r in disagree.itertuples():
            rows.append(dict(level="genome", identifier=r.accession, gene="",
                             virus_type=str(getattr(r, "database_label", "")),
                             stage="type verification",
                             reason=f"database label {getattr(r, 'database_label', '')} vs "
                                    f"computed {getattr(r, 'computational_assignment', '')} "
                                    f"({getattr(r, 'classification', '')})"))

    qc = read_tsv(os.path.join(DIR_QC, "07_cds_qc_log.tsv"))
    if qc is not None:
        bad = qc[~qc.qc_status.astype(str).str.startswith("PASS")]
        for r in bad.itertuples():
            rows.append(dict(level="CDS", identifier=r.accession, gene=r.gene,
                             virus_type=r.virus, stage="per-gene CDS QC",
                             reason=str(getattr(r, "failure_reasons", "") or r.qc_status)))

    mask = read_tsv(os.path.join(DIR_V2, "masking_report_v2.tsv"))
    if mask is not None:
        for r in mask[mask.codons_dropped > 0].itertuples():
            rows.append(dict(level="alignment column", identifier=f"{r.gene} ({r.cohort})",
                             gene=r.gene, virus_type=r.cohort, stage="occupancy masking",
                             reason=f"{r.codons_dropped} codon columns below 50% clinical "
                                    f"occupancy; ranges {r.dropped_codon_ranges}"))

    out = read_tsv(os.path.join(DIR_V2, "outlier_sequences_v2.tsv"))
    if out is not None and len(out):
        per_type = out[out.cohort.isin(["HSV1", "HSV2"])]
        for r in per_type.itertuples():
            rows.append(dict(level="sequence (flagged, retained)",
                             identifier=r.sequence_id, gene=r.gene, virus_type=r.cohort,
                             stage="alignment outlier screen",
                             reason=f"consensus distance {r.consensus_distance} exceeds "
                                    f"cutoff {r.cutoff}; retained for investigation, "
                                    "not removed"))

    return pd.DataFrame(rows)


def final_accessions():
    rows = []
    for path in sorted(glob.glob(os.path.join(DIR_READY, "*", "*_idmap.tsv"))):
        cohort = os.path.basename(os.path.dirname(path))
        gene = os.path.basename(path).replace("_idmap.tsv", "")
        df = pd.read_csv(path, sep="\t")
        for r in df.itertuples():
            rows.append(dict(virus_type=cohort, gene=gene, safe_id=r.safe_id,
                             accession=r.accession,
                             isolates_represented=r.cluster_size))
    if not rows:
        return pd.DataFrame(), pd.DataFrame()
    acc = pd.DataFrame(rows)
    per_gene = (acc.groupby(["virus_type", "gene"])
                   .agg(unique_sequences=("safe_id", "nunique"),
                        isolates_represented=("isolates_represented", "sum"),
                        distinct_accessions=("accession", "nunique"))
                   .reset_index())
    return acc, per_gene


def software_versions():
    rows = []
    try:
        import importlib.metadata as im
        for pkg in ("biopython", "pandas", "numpy", "scipy", "statsmodels",
                    "scikit-learn", "pyfamsa", "mhcflurry", "matplotlib"):
            try:
                rows.append(dict(tool=pkg, version=im.version(pkg), source="pip"))
            except Exception:
                rows.append(dict(tool=pkg, version="not installed", source="pip"))
    except Exception:
        pass

    rows.append(dict(tool="python", version=sys.version.split()[0], source="interpreter"))

    iq = os.path.join(DIR_BASE, "bin", "iqtree", "iqtree-2.4.0-Windows", "bin", "iqtree2.exe")
    rows.append(dict(tool="IQ-TREE", version="2.4.0" if os.path.exists(iq) else "missing",
                     source="bundled in bin/"))

    prov = os.path.join(DIR_BASE, "14_selection", "gard_provenance.json")
    if os.path.exists(prov):
        with open(prov, encoding="utf-8") as fh:
            p = json.load(fh)
        rows.append(dict(tool="HyPhy", version=p.get("version", ""),
                         source=f"docker {p.get('image', '')} id {p.get('image_id', '')[:19]}"))
    else:
        rows.append(dict(tool="HyPhy", version="2.5.103 (MP)", source="docker hyphy/hyphy"))

    cov = os.path.join(DIR_BASE, "18_population_coverage",
                       "allele_frequency_provenance.json")
    if os.path.exists(cov):
        with open(cov, encoding="utf-8") as fh:
            c = json.load(fh)
        rows.append(dict(tool="IEDB allele frequencies", version=c.get("tool", ""),
                         source=f"sha256 {c.get('sha256', '')[:16]}"))

    # The B-cell predictors run in their own virtual environments, so their versions are
    # not discoverable from this interpreter. step28 writes them to
    # 13_epitopes/bcell_predictor_provenance.tsv at run time; read that rather than
    # hard-coding, and fall back to NOT RUN only if the step has not been run.
    prov_path = os.path.join(DIR_BASE, "13_epitopes", "bcell_predictor_provenance.tsv")
    bcell_done = set()
    if os.path.exists(prov_path):
        prov = pd.read_csv(prov_path, sep="\t")
        for r in prov.itertuples():
            rows.append(dict(tool=r.tool, version=r.version,
                             source=f"{r.source}; {r.mode}; licence: {r.licence}"))
            bcell_done.add(r.tool)

    # NetMHCIIpan has no local install: it was run on the DTU web server, and the server
    # records its exact command line in the first line of each output file. Read the
    # version from there rather than asserting one.
    import glob as _glob, gzip as _gzip, re as _re
    c2_raw = sorted(_glob.glob(os.path.join(DIR_BASE, "13_epitopes", "classii_raw",
                                            "*.xls*")))
    c2_ver = None
    for f in c2_raw:
        op = _gzip.open if f.endswith(".gz") else open
        with op(f, "rt", encoding="utf-8", errors="ignore") as fh:
            first = fh.readline()
        m = _re.search(r"(netMHCIIpan-[\d.]+)", first)
        if m:
            c2_ver = m.group(1)
            break
    if c2_ver:
        rows.append(dict(tool="NetMHCIIpan-4.3", version=c2_ver,
                         source="DTU Health Tech web server "
                                "(services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/); "
                                "15-mers, -BA, -rankS 1 -rankW 5; licence: academic "
                                "registration required for the standalone package"))
        bcell_done.add("NetMHCIIpan-4.3")

    for tool in ("NetMHCIIpan-4.3", "BepiPred-3.0", "DiscoTope-3.0"):
        if tool in bcell_done:
            continue
        note = ("requires DTU HealthTech academic licence"
                if tool.startswith("NetMHCIIpan")
                else "free for academic use; run step28_bcell_predictors.py")
        rows.append(dict(tool=tool, version="NOT RUN", source=note))
    return pd.DataFrame(rows)


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 24: building publication supplements")

    flow = data_flow()
    flow.to_csv(os.path.join(DIR_OUT, "S1_data_flow.tsv"), sep="\t", index=False)
    log(f"S1_data_flow.tsv ({len(flow)} stages)")
    for r in flow.to_dict("records"):
        log(f"    {r['stage']:<52} {r['n']:>7}  [{r['level']}]")

    acc, per_gene = final_accessions()
    if len(acc):
        acc.to_csv(os.path.join(DIR_OUT, "S2_final_accessions.tsv"), sep="\t", index=False)
        per_gene.to_csv(os.path.join(DIR_OUT, "S4_per_gene_cohorts.tsv"),
                        sep="\t", index=False)
        log(f"S2_final_accessions.tsv ({len(acc)} rows, "
            f"{acc.accession.nunique()} distinct accessions, "
            f"{int(acc.isolates_represented.sum())} isolates represented)")
        log(f"S4_per_gene_cohorts.tsv ({len(per_gene)} gene x type cohorts)")
    else:
        log("S2/S4 skipped: no HyPhy-ready id maps found yet")

    removed = removed_sequences()
    if len(removed):
        removed.to_csv(os.path.join(DIR_OUT, "S3_removed_sequences.tsv"),
                       sep="\t", index=False)
        log(f"S3_removed_sequences.tsv ({len(removed)} rows)")
        for stage, n in removed.stage.value_counts().items():
            log(f"    {stage}: {n}")
    else:
        log("S3 skipped: no exclusion logs found")

    sw = software_versions()
    sw.to_csv(os.path.join(DIR_OUT, "S5_software_versions.tsv"), sep="\t", index=False)
    log(f"S5_software_versions.tsv ({len(sw)} tools)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
