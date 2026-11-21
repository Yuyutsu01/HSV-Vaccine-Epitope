"""
Step 20 (CORRECTION A7): Patient-level identity and sampling-bias tabulation.

Deduplication in step08 collapsed identical sequences, which is the right unit for HyPhy but
not for independence: 145 records are serial samples from the same patients (Z_d10_skin,
Z_d39_plasma, Z_d40_skin, HSV1-v54_*), so one patient can contribute several haplotypes to a
cohort (CORRECTIONS_AND_PENDING_WORK.md, A7).

This step:
  1. Derives a patient_id from the isolate naming schemes, BioSample, or the accession.
  2. Tabulates sampling bias by country, year, submitter and BioProject, so curation step 9
     can be reported honestly instead of asserted.
  3. Emits, per gene and type, a one-sequence-per-patient subcohort for use as a sensitivity
     analysis alongside the main result.
  4. States what the metadata cannot support: with isolation_source at 0.08% completeness and
     exactly one record flagged CNS/encephalitis, the planned anatomical-source and
     encephalitis comparisons are not feasible from this dataset.

Outputs:
  - 03_metadata/patient_assignments.tsv
  - 03_metadata/sampling_bias_tables.tsv
  - 03_metadata/one_per_patient_cohorts.tsv
  - 03_metadata/metadata_feasibility_statement.md
  - 03_metadata/patient_bias.log
"""

import os
import re
import sys
import glob
from datetime import datetime

import pandas as pd

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_META = os.path.join(DIR_BASE, "03_metadata")
DIR_READY = os.path.join(DIR_BASE, "06_alignments_v2", "hyphy_ready")

LOG_PATH = os.path.join(DIR_META, "patient_bias.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def derive_patient(row):
    """
    Patient identity, most reliable source first.
      - serial-sample isolate names: <patient>_d<day>_<site>[_culture]  -> <patient>
      - BioSample, when present, is one biological sample
      - otherwise the accession, i.e. treated as its own patient (conservative)
    Returns (patient_id, basis, sample_day).
    """
    iso = str(row.isolate) if pd.notna(row.isolate) else ""
    m = re.match(r"^(.*?)_d(\d+)_", iso)
    if m:
        return m.group(1), "serial_isolate_name", int(m.group(2))
    if pd.notna(row.BioSample) and str(row.BioSample).strip():
        return str(row.BioSample).strip(), "biosample", None
    if iso and re.match(r"^[A-Za-z]{1,4}\d{1,5}$", iso):
        return iso, "isolate_code", None
    return str(row.accession), "accession_singleton", None


def main():
    os.makedirs(DIR_META, exist_ok=True)
    log("=" * 78)
    log("STEP 20: patient-level identity and sampling-bias tabulation")

    meta = pd.read_csv(os.path.join(DIR_META, "curated_metadata.tsv"), sep="\t")
    log(f"Loaded {len(meta)} metadata rows")

    derived = [derive_patient(r) for r in meta.itertuples()]
    meta["patient_id"] = [d[0] for d in derived]
    meta["patient_basis"] = [d[1] for d in derived]
    meta["sample_day"] = [d[2] for d in derived]

    multi = meta.groupby("patient_id").accession.count()
    multi = multi[multi > 1].sort_values(ascending=False)
    log(f"{meta.patient_id.nunique()} distinct patients across {len(meta)} records")
    log(f"{len(multi)} patients contribute more than one record "
        f"({int(multi.sum())} records in total)")
    for pid, n in multi.head(12).items():
        log(f"    patient {pid}: {n} records")
    log("Basis for patient assignment: "
        + ", ".join(f"{k}={v}" for k, v in meta.patient_basis.value_counts().items()))

    meta[["accession", "virus", "isolate", "patient_id", "patient_basis", "sample_day",
          "country", "collection_year", "BioSample", "BioProject", "submitter"]].to_csv(
        os.path.join(DIR_META, "patient_assignments.tsv"), sep="\t", index=False)

    # Sampling bias tables.
    rows = []
    for field in ("country", "collection_year", "submitter", "BioProject", "virus"):
        vc = meta[field].fillna("UNKNOWN").astype(str).value_counts()
        for value, n in vc.items():
            rows.append(dict(field=field, value=value, records=int(n),
                             pct_of_records=round(100 * n / len(meta), 2),
                             patients=int(meta.loc[meta[field].fillna("UNKNOWN").astype(str)
                                                   == value, "patient_id"].nunique())))
    bias = pd.DataFrame(rows)
    bias.to_csv(os.path.join(DIR_META, "sampling_bias_tables.tsv"), sep="\t", index=False)
    log(f"Wrote sampling_bias_tables.tsv ({len(bias)} rows)")

    for field in ("country", "collection_year"):
        sub = bias[bias.field == field].sort_values("records", ascending=False)
        unknown = sub[sub.value == "UNKNOWN"].records.sum()
        log(f"  {field}: {round(100 * unknown / len(meta), 1)}% UNKNOWN; "
            f"top known: "
            + ", ".join(f"{r.value} ({r.records})"
                        for r in sub[sub.value != "UNKNOWN"].head(5).itertuples()))

    # One-sequence-per-patient subcohorts, per gene and type.
    acc2pat = dict(zip(meta.accession.astype(str), meta.patient_id.astype(str)))
    acc2day = dict(zip(meta.accession.astype(str), meta.sample_day))

    cohort_rows = []
    for cohort in ("HSV1", "HSV2"):
        for idmap in sorted(glob.glob(os.path.join(DIR_READY, cohort, "*_idmap.tsv"))):
            gene = os.path.basename(idmap).replace("_idmap.tsv", "")
            ids = pd.read_csv(idmap, sep="\t")
            ids["patient_id"] = [acc2pat.get(str(a), str(a)) for a in ids.accession]
            ids["sample_day"] = [acc2day.get(str(a)) for a in ids.accession]

            # Keep the earliest sampled record per patient; ties go to the larger cluster.
            ids["_day"] = ids.sample_day.fillna(10 ** 6)
            keep = (ids.sort_values(["_day", "cluster_size"], ascending=[True, False])
                       .drop_duplicates("patient_id"))
            for r in ids.itertuples():
                cohort_rows.append(dict(
                    cohort=cohort, gene=gene, safe_id=r.safe_id, accession=r.accession,
                    patient_id=r.patient_id, cluster_size=r.cluster_size,
                    in_one_per_patient_cohort=r.safe_id in set(keep.safe_id)))
            log(f"  {cohort} {gene}: {len(ids)} sequences from {ids.patient_id.nunique()} "
                f"patients -> {len(keep)} in the one-per-patient subcohort")

    cohorts = pd.DataFrame(cohort_rows)
    cohorts.to_csv(os.path.join(DIR_META, "one_per_patient_cohorts.tsv"),
                   sep="\t", index=False)
    log(f"Wrote one_per_patient_cohorts.tsv ({len(cohorts)} rows)")

    # Feasibility statement, for the limitations section.
    comp = {c: round(100 * meta[c].notna().mean(), 2)
            for c in ("country", "collection_year", "isolation_source", "passage_history",
                      "sequencing_technology", "BioSample", "BioProject")}
    cns = int((meta.is_cns_encephalitis == True).sum())
    dropped = int((cohorts.groupby(["cohort", "gene"]).in_one_per_patient_cohort
                   .apply(lambda s: (~s).sum())).sum())

    text = f"""# Metadata feasibility and sampling bias

Generated by `step20_patient_and_sampling_bias.py` on {datetime.now():%Y-%m-%d}.

## Completeness of the curated metadata (n = {len(meta)} genome records)

| Field | Complete |
|---|---|
| country | {comp['country']}% |
| collection year | {comp['collection_year']}% |
| isolation source | {comp['isolation_source']}% |
| passage history | {comp['passage_history']}% |
| sequencing technology | {comp['sequencing_technology']}% |
| BioSample | {comp['BioSample']}% |
| BioProject | {comp['BioProject']}% |

## What this supports

- **Patient-level independence.** {meta.patient_id.nunique()} distinct patients are resolvable
  across {len(meta)} records; {len(multi)} patients contribute more than one record. A
  one-sequence-per-patient subcohort is provided per gene and type
  (`one_per_patient_cohorts.tsv`); it removes {dropped} sequence entries in total across all
  gene x type cohorts. Report the main analysis on the full deduplicated cohort and this
  subcohort as a sensitivity analysis.
- **Submitter and study concentration.** The largest single submitter contributes
  {int(bias[bias.field == 'submitter'].records.max())} records. Study concentration, not
  geography, is the dominant bias in this dataset and should be stated as such.

## What this does NOT support, and must be removed from the project scope

- **Geographic stratification.** Country is known for only {comp['country']}% of records, so
  capping per country-year group or running a balanced geographic subsample would be based on a
  small, non-random fraction. Report the bias; do not model it.
- **Temporal analysis.** Collection year is known for {comp['collection_year']}% of records.
- **Anatomical-source comparison.** isolation_source is populated for
  {comp['isolation_source']}% of records.
- **The CNS / encephalitis comparison.** Exactly {cns} record is flagged
  `is_cns_encephalitis`. This analysis is not feasible from this dataset and should be dropped
  from the proposal, or served by a separately sourced CNS isolate set.
- **Passage-history sensitivity analysis.** passage_history is {comp['passage_history']}%
  complete, so the planned "analyse with and without passaged strains" check can only use the
  heuristic `lab_strain_flag`, not real passage records.

## Recommended recovery

Back-filling country and year from the {int(meta.BioProject.notna().sum())} records that carry a
BioProject, plus the source publications, is the only route to a usable geographic or temporal
analysis. Until then these analyses stay out of the paper.
"""
    with open(os.path.join(DIR_META, "metadata_feasibility_statement.md"), "w",
              encoding="utf-8") as fh:
        fh.write(text)
    log("Wrote metadata_feasibility_statement.md")
    log(f"Feasibility: country {comp['country']}%, year {comp['collection_year']}%, "
        f"isolation_source {comp['isolation_source']}%, CNS-flagged records {cns}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
