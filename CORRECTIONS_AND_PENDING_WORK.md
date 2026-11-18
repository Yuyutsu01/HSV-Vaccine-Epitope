# HSV-1 / HSV-2 Conservation, Selection and Vaccine-Epitope Mapping
## Audit of existing work, required corrections, and pending work

**Audit date:** 2026-11-19
**Audited against:** the project proposal (pipeline steps 1–11) and the fine-curation protocol (steps 0–11)
**Repository:** `C:\Users\shiva\OneDrive\Desktop\vaccine project`

---

## 0. Executive summary

| Layer | Proposal steps | Verdict |
|---|---|---|
| Sequence acquisition & genome filtering | Curation 1–2 | **KEEP** — meets protocol |
| Type verification (ANI, not DB labels) | Curation 2 | **KEEP** — meets protocol |
| Metadata curation | Curation 3 | **PARTIAL** — built, but source data too sparse; downstream claims must be dropped |
| CDS extraction | Curation 4 | **KEEP** — 92–99% success, validated against reference translations |
| Per-gene CDS QC | Curation 5 | **KEEP** — ATG/stop/frame/internal-stop/%N all enforced |
| Dedup | Curation 8 | **KEEP** — exact-duplicate collapse with cluster sizes retained |
| **Codon-aware alignment** | Pipeline 2 / Curation 6 | **REBUILD FROM SCRATCH** — not an MSA; out-of-register for all indel-bearing genes |
| **Repeat masking** | Curation 5 | **RE-DERIVE** — currently masks the symptom of the alignment bug |
| **Phylogenies** | Pipeline 4 | **RE-RUN** — built on corrupted alignments |
| **Recombination (GARD)** | Pipeline 3 / Curation 7 | **RE-RUN** — ran last instead of first, and on pooled HSV-1+HSV-2 data |
| **Selection (BUSTED/FEL/MEME/FUBAR/aBSREL)** | Pipeline 5 | **RE-RUN** — pooled cohorts, unpartitioned, no FDR |
| Conservation scoring | Pipeline 6 | **PENDING** (a diagnostic table exists; the per-residue layer does not) |
| Epitope prediction (B-cell, T-cell) | Pipeline 7 | **PENDING — not started** |
| IEDB curation & validation | Pipeline 9 / Curation 10 | **PENDING — not started** |
| Integration & ranking | Pipeline 8 | **PENDING — not started** |
| Structure mapping | Pipeline 10 | **PENDING — not started** |
| ESM zero-shot (optional) | ML section | **PENDING — not started** |
| Reproducibility (workflow, versions, supplements) | Curation 11 | **PENDING** |

**Headline:** roughly **40% of the completed work is sound and reusable** (everything up to and including deduplication), **35% must be re-run** (alignment, trees, recombination, selection), and **the remaining ~25% of the paper — the entire epitope half — has not been started.** The README's "Status: 100% Completed" badge is wrong on both counts.

---

## PART A — EXISTING WORK THAT NEEDS CORRECTION

### A1. CRITICAL — The alignments are not multiple sequence alignments

**Where:** `09_scripts/step08_dedup_and_step09_align.py`, function `codon_align_records()`; outputs in `06_alignments/{HSV1,HSV2,combined}/`.

**What it does:** each sequence is translated, aligned **pairwise against the reference protein** with `Bio.Align.PairwiseAligner` (match +2 / mismatch −1 / gap-open −4 / gap-extend −1), and its codons are projected into that pairwise alignment's coordinate frame. The variable `all_col_insertions` is declared and never used — there is no profile, no progressive alignment, and no reconciliation of indel placement between sequences. The proposal specified **MAFFT or MACSE**. Neither is installed in this repository (`bin/` contains only `datasets`, `dataformat` and IQ-TREE; `bin/hyphy/` is empty).

**Consequence:** every sequence that carries an indel relative to the reference has its gaps placed independently, so **all columns downstream of that indel are out of register with the rest of the alignment.** Direct evidence, `06_alignments/HSV1/US8_aligned.fasta`, two HSV-1 gE rows at nt 550:

```
PZ347688  GCGAGGGCGAGGACGAAAGTCTAGCCGGCACTCCCGCCAGCGGGACCCCCCGGCTCCCGCCTCCCCCCGC
ON783216  GCAGGGACGAAGATCTCGCCGGCACTCCCGCCAGCGGGGACCCG---GCTCCGCCTCCCCCCGCC---CC
```

Clearly homologous sequence, shifted by a few nucleotides. Measured impact:

| Gene | within-HSV-1 nt identity (as aligned) | columns with major allele <90% | real within-type identity |
|---|---|---|---|
| UL27 (gB) — no indels | 99.3% | 1.3% | correct |
| UL30, US6, UL22, UL1, UL10, UL53 | 97–99% | low | correct |
| **US8 (gE)** | **81.7%** | **46.3%** | should be ~98–99% |
| **US7 (gI)** | **81.2%** | **39.9%** | should be ~97–99% |
| **US4 (gG)** | **79.9%** | 15.3%, **66% of columns are gap padding** | should be ~98% |
| **RS1 (ICP4)** | **67.7%** (HSV-2: 51.8%) | — | ≥99% (confirmed by masked run: 97.8–99.1%) |
| **RL2 (ICP0)** | 64.2% | — | ≥99% (masked: 99.3%) |

Additional defect in the same function: `06_alignments/HSV1/US4_aligned.fasta` is **699 codons wide while every sequence in it is 238 aa** — two thirds of that alignment is empty gap padding, because the HSV-1 gG alignment inherited the HSV-2 gG reference length. HyPhy was run on this file.

**Verdict:** genes with no clinical indels (UL27, UL30, US6, UL22, UL1, UL10, UL53, UL48, UL54, UL23, UL49A) degenerate to a correct ungapped stack and their alignments are salvageable. Genes with indels (**US4, US7, US8, UL49, RL2, RS1, UL44, UL46, UL47**) are corrupted. Because the fix is the same script for all of them, **rebuild all 40 per-type alignments with a real codon-aware MSA.**

**Correction required:** replace step 9 with MAFFT (`--auto`, protein-level, codon back-projection) or MACSE; re-run for all 20 genes × {HSV-1, HSV-2}. Everything downstream must then be re-run.

---

### A2. CRITICAL — Pooled HSV-1 + HSV-2 "combined" cohorts were used for dN/dS and recombination

**Where:** `06_alignments/combined/*`, built by pairwise-aligning the two **reference proteins** and projecting both species onto that map (`step08...align.py`, "Reference-Anchored Merger"). Used by `step11_stage1_pipeline.py`, `step12_stage2_pipeline.py`, `step12_gard_pipeline.py`.

**Why this is invalid:** `07_qc_reports/09_conservation_metric_diagnostic.tsv` records between-type amino-acid identity of **17.9% (US4/gG) to 86.1% (UL27/gB)**. A single codon model cannot be fitted across that range, and:

- **GARD on a two-type alignment detects the HSV-1/HSV-2 topological split, not recombination.** All four "RECOMBINANT" verdicts in the README come from combined runs: UL27 nt 735 (ΔAICc 310.8), UL22 nt 604 (468.2), UL30 nt 2018 (**6,619.1**), UL46 nt 1713 (2,385.7). A ΔAICc of 6,619 is orders of magnitude beyond what intra-species recombination in a herpesvirus produces — it is the species split.
- **MEME on combined UL27** returns 24 sites at p<0.05 concentrated at codons 2, 3, 13, 14, 15, 19, 21, 23, 26, 27, 33, 36, 40, 44, 49, 58 — precisely the N-terminal region where gB differs between types. Type divergence is being scored as episodic positive selection.
- `US4_combined` was analysed at 17.88% between-type identity.

**Correction required:**
1. Delete all `*_combined` results from Stage 1, Stage 2 and GARD (22 Stage 1 units, 12 Stage 2 units, 10/10 GARD runs).
2. Re-run GARD and all selection analyses **per type only**.
3. Rebuild combined alignments with a true MSA and use them for **one purpose only**: between-type conservation comparison and cross-protective epitope identification (proposal pipeline step 6). Never for dN/dS or GARD.
4. Retract the four domain-architecture interpretations in the README until per-type GARD results exist.

---

### A3. CRITICAL — Recombination screening ran *after* selection, and breakpoints were never applied

**Evidence (progress logs):** Stage 1 selection 2026-10-28 → 11-06; Stage 2 selection 11-11 → 11-16; **GARD 11-11 → 11-18**.

The proposal (pipeline step 3, curation step 7) requires recombination screening **first**, then splitting alignments at breakpoints and analysing segments separately. `step12_gard_pipeline.py` has no partition-export logic, so no alignment was ever split. Every selection result in the repository was computed on an unpartitioned alignment.

**Correction required:** reorder the pipeline to MSA → QC → trees → GARD (per type) → partition → selection. Export GARD partitions as separate FASTA segments with their own trees.

---

### A4. CRITICAL — No multiple-testing correction anywhere; p<0.10 is being reported

A repository-wide search for `fdr`, `benjamini`, `multipletests`, `holm` across every `.py`, `.md` and `.tsv` returns **zero hits**. `stage2_summary.tsv` publishes `fel_sites_p10` and `meme_sites_p10` columns; `stage1_summary.tsv` triggers Stage 2 on raw `busted_p` and uncorrected FUBAR posteriors.

Scale of the problem: 198 Stage 1 runs + 62 Stage 2 runs, each testing 91–1,318 codons. `UL27_combined` alone reports 24 MEME sites at "p<0.05".

**Correction required:**
- Benjamini–Hochberg FDR within each alignment for FEL and MEME; report q-values, headline q<0.05.
- FUBAR: report posterior probability ≥0.95 only (drop the ≥0.90 column from headline results).
- aBSREL: use the Holm–Bonferroni-corrected branch p-values that HyPhy already emits.
- BUSTED: correct across the gene × cohort family.
- Require **concordance across ≥2 methods** before a site enters the epitope shortlist, as the proposal specified.
- Install `statsmodels` (currently missing) or implement BH directly.

---

### A5. MAJOR — Repeat masking currently hides the alignment bug rather than fixing it

`11b_tier_b_masking_investigation.tsv` correctly identifies low-complexity tracts in RL2/RS1 and reports `CONFIRMED_REPEAT_ARTIFACT`, and masking does raise within-type identity to 97.8–99.4%. But it achieves this by **deleting the misaligned columns** (RL2 HSV-1 2568→1935 nt, RS1 4167→3375 nt, US4 2097→654 nt — 69% of gG discarded). For US8 the masked file is **byte-identical to the unmasked file** — nothing was masked, and gE remains corrupted at 81.7% identity while its Stage 1 tier-C results are reported.

**Correction required:** re-derive masking **after** a correct MSA. Mask only genuinely unalignable low-complexity tracts, document coordinates in reference numbering, and re-check whether RL2/RS1/US4 still need masking at all once aligned properly. Expect far less masking to be necessary.

---

### A6. MAJOR — Phylogenies must be rebuilt

All 60 IQ-TREE runs in `06_alignments/trees/` (HSV1 / HSV2 / combined) were built on the corrupted alignments, and the combined trees additionally span two species. HyPhy consumed these trees for every FEL/MEME/FUBAR/aBSREL run. IQ-TREE 2.4.0 is present locally at `bin/iqtree/iqtree-2.4.0-Windows/bin/iqtree2.exe`, so this is cheap to redo.

**Correction required:** rebuild per-type, per-gene (and per-GARD-partition) trees with model selection and 1,000 ultrafast bootstraps. Keep the IQ-TREE name-mangling map (`parse_iqtree_mapping` already handles this correctly — reuse it).

---

### A7. MAJOR — Metadata is too sparse to support several planned analyses

`07_qc_reports/05_metadata_completeness.tsv` (n = 1,305):

| Field | Completeness |
|---|---|
| host | 100% |
| isolate | 97.1% |
| country | **17.6%** |
| collection year | **8.4%** |
| BioProject | 3.3% |
| BioSample | 2.8% |
| isolation source | **0.08%** (1 record) |
| passage history | **0%** |
| sequencing technology | **0%** |

**Consequences:** curation step 9 (cap per country-year, balanced subsample sensitivity check) cannot be run meaningfully; the anatomical-source analysis cannot be run; **the CNS/encephalitis comparison is impossible** (`is_cns_engephalitis` is False for essentially every record). Deduplication is identity-based, so patient-level independence was never enforced — 145 records are serial samples from ~40 patients (`Z_d40_skin`, `C_d10_plasma`, `HSV1-v54_*`).

**Correction required:**
1. Back-fill country/year from the source publications and BioSample records where recoverable (one pass over the ~43 BioProjects).
2. Build an explicit `patient_id` column from isolate-name parsing plus BioProject, and enforce **one sequence per patient per gene** in a sensitivity cohort; report main results on both.
3. Delete the encephalitis/CNS analysis from the project scope, or source a dedicated CNS isolate set, and remove the claim from the proposal.
4. State the sampling-bias limitation explicitly with these numbers.

---

### A8. MODERATE — Low power for several HSV-2 genes; report it rather than interpret it

Post-dedup unique sequence counts (`08_dedup_log.tsv`): HSV-2 **UL1 = 32**, US6 = 54, UL53 = 74, RS1 = 57; HSV-1 UL49A = 32. The proposal warned that <100 per type is a power problem. "STRICTLY CLONAL" / "high structural conservation" conclusions for gD and gL in HSV-2 are being drawn from cohorts where no test could have detected anything.

**Correction required:** report per-cohort n alongside every result; replace "strictly clonal" with "no breakpoint detected at n = X"; add a power statement.

---

### A9. MODERATE — README contains countable factual errors and unsupported conclusions

| README claim | Reality |
|---|---|
| "342 / 342 Runs Complete"; "Stage 1 = 270/270" | Stage 1 = **198** runs (the log itself says "66 alignments x 3 methods = 198 runs"). 198 + 62 + 10 = **270** total. The 342 double-counts Stage 1. |
| "33 targets" | **20 genes** (`stage1_summary.tsv`) |
| "Tier B (Expanded Cohorts across 33 targets × 3 methods): 99/99" | tier_b = **RL2 and RS1 only**, masked vs unmasked, **36** runs |
| "Status: 100% Completed" | the epitope half of the project does not exist |
| gD "confirmed as an optimal monophyletic subunit candidate"; gL "stable co-immunogen" | not supportable from a dN/dS screen; requires the missing epitope and validation work |
| "STRICTLY CLONAL" from ΔAICc = 0.00 | absence of a detected breakpoint ≠ clonality |
| Directory tree | omits `10_hyphy_outputs/`, `bin/`, `logs/`; lists 2 of 28 scripts |
| "HyPhy v2.5.103 (MP)" | `bin/hyphy/` is empty; HyPhy runs via `docker run hyphy/hyphy` with no pinned tag, so the version is unverifiable |

**Correction required:** rewrite the README after the re-runs, with corrected counts, per-type results only, q-values, n per cohort, and no conclusions beyond the evidence.

---

### A10. MODERATE — Reproducibility gaps

- Git: a single `Initial commit`; the entire dataset, all 28 scripts and all results are **uncommitted**.
- No Snakemake/Nextflow, no pinned environment, no `requirements.txt` with versions.
- `bin/*.exe` is gitignored, so tool provenance is lost.
- HyPhy is only reachable through Docker, and the Docker daemon is not currently running.
- `06_alignments/HSV1/test_UL27.*` and `hyphy/test_UL27_HSV1_fubar*.json` are scratch files sitting in the results tree.

**Correction required:** commit a snapshot now; add `requirements.txt` and a `workflow/Snakefile`; pin the HyPhy image by digest; move scratch files out of the results tree; publish the final accession list and the removed-sequence table as supplements.

---

## PART B — PENDING WORK (NOT STARTED)

Verified absent: a repository-wide search for `iedb`, `epitope`, `bepipred`, `netmhc`, `discotope`, `alphafold`, `esm` returns **zero hits** in any script, document or table.

### B1. Conservation scoring layer (proposal pipeline step 6)
Per-residue Shannon entropy and percent identity, computed **within HSV-1, within HSV-2, and between types**, in reference residue numbering, with cluster-size weights from `08_dedup_log.tsv` so population frequencies are preserved. `09_conservation_metric_diagnostic.tsv` has gene-level means only — the per-residue layer that the epitope ranking depends on does not exist.

### B2. IEDB epitope curation (pipeline step 9, curation step 10)
Download HSV-1 and HSV-2 entries; keep human host and positive assays; record assay type and HLA allele; map each peptide to the reference protein by exact match and drop those that do not map cleanly; **do not treat untested regions as negatives.**

### B3. B-cell epitope prediction (pipeline step 7)
- Linear: BepiPred-3.0.
- Conformational: DiscoTope-3 on PDB structures where available (gD, gB, gH/gL) and AlphaFold models otherwise.

### B4. T-cell epitope prediction (pipeline step 7)
- CD8: NetMHCpan-4.1 or MHCflurry across a defined panel of common HLA class I alleles.
- CD4: NetMHCIIpan.
- Record the allele panel and thresholds up front.

### B5. Predictor benchmarking against IEDB (ML section)
AUC and precision-at-top-k for each predictor's recovery of known epitopes, with **gene-level grouping** so the evaluation is honest.

### B6. Integration and ranking (pipeline step 8)
Per candidate epitope: predicted binding/antigenicity, conservation (within and between type), overlap with positively selected sites (q<0.05), overlap with masked or low-complexity regions, predicted N-glycosylation shielding. Keep epitopes that are predicted, highly conserved and **not** positively selected. Highlight those shared by both types. Add HLA population coverage via the IEDB tool.

### B7. Validation against known biology (pipeline step 9)
Check that the pipeline recovers published neutralising antibody epitopes on gD and gB and known immunodominant T-cell epitopes from IEDB. This is the substitute for wet-lab validation and is what makes the ranking credible.

### B8. Structure mapping (pipeline step 10)
Map conserved epitopes and selected sites onto 3D structures (PyMOL or ChimeraX); produce the figures.

### B9. ESM zero-shot analysis (optional)
ESM-2 per-position substitution scores; test correlation against observed isolate diversity and per-site dN/dS. `torch` is already installed. Nothing is trained, so there is no leakage.

### B10. Writing and supplements
Methods section (`08_documentation/methods_draft.md` is currently 2 paragraphs), data-flow table with counts surviving each filter, removed-sequence table with reasons, accession list.

---

## PART C — EXECUTION ORDER

Dependencies are strict: nothing in Phase 4+ is meaningful until Phase 1 is done.

| Phase | Work | Tool gate | Status |
|---|---|---|---|
| **0** | Git snapshot; `statsmodels`; obtain MAFFT; start Docker for HyPhy | Docker daemon **down**; MAFFT **not installed** | |
| **1** | Rebuild all 40 per-type codon MSAs + 20 combined (conservation use only) | MAFFT | |
| **2** | Alignment QC (per-column identity, outlier/long-branch check); re-derive masking | — | |
| **3** | Rebuild per-type trees | IQ-TREE 2.4.0 ✓ local | |
| **4** | GARD per type; export breakpoint partitions + partition trees | HyPhy (Docker) | |
| **5** | BUSTED / FEL / MEME / FUBAR / aBSREL per type per partition; BH-FDR pass; concordance table | HyPhy (Docker) | |
| **6** | Per-residue conservation layer (within HSV-1, within HSV-2, between types) | — | |
| **7** | IEDB download and curation | network | |
| **8** | B-cell and T-cell epitope prediction; predictor benchmarking | web servers / MHCflurry | |
| **9** | Integration, ranking, population coverage, known-epitope validation | — | |
| **10** | Structure mapping; ESM zero-shot; README and methods rewrite; Snakemake; supplements | — | |

**Blocking prerequisites for the user:**
1. **Start Docker Desktop** — HyPhy is only reachable via `docker run hyphy/hyphy` and the daemon is not running. Phases 4–5 cannot execute without it.
2. **NetMHCpan-4.1 / NetMHCIIpan / BepiPred-3.0 / DiscoTope-3** require academic registration or are web-server only. Decide per tool: local licensed install, web submission, or substitute (MHCflurry is pip-installable and needs no licence).

---

## PART D — WHAT TO SALVAGE, FILE BY FILE

**Keep as-is:**
`00_raw_downloads/`, `01_reference_genomes/`, `02_filtered_genomes/`, `03_metadata/`, `04_cds_unaligned/`, `05_cds_qc_passed/`, and all of `07_qc_reports/04*`, `05*`, `06*`, `07*`, `08*` (genome filtering, type verification, extraction, CDS QC, dedup logs). Scripts `step02`–`step07` are sound.

**Rebuild:**
`06_alignments/**` (all alignments, masked sets and trees); `09_scripts/step08_dedup_and_step09_align.py` step-9 half.

**Delete / supersede after re-run:**
`07_qc_reports/hyphy/stage1/**`, `stage2/**`, `gard/**` and the three summary TSVs; `07_qc_reports/09_conservation_metric_diagnostic.tsv` (superseded by the per-residue layer); `11b`/`11c` masking tables (re-derive); all `test_UL27*` scratch files.

**Rewrite:**
`README.md`; `08_documentation/methods_draft.md`.
