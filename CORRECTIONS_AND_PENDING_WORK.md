# HSV-1 / HSV-2 Conservation, Selection and Vaccine-Epitope Mapping
## Audit of existing work, required corrections, and pending work

**Audit date:** 2026-11-19
**Audited against:** the project proposal (pipeline steps 1–11) and the fine-curation protocol (steps 0–11)
**Repository:** `C:\Users\shiva\OneDrive\Desktop\vaccine project`

---

## 0. Executive summary

> **This section is the audit snapshot of 2026-11-19 and is kept unchanged as the record of
> what was found.** Every "RE-RUN", "REBUILD" and "PENDING" verdict below has since been
> acted on — see the execution log in Part B and "What is left" at the end of this file for
> current status. As of 2026-11-21 the class II and B-cell arms are complete; the former
> only outstanding items were the predictors once recorded as licence-gated
> predictors (user action) and ESM-2 zero-shot (dropped by decision).

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

| Phase | Work | Tool gate | Status | Script |
|---|---|---|---|---|
| **0** | `statsmodels`, `pyfamsa`, `mhcflurry` installed; Docker Desktop started | — | **DONE** | — |
| **1** | Rebuild all 40 per-type codon MSAs + 20 combined (conservation use only) | FAMSA | **DONE** | `step09b_rebuild_codon_msa.py` |
| **2** | Alignment QC, occupancy masking, reference coordinate maps, outlier flags | — | **DONE** | `step09c_mask_and_coordinate_map.py` |
| **3** | Rebuild per-type trees | IQ-TREE 2.4.0 local | **DONE** (40/40, 4 needed `-safe`) | `step10b_prep_and_trees.py` |
| **4** | GARD per type; export breakpoint partitions + partition trees | HyPhy 2.5.103 | **DONE** (40/40; 16 with breakpoints, 37 partitions) | `step16_gard_per_type.py` |
| **5** | BUSTED / FEL / MEME / FUBAR / aBSREL per type per partition | HyPhy 2.5.103 | **DONE** (305/305 runs) | `step17_selection_per_type.py` |
| **5b** | BH-FDR, direction, cross-method concordance, reference numbering | — | **DONE** (no `--allow-incomplete`; provisional marker removed) | `step18_fdr_and_concordance.py` |
| **6** | Per-residue conservation layer (within HSV-1, within HSV-2, between types) | — | **DONE** | `step13_conservation_layer.py` |
| **7** | IEDB download and curation | network | **DONE** | `step14_iedb_curation.py` |
| **8a** | CD8 class I prediction + conservation join | MHCflurry 2.3 | **DONE** | `step15_tcell_class1_prediction.py` |
| **8b** | Predictor benchmarking against IEDB (B5 + B7) | — | **DONE** | `step19_benchmark_predictors.py` |
| **8c** | CD4 class II, B-cell linear, B-cell conformational | was recorded as licence-gated | **DONE 2026-11-21** | BepiPred-3.0 and DiscoTope-3.0 are not licence-gated at all and were run locally; NetMHCIIpan-4.3 is, and was run on the DTU web server |
| **9a** | Integration and ranking | — | **SUPERSEDED** by 9c | `step21_integrate_and_rank.py` |
| **9b** | HLA class I population coverage, 239 populations | — | **DONE** | `step23_population_coverage.py` |
| **9c** | Re-run 9a with the selection filter applied | — | **DONE** (`selection_filter_applied=True`, 9 peptides excluded) | `step21_integrate_and_rank.py` |
| **10a** | Structure mapping (2C36 gD, 2GUM gB, 3M1C gH-gL) | — | **DONE**, re-run after 5b | `step22_structure_mapping.py` |
| **10b** | Pinned environment + Snakemake workflow | — | **DONE** | `requirements.txt`, `Snakefile` |
| **10c** | Supplements S1-S5 (data flow, accessions, removals, cohorts, versions) | — | **DONE**, rebuilt after 5b | `step24_build_supplements.py` |
| **10d** | README rewrite; methods draft rewrite | — | **DONE** | — |
| **10e** | ESM-2 zero-shot | — | **DROPPED** by decision | — |

### Results of the corrections completed so far

**Phase 1 — the alignment rebuild is the single biggest fix.** Within-type nucleotide identity,
before and after, with masked lengths now landing exactly on the true reference protein lengths:

| Gene | v1 identity | v2 identity | v2 masked codons | reference aa |
|---|---|---|---|---|
| UL27 (gB) | 99.3% | 99.35% | 904 | 904 |
| UL30 (pol) | 99.6% | 99.60% | 1235 | 1235 |
| US6 (gD) | 99.4% | 99.36% | 394 | 393 |
| **US8 (gE)** | **81.7%** | **99.21%** | 550 | 550 |
| **US7 (gI)** | **81.2%** | **98.56%** | 390 | 390 |
| **US4 (gG)** | **79.9%** (2097 nt, 66% padding) | **97.56%** | 238 | 238 |
| **RS1 (ICP4)** | **67.7%** | **99.48%** | 1299 | 1298 |
| **RL2 (ICP0)** | **64.2%** | **98.87%** | 775 | 775 |

Genes I predicted were unaffected (no clinical indels) are unchanged, which confirms the
mechanism. All 60 alignments pass the rectangular / in-frame / no-internal-stop invariants.

**Phase 6 — conservation.** 11,801 residue rows; within-type identity 97.2–99.9%; between-type
identity 37.5% (gG) to 90.8% (UL30); **283 windows of ≥9 aa conserved in both types** and
identical between them, the cross-protective candidate set. gG/US4 is correctly exposed as a
non-target: HSV-1 gG occupies only 238 of 698 alignment columns.

**Phase 7 — IEDB.** 4,806 records fetched; **606 unique human-host, positively assayed peptides**
mapped to priority-gene reference proteins (289 HSV-1, 317 HSV-2), best represented in gD and gB.

**Phase 8a — class I prediction.** ~24,000 9/10-mers across 40 reference proteins × 27-allele
HLA class I panel, with conservation joined per peptide. Binder rates 22–38% per gene.

**Phase 8b — validation against known biology (the substitute for wet-lab work).** MHCflurry
recovers curated IEDB class I epitopes with **mean AUC 0.809 (median 0.892) and 6.0× top-decile
enrichment** over 27 gene x type groups spanning 17 genes; **HSV-1 mean AUC 0.931**, HSV-2
0.600 (HSV-2 class I data in IEDB is sparser and CD4-weighted). Conservation alone scores **AUC 0.565** — near chance, which is the
direct evidence that conservation is not a proxy for immunogenicity and that the proposal's
two-filter design (predicted **and** conserved **and** not positively selected) is the right one.

**Correction A4, quantified retrospectively.** Applying Benjamini–Hochberg to the 62 old Stage 2
runs: of **5,126 sites reported at p<0.10, only 1,332 survive at q<0.05** — a 74% reduction — and
19 of the 46 surviving analysis units are the invalid pooled cohorts. Written to
`07_qc_reports/hyphy/retrospective_fdr_on_old_stage2.tsv`.

**Amendment to A10.** The cached `hyphy/hyphy` image does report **HYPHY 2.5.103(MP)**, so the
README's version claim was accurate; the real defect was that the image was never pinned and the
version was never recorded. `hyphy_runner.py` now writes version and image id to
`14_selection/gard_provenance.json` on every stage.

**Two defects found while re-running, both fixed in `hyphy_runner.py`:** the HyPhy container's
working directory is `/`, not the `/data` mount, so every relative path silently failed; and the
low-divergence per-type cohorts trigger HyPhy's `ComputeBranchCache` internal error, which needs
`ENV=TOLERATE_NUMERICAL_ERRORS=1`. Both would have blocked any per-type re-run.

**Blocking prerequisites for the user:**
1. **NetMHCpan-4.1 / NetMHCIIpan-4.x / BepiPred-3.0 / DiscoTope-3** require academic registration
   (DTU HealthTech) and cannot be installed from here. MHCflurry 2.3 substitutes for the CD8 arm
   and is already running, but **CD4 class II, B-cell linear and B-cell conformational prediction
   cannot proceed without those licences.** Register, download, and the integration step will
   consume their outputs.
2. **Keep Docker Desktop running** while phases 4–5 execute. They resume from cache if interrupted.

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


---

## Addendum — 2026-11-19, phases 9, 10 and 13 complete

Everything that does not depend on the selection analysis is now done. The two HyPhy-gated
phases (5 and 5b) remain. The three predictors once recorded as licence-gated are all
complete: BepiPred-3.0 and DiscoTope-3.0 were never gated and were run locally, and
NetMHCIIpan-4.3 was run on the DTU web server.

**Population coverage (new, phase 9b).** 27-allele class I panel over 239 IEDB populations:
World coverage 99.71%, median 97.5% across populations, 73.6% of populations above 90%. A
greedy minimal panel reaches 90.58% World coverage with one peptide (`ILDNPGFVTF`) and
99.22% with two (`+EVIDKINAK`). Three caveats are recorded in the methods draft and the
step log: HLA-C is absent from the panel so the figure is conservative; the metric saturates
once every panel allele is hit, so it discriminates panels rather than shortlist depth; and
populations with sparse frequency records in the IEDB dataset score low for that reason and
not biologically (Central America, 7.76%). The published Bui et al. (2006) formula was
reimplemented rather than the downloaded tool executed, and the downloaded data file was
read with a restricted unpickler that admits only `collections.OrderedDict`.

**Benchmark reproduced after the per-allele fix (phase 8b).** Re-running step15 to predict
per allele rather than best-allele-only changed no benchmark conclusion: affinity mean AUC
0.809, median 0.892, 6.03-fold top-decile enrichment over 27 scorable gene x type groups;
conservation alone 0.565 and 1.35-fold, i.e. near chance. That contrast is the result that
justifies treating conservation and predicted binding as two independent filters.

**Shortlist is A2-skewed (new observation, phase 9a).** With per-allele predictions, the
top 15 cross-protective candidates are all HLA-A*02:03. This is a known property of class I
predictors on hydrophobic 9-mers, not a biological finding, and is now stated as such in the
methods draft. The per-allele coverage analysis is what establishes that the shortlist as a
whole is not A2-restricted.

**Supplements (phase 10c), generated from the QC logs rather than transcribed.** The
data-flow table resolves the counts the curation protocol asked for: 1,305 records retrieved,
417 passing genome filters, 11,504 CDS extracted, 11,340 passing CDS QC, 4,163 unique
sequences after deduplication representing 12,054 isolates across 558 distinct accessions
and 40 gene x type cohorts. `S3_removed_sequences.tsv` has 1,329 rows: 888 genome-level
exclusions, 164 CDS QC failures, 80 type-verification disagreements, 34 occupancy-masking
column drops, and 163 alignment outliers that were **flagged and retained**, not removed.

**Workflow (phase 10b).** The `Snakefile` declares 16 stage rules plus partial targets, with
stages 1-7 as frozen inputs rather than rules, because re-downloading from NCBI would change
the cohort. Its rule graph and all 15 script references were checked; it has **not** been
dry-run, because snakemake is not installed in the environment the pipeline ran in. That is
stated in `requirements.txt` next to the commented-out pin rather than left implicit.

**GARD settles the pooled-cohort defect (A6/A7).** All 40 cohorts completed. US6/gD — the
gene for which the earlier pooled HSV-1+HSV-2 run returned a RECOMBINANT verdict — reports
**no breakpoint in either type** (dAICc = 0.0). A6 predicted exactly this: those verdicts
were the HSV-1/HSV-2 species split being read as recombination. Across the screen, 16 of 40
cohorts carry at least one real within-type breakpoint, yielding 37 partitions; with the 24
unpartitioned cohorts that is 61 analysis units.

Two cohorts are flagged rather than smoothed over: **UL44/HSV-1** is retained at a marginal
dAICc of 3.70 (partitioning is the conservative error here), and **RS1/HSV-1** is a
**truncated search** stopped by decision at its confirmed two-breakpoint model after ~16 h
and 5,215 models, recorded in the `search_truncated` column and in the first line of its
`.done` marker. Methods §7 states both.

### What is left

Everything that was computable in this environment is done. Two items remain, one of them
by decision:

| Item | Blocked on | Note |
|---|---|---|
| ESM-2 zero-shot | — | dropped by decision |

Nothing else is outstanding. The CD4 class II arm, the last licence-gated item, was
completed on 2026-11-21 via the DTU web server (see below).

**Resolved 2026-11-21 — the B-cell arm was never actually blocked.** An earlier audit
recorded BepiPred-3.0 and DiscoTope-3.0 as requiring the same DTU academic registration as
NetMHCIIpan. That was wrong: both are published as source under academic non-commercial
licences (BepiPred: DTU Health Tech terms; DiscoTope: CC BY-NC 4.0). Both were installed
and run locally on GPU. See methods §16, `09_scripts/step28_bcell_predictors.py`, and
`13_epitopes/bcell_predictor_provenance.tsv`. Outputs:
`13_epitopes/bcell_linear_scores.tsv` (23,184 residues, all in register) and
`13_epitopes/bcell_conformational_scores.tsv` (3,145 residues, 0 unmapped, 0 mismatches).
Neither predictor's published threshold behaves as documented on this data — BepiPred's
flags 37.3% of all residues, DiscoTope's flags none — so both are used as continuous
annotations and the thresholds are reported as unmet rather than dropped.

### Final selection results

| Level | Tested | Called under selection |
|---|---|---|
| Gene-wide (BUSTED) | 61 runs | 32 runs at q < 0.05 (= 24 of 40 cohorts, 18 of 20 genes) |
| Site (FEL, MEME, FUBAR) | 23,155 codons x 3 = 69,465 rows | 92 by >= 1 method |
| Site, primary (>= 2 methods) | 92 candidates | **2** |
| Branch (aBSREL) | 8,789 branches | 75 at Holm p < 0.05 |

The two concordant sites are HSV-1 RL2 (ICP0) residue 308 (FUBAR;MEME, q = 0.0246) and
HSV-1 US7 (gI) residue 155 (FUBAR;MEME, q = 0.0051). **Both fail the variant-support screen
of 2026-11-21 and must not be cited as positive-selection findings** — see finding V12 and
methods section 8.1. The defensible summary is that no site is both multi-method concordant
and robustly supported, while 76 of the 92 single-method sites are well supported. Nine of the 12,771 ranked peptides
overlap them and are excluded; the final shortlist is 9,450 rows (8,070 distinct peptides),
3,202 rows (2,120 distinct) of them cross-protective. Neither residue lies in one of the four structurally mapped proteins, so
`17_structures/structure_residue_annotation.tsv` has `positively_selected=False` for all
3,145 mapped residues — a structure-coverage limit, not a finding about those residues.

### Verification pass, 2026-11-21

Every number in the methods draft, the README and this file was re-derived from the source
files. What held, and what did not:

**Verified exactly.** Benjamini–Hochberg q-values recomputed independently for all 122
FEL/MEME families and the 61-unit BUSTED family (max absolute difference 5.3e-14 and
2.4e-16); aBSREL Holm-corrected p-values consistent with the `significant` flag on all
8,789 branches; the 92-site union and 2-site concordant set reproduced exactly from the
per-method flags (94 positively selected rows, 2 sites shared by MEME and FUBAR); all 37
GARD partition spans tile their alignments exactly, on codon boundaries, with the parent
sequence set, and all 37 partition-tree content hashes match their current alignments; all
40 whole-cohort trees have matching tip counts, a completed `.iqtree` and bootstrap labels;
masked codon counts equal the analysis alignment lengths for all 40 cohorts; 4,037 dropped
codon columns, 163 flagged outliers, 283 conserved windows, 688/4,118 IEDB records, 239
populations, 3,145 structural residues / 7 chains / 33 patches, S4 sums (4,163 and 12,054),
and the retrospective 5,126 -> 1,332 figures all match their files. The `RMLGDVMAV` positive
control sits at rank 11 of 3,202 in the cross-protective shortlist.

**Corrected.** Nine discrepancies were found and fixed in the documents:

| # | Finding | Resolution |
|---|---|---|
| V1 | `18_population_coverage/` was computed on 2026-11-19 from the *provisional* shortlists; step23 reads the shortlist files directly and had not been re-run | step23 re-run; `shortlist_primary` 8,075 -> 8,070 distinct peptides, coverage figures unchanged (they are allele-set driven) |
| V2 | The CDS counts 11,504 / 11,340 / 164 describe the QC log, not the dataset: **12,054** sequences were carried forward. 697 are absent from the QC log (a later extraction top-up), 19 carry `FAIL` | Methods §2 rewritten to cite 12,054 and explain the gap; §3's "11,340 -> 4,163" corrected to "12,054 -> 4,163". All 686 traceable extras are ATG-start, in-frame, valid-stop, 92.9-102.4% of reference, so the data is sound and the bookkeeping was not |
| V3 | 8 `FAIL` sequences reach the final analysis, all HSV-1 US7, all failing only the 0.95-1.05 length-ratio window | Methods §2 now states this and the evidence that it is a real in-frame +84 nt gI polymorphism (`07b`, `07c`), and that the length-ratio rule is a mis-annotation heuristic, not a correctness test |
| V4 | Methods §1 claimed the 229-genome HSV-2 tier-2 set was "analysed separately". It is not: **178 of the 286 HSV-2 accessions (62%)** in the final cohorts are tier-2 genomes, readmitted from the excluded pool for 1-2% ambiguity | §1 rewritten to state the pooling, the 62% figure, and that HSV-2 inference therefore rests on lower-quality assemblies than HSV-1. This is the most consequential of the nine |
| V5 | S1/§1 "888 excluded, 417 passed" lumps 29 flagged-but-kept records in with 859 exclusions, and cannot be reconciled with the 649 genomes that entered extraction | §1 now gives 417 PASS / 29 FLAGGED / 859 EXCLUDED, notes that 229 of the excluded were readmitted (so 630 are excluded from all analysis), and says `02_filtered_genomes/` defines what entered extraction |
| V6 | README said "32 genes" | 32 is analysis units; = 24 of 40 cohorts, 18 of 20 genes. Fixed in the README and stated in methods §8 |
| V7 | "9,450 peptides" / "3,202 cross-protective" are row counts; a row is peptide x virus type | Distinct peptide counts are **8,070** and **2,120**; both forms now given, with the distinct counts named as the ones to quote |
| V8 | Methods §12 said "17 scorable gene x type groups" and "HSV-1 mean AUC 0.920" | There are **27** groups (17 HSV-1 + 10 HSV-2) spanning 17 distinct genes; HSV-1 mean AUC is **0.931**, HSV-2 0.600. `genes_scored=17` in `benchmark_summary.tsv` counts genes, not groups |
| V9 | The 37 partition trees use a fixed GTR+F+G4 with no ModelFinder and no bootstrap, unlike the 40 whole-cohort trees of §6 | Methods §7 now states this and why it is sound (HyPhy does not read branch support) |

**Found on the second pass (2026-11-21), after the nine above.** Re-deriving the QC record
over the full carried-forward set turned up three further items, one of which changes what
this study can claim:

| # | Finding | Resolution |
|---|---|---|
| V10 | `step07_gaps_abc.py` rebuilds the RL2/RS1 QC-passed FASTAs and its comment claimed "identical Step 7 criteria", but it omitted the **ambiguity** test (`%N < 0.5`) and the length-ratio window. Those two genes passed a weaker test than the other eighteen, and **81 sequences with 0.5-16.8% ambiguity (median 2.2%) reached the final cohorts** | Script corrected to apply the full criteria through a shared `qc_ok()`. The cohort is **not** re-filtered: masking removes the repeat regions where the ambiguity sits, so in the analysed alignments ambiguous bases are 0.007-0.042% of cells and the worst column has `N` in 7.0% of sequences, while re-filtering would shrink the two smallest HSV-2 cohorts and invalidate a 16-hour GARD search. Documented with the measured impact in methods 2.2, and the cohort is pinned by accession in S2 |
| V11 | 16.2% of the primary shortlist (1,527 of 9,450 rows) are low-complexity peptides from the ICP0/ICP4 tandem repeats -- `SSSSSSSSA`, `GAGGGGGGR` -- which are perfectly conserved, score as binders, are not credible epitopes, and cannot be assigned a unique position | `low_complexity` flag added to step21 with the entropy / max-residue-fraction / distinct-residue rule; counts logged. Flagged rather than dropped so the shortlist size is not changed silently; one line in step21 excludes them instead. Only 1 of the top 50 is affected |
| **V12** | **The two concordant positively selected sites rest on 4 and 2 isolates of ~290.** FUBAR posteriors are 0.9560 and 0.9523 against a 0.95 cut; MEME beta estimates are 1,024 and 1,842, the unbounded values a near-invariant column produces; FEL supports neither (p = 0.169, q = 0.637). RL2/308's variants include **OR771687, which fails the V10 ambiguity test**; US7/155's only two variants are **PP379031 and PP379033, adjacent accessions from one submission batch**, and one is a 3-nucleotide codon change seen once | `step26_site_robustness.py` added: it annotates all 92 sites with variant sequence and isolate support, distinct variant codons, maximum nucleotide change, QC status of supporting sequences and accession adjacency, under a threshold fixed before inspecting the outcome. Both concordant sites are `NOT_ROBUST`; 76 of 92 single-method sites are `ROBUST`. `selected_sites_final.tsv` is left unchanged so the pre-registered rule's output stands, and `site_robustness.tsv` qualifies it. Methods section 8.1 now states that no site is both concordant and robust |

V12 is the finding that matters. The pre-registered ">= 2 methods" rule was intended as the
conservative choice, and on these data it did the opposite: it selected the two sites most
likely to be artifacts, because MEME fires on singleton non-synonymous changes at invariant
columns and FUBAR's posterior grazes the threshold there. Concordance and variant support
are independent checks; this study now applies both, and reports that nothing passes both.

**Also noted, not defects.** FEL calls 154 sites at q < 0.05 but none with beta > alpha, so it
contributes no positively selected site — these cohorts carry pervasive purifying selection,
which methods §8 now says. Ten IUPAC ambiguity characters beyond `N` (R, Y, S, B, D, M)
occur across the 4,163 analysis sequences, in 5 HSV-2 cohorts; HyPhy handles IUPAC codes and
the quantity is negligible. Population coverage is identical across all three shortlist
tiers because every one of the 27 panel alleles is hit by all three, which methods §13
already explains. `S3_removed_sequences.tsv` correctly labels the 163 outlier rows
"sequence (flagged, retained)" even though the file is named for removals. Minor cosmetic
items left alone: `cohort_tier=Tier` is truncated in FASTA headers, 2 shortlist rows
duplicate on (peptide, gene, virus_type), and `gard_provenance.json` records `cpus: 4`
rather than the per-worker setting actually used.

### The one error class behind six bugs

Every caching and validation bug found during execution was the same mistake: treating
**"a file exists"** as **"the run finished"**. HyPhy rewrites its JSON after every
breakpoint stage and IQ-TREE writes `.treefile` before bootstrapping ends, so a killed run
leaves a file that parses. Six instances were found and fixed (GARD JSON cache, tree cache,
step17 tree validation, `--allow-partial-gard` pairing, partition-tree reuse, HyPhy
filename-keyed grid caches). Each now keys on an artifact written only on completion:
`.done` for HyPhy, `.iqtree` for IQ-TREE, a content hash `.alnmd5` for partition trees.
This is the one correction most worth carrying into any future stage of this project.

---

## Addendum — 2026-11-21, repository cleanup and data dictionary

### What was removed

Recoverable from git history at commit `4f56731`:

| Path | Size | Why |
|---|---|---|
| `hsv_curated_dataset/06_alignments/` | 53 MB, 519 files | Withdrawn v1 alignments. The two virus types were pooled before aligning, so codons are out of register. Superseded by `06_alignments_v2/`. No live script read it. |
| `_archived_discarded/` | 24 MB, 27 files | Pre-screening prototypes and superseded audit scratch. |

Not recoverable — never tracked in git:

| Path | Size | Why |
|---|---|---|
| `hsv_curated_dataset/07_qc_reports/hyphy/_tmp/` | 410 MB | HyPhy scratch and cache dumps. Pure intermediate; already matched by the `_tmp/` and `*.cache` entries in `.gitignore`. |
| `.../hyphy/test_UL27_HSV1_fubar*.json` | 16 MB | Two stray single-gene probe runs from pipeline development. |
| `.../hyphy/tier_{a,b,c}/`, `10_hyphy_outputs/`, `logs/` | 0 | Empty directories and a dead `DEPRECATED` marker. |
| `09_scripts/fix_strings.py`, `find_datasets.py`, `__pycache__/` | <1 MB | A one-off string-repair hack, an exploratory helper, and bytecode. |

### What was retained deliberately

- **`07_qc_reports/hyphy/stage1/`, `stage2/`, `gard/` — 1.3 GB of withdrawn v1
  HyPhy JSON.** Read by nothing. Retained because those runs took roughly 16
  hours and deciding what analysis output leaves the record is not a cleanup
  decision. Flagged here as the open question.
- **The 14 small v1 summary tables and logs** in `07_qc_reports/hyphy/`, and the
  four `11_*` v1 summaries in `07_qc_reports/` — these are what this document
  cites, so they are the part of the withdrawn arm that matters.
- **`09_scripts/deprecated_v1/`** — the 10 withdrawn v1 selection scripts, moved
  into a subdirectory rather than deleted. They are the evidence for the
  corrections in Part A. Do not run them.

### Stale references found and fixed

Both were discovered by grepping for the deleted paths, not by assumption:

1. `hyphy_runner.py:89` pointed at `10_hyphy_outputs/DEPRECATED`, which no longer
   exists. Repointed at `07_qc_reports/hyphy/`.
2. `step08_dedup_and_step09_align.py` still writes to the deleted
   `06_alignments/`. The script is retained because its **step-8 half** produced
   `07_qc_reports/08_dedup_log.tsv`, a frozen input that cannot be regenerated
   independently of it. Its code was left unchanged — it is frozen-stage
   provenance — but a header now states that the step-9 half is withdrawn and
   must not be run.

### New: `08_documentation/DATA_DICTIONARY.md`

Every directory and file in the repository, with contents, producing script and
consuming script. Files with no consumer are marked terminal explicitly rather
than left blank.

### New: `09_scripts/step30_provenance_audit.py`

Enforces the dictionary. Walks every data file and fails, naming the gaps, if any
has no entry. Read-only. Available as `snakemake provenance_audit`; **not** in
`rule all`, because it audits the tree rather than producing a result and must be
able to fail without blocking a run.

Current state: **2,511 files audited, 2,511 accounted for, 0 unaccounted.**
Tool binaries under `bin/` are excluded, and the audit prints that exclusion and
its reason rather than omitting it silently.

**The first version of this audit was wrong in a way worth recording.** It
accepted any file sitting in a documented *directory*, so a canary file dropped
into `11_conservation/` passed. It reported 2,511/2,511 while being incapable of
failing. The fix was to restrict directory-level accounting to an explicit
`BULK_DIRS` list — the per-gene and per-allele directories where describing 610
JSON files collectively is the right granularity — and to require every other
directory to name each of its files. Tightening it surfaced 7 genuine gaps that
the loose version had hidden: the two `.gbk` reference annotations, and all five
files in `02_filtered_genomes/`, which had been described as "five FASTA files"
without naming them or distinguishing `_filtered` from `_flagged` from
`_filtered_tier2`. A passing audit that cannot fail is worse than no audit, and
this one passed for a while before it could.

---

## Addendum — 2026-11-21, three Windows line-ending defects

Recorded here because all three were silent: each returned success while doing
the wrong thing, and none would have been noticed by reading output alone.

### 1. Corrupted allele-frequency data

`18_population_coverage/population_genotype_map.p` could not be loaded; `step23`
aborted. Git's heuristic did not classify the pickle as binary, so LF became CRLF
on checkout and 380,776 CRLF pairs were written into it — 4,000,753 bytes on disk
against 3,619,977 committed, a difference of exactly the CRLF count. The `GLOBAL`
opcode read `collections\r` / `OrderedDict\r`.

Restored from the committed blob (SHA-256 verified identical) and confirmed by
regeneration: every coverage number reproduces exactly as previously published.
All 136 other binary files were checked and are intact. `.gitattributes` now pins
the binary formats.

### 2. Figure 4 was not reproducible

`step27_figures.py` seeded its jitter with `hash((verdict, cohort))` under a
comment reading "deterministic jitter". Python randomises string hashing per
process (`PYTHONHASHSEED`), so the seed changed on every run and the figure's
point positions moved each time it was regenerated — a reviewer regenerating
Figure 4 would get a different picture. Now seeded from `zlib.crc32`; verified by
running the script twice in separate processes and comparing SHA-256 of the PNG.

### 3. Git's stdin, during the staged upload

While splitting the initial push into size-bounded commits, `subprocess.run(...,
text=True)` translated `\n` to `\r\n` on stdin, so every path handed to
`git update-index --index-info` arrived with a trailing carriage return and
matched nothing. Git printed `Ignoring path ...` to stderr and **returned exit
code 0**, so the first attempt produced ten commits whose trees were byte-identical
to the base. It was caught only by comparing the final tree SHA against the
intended one. Fixed by passing raw bytes.

**The common thread is that exit code 0 meant nothing in all three cases.** The
same mistake appeared twice more in shell: `if git push ... | tail -4` tests
`tail`'s status, not git's, and so reported a failed push as successful — once on
the first push attempt and again in the staged-upload loop, which printed
"push 5 OK" directly beneath push 5's own error. Check the artifact, not the
status code.


## Addendum — 2026-11-21, three defects in the conservation layer

Found while assembling the candidate-targets section of the manuscript, which is the first
thing in this project to read the shortlist as a *list of peptides someone would synthesise*
rather than as a table of counts. Two of the five peptides in the published minimal panel
turned out not to exist in HSV-1.

### How it surfaced

`AEYGGSFRF` (gB) and `AAFDGGSTA` (UL47) were flagged cross-protective and sat in the greedy
minimal panel. Neither occurs anywhere in the HSV-1 reference protein. Both reference
sources — the GenBank translations and `13_epitopes/bcell_inputs/reference_proteins.fasta` —
agree, so it was not a FASTA problem:

```
HSV-2 gB 373-381   AEYGGSFRF
HSV-1 gB 373-381   DEMLRSEYG      "AEYGGSFRF" in HSV-1 gB -> False
```

`11_conservation/per_residue_conservation.tsv` said HSV-2 gB residue 373 was **G**; the
reference protein says **A**. The table's residue 373 was the protein's residue 376.

### Defect 1 — reference residues numbered on the masked alignment

`step13_conservation_layer.py` numbered reference residues by counting non-gap positions in
the **masked** combined alignment. Masking removes columns, so a reference residue in a
dropped column was never counted and every number after it was short by one, drifting
further with each drop.

Audited by checking every row's residue letter against the reference protein at the number
the row claimed: **7,343 of 22,955 rows were wrong** before the fix, across 17 of 40 gene ×
type combinations — RS1 HSV-2 1,125, UL30 HSV-2 1,122, UL27 HSV-1 817 and HSV-2 802, RL2
HSV-2 653, UL47 HSV-2 607. After the fix, 0. Comparing the two tables, 8,072 residue numbers
changed (1,860 HSV-1, 6,212 HSV-2).

Fixed by numbering on the unmasked combined alignment and indexing through the retained
column list that `step09c` already records in `masking_report_v2.tsv` — the authoritative
record, reconstructed from `dropped_codon_ranges` rather than inferred by matching column
content, which is ambiguous wherever two columns are identical and so unusable in the
tandem repeats of RL2 and RS1.

**`conserved_both_types` and `between_type_identical` changed in zero columns.** The
column-level conservation was right all along; what was wrong was which reference residue
each column was labelled with. Every downstream join therefore picked up the correct column
property and attached it to the wrong residue.

### Defect 2 — flags set from a partial residue set

With the numbering fixed, 9 non-existent peptides survived. `step15`'s
`annotate_conservation` evaluated `all(...)` over however many of a peptide's residues had
a conservation row. Masked columns have none, so `VVSPRQLAL` (RS1, HSV-2 108–116) was
flagged fully identical between types on the strength of residues 109–116: residue 108's
column was masked and was never examined. `residues_with_conservation_data` was already
recorded and never consulted.

Fixed: a flag that cannot be evaluated for every residue of the peptide is `False`.

### Defect 3 — consensus residues answering a reference-residue question

Three peptides still survived. `step13` set `between_type_identical` from
`hsv1_consensus_aa == hsv2_consensus_aa` — the cohort consensus. Candidate peptides are cut
from *reference* proteins, so the claim the flag was used to support was about reference
residues. The two disagree in **55 of 11,801 columns** (45 where the consensus agrees but
the reference residues differ, 10 the reverse). At UL10 column 105 the HSV-1 reference is V
and the HSV-2 reference is A, while both cohort consensuses are A — which is how
`LLLAAAVYA` (UL10, HSV-2 101–109) became cross-protective although HSV-1 UL10 does not
contain it.

Fixed: `between_type_identical` keeps its consensus meaning and is documented as such;
`between_type_reference_identical` was added and is what the shortlist uses;
`conserved_both_types` requires both.

### What the three had in common

Each was a quantity computed correctly and then used to answer a slightly different
question — residue numbers from one alignment applied to another, `all()` over a subset read
as `all()` over the whole, consensus identity read as reference identity — with nothing
asserting the difference. That is the same error class as the six bugs recorded above, and
the remedy was the same: assert the thing you actually mean.

Two assertions now stand in the pipeline, either of which would have caught all three:

- `step13` aborts unless every emitted `(ref_residue, ref_aa)` pair matches the GenBank
  translation at that number. This is the register check that `step28` and `step29` already
  applied to the B-cell and class II predictors; the conservation layer was the one place it
  was missing.
- `step21` confirms every cross-protective peptide against the other type's reference
  protein directly, and that check does not pass through the alignment at all, so it cannot
  fail the same way. Peptides it rejects go to `16_shortlist/excluded_not_in_other_type.tsv`
  with their row intact rather than being dropped silently. It currently rejects 2 — genuine
  indel-inside-the-window cases, where every residue corresponds but the contiguous peptide
  exists in one type only.

### What changed in the results

Re-ran `step15`, `step19`, `step21`, `step22`, `step23`, `step24`, `step26` and `step27`.
`step17` was not re-run: HyPhy does not read this table.

| Quantity | Before | After |
|---|---|---|
| Primary shortlist | 9,450 rows / 8,070 peptides | 9,488 / 8,010 |
| Cross-protective | 3,202 rows / **2,120** peptides | 3,207 / **1,729** |
| Low-complexity primary rows | 1,527 (16.2%) | 1,460 (15.4%) |
| Rows the selection filter actually removes | 9 | 5 |
| Conserved windows ≥ 9 aa | 283 | 281 |
| Conserved exposed patches | 33 | 32 |
| Minimal panel | **4 peptides** | **5 peptides** |
| Panel members | ILDNPGFVTF, EVIDKINAK, AEYGGSFRF, AAFDGGSTA | ILDNPGFVTF, EVIDKINAK, AMAAGRFGW, AQTGVRVAL, AAQYRNAVV |
| Low-complexity barred / eligible | 1,404 / 1,817 of 2,120 | 1,307 / 1,559 of 1,729 |
| `RMLGDVMAV` control rank | 11 of 3,202 | 13 of 3,207 |
| Conservation-alone AUC / enrichment | 0.565 / 1.35× | 0.566 / 1.38× |

**Unchanged:** the entire selection screen (154 FEL sites, 92 candidates, 2 concordant, 76
ROBUST / 7 WEAK / 9 NOT_ROBUST, BUSTED 32 of 61, aBSREL 75 of 8,789); the affinity benchmark
(mean AUC 0.8088, median 0.8915, HSV-1 0.9313, HSV-2 0.6003, 428 positives, 27 groups); all
population-coverage percentages (99.71% World, median 97.5%, 73.6% ≥ 90%, 239 populations,
27 alleles) and the whole class II bracket; 3,145 structure residues and 1,099 exposed; and
the within- and between-type conservation percentages.

**gB `SRVEAFHR` survives.** It is still the longest fully conserved, fully epitope-covered
exposed patch, in all three chains of 2GUM at 100% minimum identity, and it is contiguous
and identical in both types — HSV-1 257–264, HSV-2 252–259. The manuscript's strongest
structural claim is unaffected.

### The uncomfortable part

The defect put two non-existent peptides into the headline deliverable of the vaccine-design
arm, and it was caught only when someone asked what the peptides *were* rather than how many
there were. A count can be wrong in ways that reading the list cannot be. That is worth
remembering for the next table of candidates: check a few by hand against the source
sequence, because the assertion you did not write is the one you need.
