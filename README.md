# Herpes Simplex Virus (HSV-1 & HSV-2): Evolutionary Constraint and Vaccine-Epitope Mapping

[![Status](https://img.shields.io/badge/Status-Corrected_pipeline_complete-success.svg)](CORRECTIONS_AND_PENDING_WORK.md)
[![Alignments](https://img.shields.io/badge/Alignments-rebuilt_(FAMSA_codon--aware)-success.svg)](#)
[![Selection](https://img.shields.io/badge/Selection-complete_(305%2F305_runs%2C_per_type)-success.svg)](#)
[![Epitopes](https://img.shields.io/badge/Epitopes-CD8_%2B_IEDB_done%3B_CD4%2FB--cell_licence--gated-yellow.svg)](#)

Conservation, selection-pressure and vaccine-epitope analysis of HSV-1 and HSV-2 clinical
isolates, across the glycoproteins, replication machinery and regulatory tegument genes.

> **Read [`CORRECTIONS_AND_PENDING_WORK.md`](CORRECTIONS_AND_PENDING_WORK.md) first.**
> An audit on 2026-11-19 found that the alignments underpinning the earlier results were not
> true multiple sequence alignments, and that the recombination and selection results had been
> computed on pooled HSV-1 + HSV-2 data, after rather than before recombination screening, with
> no multiple-testing correction. Those results have been withdrawn and the pipeline is being
> re-run. That document records every defect, what was kept, and what is still outstanding.

---

## Current status

| Stage | State |
|---|---|
| Genome acquisition, filtering, ANI type verification | **Complete and validated** |
| Metadata curation | Complete; feasibility limits documented (`03_metadata/metadata_feasibility_statement.md`) |
| CDS extraction and per-gene QC | **Complete and validated** |
| Deduplication | **Complete** |
| Codon-aware alignment | **Rebuilt** (`06_alignments_v2/`) — superseded `06_alignments/` |
| Masking and reference coordinate maps | **Rebuilt** (`06_alignments_v2/masked/`, `coordinate_maps/`) |
| Per-type phylogenies | **Complete** (40/40, all with finished bootstraps; 4 needed IQ-TREE `-safe`) |
| Recombination screening (GARD, per type) | **Complete** (40/40; 16 cohorts with breakpoints -> 37 partitions; no breakpoint in US6 in either type) |
| Selection (BUSTED / FEL / MEME / FUBAR / aBSREL, per type, breakpoint-aware) | **Complete** (305/305 runs over 61 analysis units) |
| FDR correction and cross-method concordance | **Complete** (`14_selection/`) — 32 of 61 analysis units gene-wide (24 of 40 cohorts, 18 of 20 genes), 75 branches, 92 sites by ≥1 method (76 robust), 2 concordant (**both non-robust**) |
| Site variant-support screen | **Complete** (`14_selection/site_robustness.tsv`) — no site is both multi-method concordant and robustly supported |
| Per-residue conservation layer | **Complete** (`11_conservation/`) |
| IEDB curation | **Complete** (`12_iedb/`) |
| CD8 class I prediction | **Complete** (`13_epitopes/`) |
| Predictor validation against IEDB | **Complete** (`15_benchmark/`) |
| Candidate integration and ranking | **Complete** (`16_shortlist/`) — selection filter applied; 9,450 rows / 8,070 distinct peptides primary, 3,202 rows / 2,120 distinct cross-protective |
| CD4 class II, B-cell linear, B-cell conformational | **Blocked** on academic licences (DTU HealthTech) |
| Structure mapping (2C36 gD, 2GUM gB, 3M1C gH-gL) | **Complete** (`17_structures/`) |
| HLA class I population coverage, 239 populations | **Complete** (`18_population_coverage/`) |
| Pinned environment and Snakemake workflow | **Complete** (`requirements.txt`, `Snakefile`) |
| Supplements S1-S5 (data flow, accessions, removals, cohorts, versions) | **Complete** (`19_supplements/`) |
| Methods draft | **Rewritten and filled in** (`08_documentation/methods_draft.md`) — no section left pending |
| Results draft | **Written** (`08_documentation/results_draft.md`) — headline is the site-level negative result |
| Main-text figures | **Complete** (`20_figures/`, 6 figures, PNG + vector PDF) |
| ESM-2 zero-shot | Dropped by decision |

The corrected pipeline is complete as of **2026-11-20**, and every number in it was
re-verified against the source files on **2026-11-21** (see "Verification pass" in the
corrections document for what that found and changed). Results live in `14_selection/`
(selection), `16_shortlist/` (candidates), `17_structures/` (structural mapping) and
`19_supplements/` (S1-S5); the methods draft states every number and the file it comes from.

**The v1 outputs are withdrawn, not superseded in place.** `06_alignments/` and
`10_hyphy_outputs/` each carry a `DEPRECATED_DO_NOT_USE.md` explaining why and naming the
replacement. The earlier "342/342 runs, 100% complete" badge was both arithmetically wrong
(the real total was 270) and describing a pipeline whose inputs were invalid. Two results
carry explicit caveats that should travel with any citation: the marginal UL44/HSV-1
breakpoint (dAICc 3.70, retained as the conservative choice) and the truncated RS1/HSV-1
GARD search (stopped by decision at two confirmed breakpoints). The licence-gated
predictors (CD4 class II, B-cell linear, B-cell conformational) were never run, so this work
covers the CD8 arm only.

---

## What the corrections changed, so far

**The alignments.** Step 9 of the original pipeline aligned each sequence *pairwise against the
reference protein* and projected codons into that sequence's own coordinate frame. Indel
placement was never reconciled between sequences, so every column downstream of a clinical
indel was out of register. Rebuilt with a genuine progressive codon-aware MSA:

| Gene | old within-HSV-1 identity | rebuilt | masked codons | reference aa |
|---|---|---|---|---|
| UL27 (gB) | 99.3% | 99.35% | 904 | 904 |
| UL30 (pol) | 99.6% | 99.60% | 1235 | 1235 |
| US6 (gD) | 99.4% | 99.36% | 394 | 393 |
| US8 (gE) | **81.7%** | **99.21%** | 550 | 550 |
| US7 (gI) | **81.2%** | **98.56%** | 390 | 390 |
| US4 (gG) | **79.9%** | **97.56%** | 238 | 238 |
| RS1 (ICP4) | **67.7%** | **99.48%** | 1299 | 1298 |
| RL2 (ICP0) | **64.2%** | **98.87%** | 775 | 775 |

Genes with no clinical indels were unaffected, which confirms the mechanism. Masked lengths now
land on the true protein lengths.

**Multiple testing.** Applying Benjamini–Hochberg retrospectively to the 62 old Stage 2 runs:
of **5,126 sites reported at p<0.10, 1,332 survive at q<0.05** (`07_qc_reports/hyphy/
retrospective_fdr_on_old_stage2.tsv`). Nineteen of the 46 surviving analysis units were pooled
HSV-1 + HSV-2 cohorts and are invalid regardless.

**Selection — the headline, stated carefully.** The screen is complete (305/305 HyPhy runs over
61 analysis units) and gene-level and branch-level signal is substantial: **32 of 61 runs** show
gene-wide episodic diversifying selection at q<0.05 (24 of 40 cohorts, 18 of 20 genes), and
**75 of 8,789 branches** at Holm-corrected p<0.05. At site level the result is negative, and
this is the claim to be careful with. 92 sites are called by at least one method and only **2**
by the pre-registered ≥2-method rule — and a variant-support screen added on 2026-11-21
([`site_robustness.tsv`](hsv_curated_dataset/14_selection/site_robustness.tsv)) shows **both of
those 2 rest on 4 and 2 isolates out of ~290**, with FUBAR posteriors a thousandth above the
0.95 cut, MEME β estimates of 1,024 and 1,842, no FEL support, one variant from a sequence that
fails CDS QC and the other pair from a single submission batch. **Neither is reported as a
positive-selection finding.** Meanwhile 76 of the 92 single-method sites *are* well supported
(median 44 isolates) but rest on one method. So: **no site in this study is both multi-method
concordant and robustly supported.** That is the expected outcome at 99.4%/99.7% within-type
identity, and it is the correct result to publish — the ≥2-method rule, meant to be
conservative, in fact selected the two most likely artifacts. See methods §8.1.

**Conservation.** 11,801 residue rows; within-type identity 97.2–99.9%; between-type identity
37.5% (gG) to 90.8% (UL30); **283 windows of ≥9 aa conserved in both types and identical
between them.** gG/US4 is confirmed as a poor cross-protective target — HSV-1 gG occupies only
238 of 698 alignment columns against HSV-2 gG.

**Validation against known biology.** MHCflurry recovers curated IEDB class I epitopes with
**mean AUC 0.809 (median 0.892), 6.0× top-decile enrichment** over 27 gene × type groups
spanning 17 genes; **HSV-1 mean AUC 0.931**, HSV-2 0.600. Conservation alone scores **AUC 0.565**,
near chance — direct evidence that conservation is not a proxy for immunogenicity, and that both
filters are needed. `RMLGDVMAV` (gB 559–567 in HSV-2 numbering), a published HLA-A\*02 gB epitope,
is recovered at rank 11 of 3,202 in the cross-protective shortlist and is
recovered and flagged as IEDB-known.

---

## Repository layout

**Every file in the tree below is described individually in
[`08_documentation/DATA_DICTIONARY.md`](hsv_curated_dataset/08_documentation/DATA_DICTIONARY.md)** —
contents, producing script, and consuming script. `step30_provenance_audit.py`
enforces that: it fails if any data file on disk has no dictionary entry.

```
vaccine project/
├── README.md
├── CORRECTIONS_AND_PENDING_WORK.md       # audit, corrections, pending work, execution log
├── Snakefile                             # 26 rules, reproduces stage 13 onward
├── requirements.txt
├── hsv_curated_dataset/
│   ├── 00_raw_downloads/                 # NCBI Datasets genome dumps        [not in git]
│   ├── 01_reference_genomes/             # NC_001806.2 (HSV-1 strain 17), NC_001798.2 (HSV-2 HG52)
│   ├── 02_filtered_genomes/              # quality-filtered genomes, flagged sets kept  [not in git]
│   ├── 03_metadata/                      # provenance, patient assignments, sampling-bias tables
│   ├── 04_cds_unaligned/                 # extracted CDS, all genes
│   ├── 05_cds_qc_passed/                 # ATG/stop/frame/%N validated, deduplicated
│   ├── 06_alignments_v2/                 # rebuilt codon MSAs, masked sets, coordinate maps,
│   │                                     #   HyPhy-ready alignments, per-type trees
│   ├── 07_qc_reports/                    # genome/CDS/dedup QC logs; hyphy/ = withdrawn v1 output
│   ├── 08_documentation/                 # methods, results, DATA_DICTIONARY
│   ├── 09_scripts/                       # pipeline step02-step30; deprecated_v1/ = withdrawn
│   ├── 11_conservation/                  # per-residue conservation, conserved windows
│   ├── 12_iedb/                          # curated IEDB epitopes and residue flags
│   ├── 13_epitopes/                      # class I, class II and B-cell predictions + raw output
│   ├── 14_selection/                     # per-type GARD, partitions, selection, FDR outputs
│   ├── 15_benchmark/                     # predictor benchmarking against IEDB
│   ├── 16_shortlist/                     # integrated, ranked epitope candidates
│   ├── 17_structures/                    # PDB structures, residue annotation, exposed patches
│   ├── 18_population_coverage/           # class I and class II coverage, minimal allele panel
│   ├── 19_supplements/                   # S1-S5 manuscript supplementary tables
│   ├── 20_figures/                       # fig1-fig6, PDF and PNG
│   └── bin/                              # NCBI datasets, IQ-TREE 2.4.0            [not in git]
```

Deleted in the 2026-11-21 cleanup, recoverable from git history at commit
`4f56731`: `06_alignments/` (withdrawn v1 alignments, out of codon register) and
`_archived_discarded/` (pre-screening prototypes).

---

## Reproducing the corrected pipeline

Order matters: alignment → QC/masking → trees → recombination → selection → FDR.

```powershell
python hsv_curated_dataset/09_scripts/step09b_rebuild_codon_msa.py
```
```powershell
python hsv_curated_dataset/09_scripts/step09c_mask_and_coordinate_map.py
```
```powershell
python hsv_curated_dataset/09_scripts/step10b_prep_and_trees.py
```
```powershell
python hsv_curated_dataset/09_scripts/step16_gard_per_type.py
```
```powershell
python hsv_curated_dataset/09_scripts/step17_selection_per_type.py
```
```powershell
python hsv_curated_dataset/09_scripts/step18_fdr_and_concordance.py
```

Epitope arm (independent of the HyPhy stages until integration):

```powershell
python hsv_curated_dataset/09_scripts/step13_conservation_layer.py
```
```powershell
python hsv_curated_dataset/09_scripts/step14_iedb_curation.py
```
```powershell
python hsv_curated_dataset/09_scripts/step15_tcell_class1_prediction.py
```
```powershell
python hsv_curated_dataset/09_scripts/step19_benchmark_predictors.py
```
```powershell
python hsv_curated_dataset/09_scripts/step23_population_coverage.py
```
```powershell
python hsv_curated_dataset/09_scripts/step21_integrate_and_rank.py
```
```powershell
python hsv_curated_dataset/09_scripts/step22_structure_mapping.py
```
```powershell
python hsv_curated_dataset/09_scripts/step24_build_supplements.py
```

```powershell
python hsv_curated_dataset/09_scripts/step25_qc_reconciliation.py
```

```powershell
python hsv_curated_dataset/09_scripts/step26_site_robustness.py
```

```powershell
python hsv_curated_dataset/09_scripts/step27_figures.py
```

Every stage caches completed work and is safe to re-run. Caching keys on an artifact the
tool writes only when it finishes — `.done` for HyPhy, `.iqtree` for IQ-TREE, a content hash
`.alnmd5` for partition trees — never on the output file merely existing, because HyPhy
rewrites its JSON at every breakpoint stage and IQ-TREE writes `.treefile` before
bootstrapping ends. `step21` and `step22` must be run after `step18` so the
positive-selection filter is applied rather than recorded as `NOT_YET_APPLIED`; both scripts
say so in their own logs, and that re-run has been done. `step25` and `step26` are the
verification-pass additions: `step25` rebuilds the QC record over the CDS actually carried
forward, and `step26` needs `step25`'s output plus `step18`'s.

The same dependency order is declared machine-readably in the [`Snakefile`](Snakefile)
(`snakemake -n` to see the plan). Stages 1-7 appear there as frozen inputs, not rules:
re-downloading from NCBI would change the cohort, which is pinned by accession in
`19_supplements/S2_final_accessions.tsv`. The rule graph and script paths have been checked,
but the workflow has not been dry-run, because snakemake is not installed in the environment
the pipeline was executed in.

### Requirements

* **Python 3.11** with `biopython`, `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn`,
  `pyfamsa` (codon-aware MSA), `mhcflurry` (class I prediction, models fetched via
  `mhcflurry-downloads fetch models_class1_presentation`)
* **IQ-TREE 2.4.0** — bundled in `hsv_curated_dataset/bin/`
* **HyPhy 2.5.103 (MP)** — via Docker (`hyphy/hyphy`); **Docker Desktop must be running**.
  The version and image id are recorded to `14_selection/gard_provenance.json` on every run.
* **BepiPred-3.0** (linear B-cell, ESM-2 650M + 5-fold FFNN ensemble) and **DiscoTope-3.0**
  (conformational B-cell, ESM-IF1 + 100-model XGBoost ensemble) — run locally on GPU in
  their own virtual environments. Both are free for academic / non-commercial use
  (BepiPred: DTU Health Tech terms; DiscoTope: CC BY-NC 4.0). An earlier version of this
  README listed both as licence-gated; that was wrong — only NetMHCIIpan is. Raw output is
  archived in `13_epitopes/bcell_raw/`, versions in
  `13_epitopes/bcell_predictor_provenance.tsv`.
* **NetMHCIIpan-4.3** (CD4 class II) — licence-gated, with no public source repository, so
  the predictions were produced on the DTU web server in two submissions over the 27-allele
  IEDB class II reference set and adapted by `step29`. Raw output archived gzipped in
  `13_epitopes/classii_raw/`; 610,848 predictions, all coordinate-verified. No class II
  predictor is benchmarked here, so no CD4 epitope is claimed.

---

## Methods notes

* **Aligner.** The proposal specified MAFFT or MACSE. MAFFT's distribution site was unreachable
  and neither tool has a Windows build available here, so alignment uses **FAMSA**
  (Deorowicz, Debudaj-Grabysz & Gudyś, *Sci Rep* 2016) at the protein level with codon
  back-projection — the same MACSE-style strategy, with a genuine progressive profile.
* **Pooled cohorts.** Combined HSV-1 + HSV-2 alignments are built and used **only** for
  between-type conservation and cross-protective epitope identification. They are never used for
  dN/dS or recombination inference.
* **Significance.** FEL and MEME site p-values receive Benjamini–Hochberg FDR within each
  analysis unit; BUSTED receives BH across the gene × cohort family; aBSREL uses HyPhy's
  Holm–Bonferroni corrected branch p-values; FUBAR is reported at posterior ≥ 0.95. A site is
  called positively selected only with q<0.05 support from ≥2 methods.
* **Power.** Post-deduplication cohorts range from 29 to 216 unique sequences; HSV-2 UL1 (n=32),
  US6 (n=54) and UL53 (n=74) are underpowered for site-level selection tests. Absence of a
  detected breakpoint or selected site is reported as such, never as evidence of clonality.
* **Negatives.** Untested regions are never treated as negative epitope labels.

## Citation

HyPhy methods per Kosakovsky Pond et al.; IQ-TREE per Minh et al.; FAMSA per Deorowicz et al.;
MHCflurry per O'Donnell et al.; epitope evidence from IEDB.
