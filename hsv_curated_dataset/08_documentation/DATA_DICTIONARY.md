# Data Dictionary

Every directory and every named file in this repository: what it contains, which
script produced it, and which script consumes it.

Written 2026-11-21, after the repository cleanup of the same date. Companion to
`methods_draft.md` (what was done and why) and `../../CORRECTIONS_AND_PENDING_WORK.md`
(what was withdrawn and what remains open).

**How to read the tables.** *Producer* is the script that writes the file;
*Consumer* is the script that reads it. A file with no consumer is a terminal
deliverable (a supplement, a figure, a QC record) and that is noted explicitly, not
left blank. Paths are relative to `hsv_curated_dataset/` unless stated otherwise.

---

## 0. Repository root

| File | Contents |
|---|---|
| `README.md` | Project overview, cohort summary, stage-by-stage pipeline map, how to re-run. |
| `CORRECTIONS_AND_PENDING_WORK.md` | The error register. Every finding that was withdrawn or corrected during the project, with the reason. Read this before trusting any number quoted elsewhere. |
| `Snakefile` | 25 rules reproducing the pipeline from stage 13 onward. Earlier stages are frozen (see §"Frozen stages" below). |
| `requirements.txt` | Python dependencies with pinned versions. |
| `.gitignore` | Excludes tool binaries, bytecode, HyPhy scratch, and the two bulk genome directories (see §1). |
| `.gitattributes` | Pins line-ending handling per file type. Added 2026-11-21 after `population_genotype_map.p` was found corrupted in the working tree by LF→CRLF conversion; see §"Data integrity" below. |

---

## 1. Frozen stages (00–08): cohort assembly

These stages ran once, early, and their outputs are inputs to everything after
them. They are **frozen**: the scripts are retained as the provenance record, but
re-running them against today's NCBI would return a different cohort. The
accession list in `19_supplements/S2_final_accessions.tsv` is what pins the
analysed set.

### `00_raw_downloads/` — 222 MB, **not tracked in git**

Raw NCBI `datasets` output per virus type: the download `.zip`, the unpacked
`ncbi_dataset/` tree, and the concatenated `.fasta`.

Excluded from git because two files exceed GitHub's 100 MB hard limit and the
content is exactly reproducible: `S2_final_accessions.tsv` pins every accession,
and `step03_download_full_cohort.py` re-fetches them. Everything *derived* from
these files is tracked, so the analysis is checkable without re-downloading.

| File | Producer | Consumer |
|---|---|---|
| `HSV{1,2}/*.zip` | `step03_download_full_cohort.py` | `step04_genome_filtering.py` |
| `HSV{1,2}/*.fasta` | `step03_download_full_cohort.py` | `step04_genome_filtering.py` |
| `unpacked/md5sum.txt` | NCBI `datasets` CLI | checksum verification |
| `unpacked/ncbi_dataset/data/genomic.fna` | NCBI `datasets` CLI | `step04_genome_filtering.py` |
| `unpacked/ncbi_dataset/data/data_report.jsonl` | NCBI `datasets` CLI | `step03_query_complete_genomes.py` |
| `unpacked/ncbi_dataset/data/dataset_catalog.json` | NCBI `datasets` CLI | manifest, unused downstream |

### `01_reference_genomes/` — 1.4 MB

| File | Contents |
|---|---|
| `NC_001806.2.fasta`, `NC_001806.2.gbk` | HSV-1 reference genome and GenBank annotation (strain 17). |
| `NC_001798.2.fasta`, `NC_001798.2.gbk` | HSV-2 reference genome and GenBank annotation (strain HG52). |
| `JQ673480.1_KOS.fasta`, `JX142173.1_McKrae.fasta`, `GU734771.1_F.fasta` | Three additional HSV-1 lab strains, used for the lab-strain audit. |
| `reference_cds_table.tsv` | Every CDS feature parsed out of the two `.gbk` files: gene name, coordinates, strand, product, protein length. **The coordinate authority for the whole project.** |

Producer: `step02_download_reference_genomes.py`, with
`step04_5a_fix_reference.py` repairing one malformed record.
Consumer: `step06_cds_extraction.py`, and every downstream stage that maps a
residue back to a reference position.

### `02_filtered_genomes/` — 102 MB, **not tracked in git** (same reason as `00_`)

The genome sets surviving each filtering tier. Producer:
`step04_genome_filtering.py` (with `step04_5_verification_and_recovery.py`
recovering borderline records). Consumer: `step06_cds_extraction.py`.

| File | Contents |
|---|---|
| `HSV1_filtered.fasta` | HSV-1 genomes passing all filters. **The HSV-1 analysis set.** |
| `HSV2_filtered.fasta` | HSV-2 genomes passing all filters. **The HSV-2 analysis set.** |
| `HSV2_filtered_tier2.fasta` | HSV-2 genomes passing only the relaxed tier-2 criteria, audited separately in `07_qc_reports/04f_tier2_qa_report.tsv`. |
| `HSV1_flagged.fasta`, `HSV2_flagged.fasta` | Genomes flagged and excluded, retained so every exclusion can be re-examined. Reasons are in `07_qc_reports/04_genome_filtering_log.tsv`. |

### `03_metadata/` — 1.4 MB

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `HSV1_raw_summary.jsonl`, `HSV2_raw_summary.jsonl` | Raw NCBI metadata, one JSON object per genome, unmodified. | `step03_query_complete_genomes.py` | `step05_metadata_and_verifications.py` |
| `raw_metadata_all_genomes.tsv` | The same records flattened to a table, before curation. | `step05_metadata_and_verifications.py` | `step05` (self) |
| `curated_metadata.tsv` | **The metadata table used throughout.** Per accession: country, collection year, host, isolation source, clinical context, with normalised country and year fields. | `step05_metadata_and_verifications.py` | `step20`, `step21`, `step27` |
| `patient_assignments.tsv` | Inferred patient groupings, for detecting multiple isolates from one individual. | `step20_patient_and_sampling_bias.py` | `step20` (self), `step21` |
| `sampling_bias_tables.tsv` | Geographic and temporal skew of the cohort, quantified. | `step20_patient_and_sampling_bias.py` | terminal — reported in `results_draft.md` §3 |
| `patient_bias.log` | Run log for `step20`. | `step20` | terminal |

### `04_cds_unaligned/` — 72 MB, 148 files

Per-gene nucleotide FASTA, one file per gene per virus type, extracted from the
filtered genomes by GenBank coordinate projection. **Unfiltered** — includes
sequences that later fail QC.
Producer: `step06_cds_extraction.py`. Consumer: `step07_part1_clean_qc.py`.

### `05_cds_qc_passed/` — 28 MB, 82 files

The same per-gene FASTA after coding-sequence QC: length within tolerance of the
reference, no internal stop codons, clean translation, length a multiple of three.
41 genes per virus type survive, from 74 attempted.
Producer: `step07_part1_clean_qc.py`, with `step07_gaps_abc.py` handling the
gap-tolerance cases. Consumer: `step08_dedup_and_step09_align.py`.

### `06_alignments_v2/` — 198 MB, 713 files

**The live alignment set.** The v1 alignments were withdrawn (out-of-register
codons from pooling the two virus types before aligning) and have been deleted;
they remain in git history at commit `4f56731` and earlier. Everything here is
per-virus-type, built from the QC-passed CDS.

| Path | Contents | Producer | Consumer |
|---|---|---|---|
| `HSV1/`, `HSV2/` | Codon-aware per-gene MSA, 20 genes each. | `step09b_rebuild_codon_msa.py` | `step09c`, `step13` |
| `combined/` | HSV-1 and HSV-2 aligned together per gene, for cross-type conservation only — never for selection. | `step09b_rebuild_codon_msa.py` | `step13_conservation_layer.py` |
| `masked/` | The same alignments with low-confidence columns masked. | `step09c_mask_and_coordinate_map.py` | `step10b`, `step13` |
| `coordinate_maps/HSV{1,2}/<GENE>_colmap.tsv` | Per gene, the mapping from alignment column → reference codon position. **The register authority**: every residue-level annotation downstream resolves through these. One file per gene per type, 20 each. | `step09c_mask_and_coordinate_map.py` | `step13`, `step15`, `step17`, `step21`, `step22`, `step28`, `step29` |
| `hyphy_ready/` | Alignments in the form HyPhy consumes, plus `*.best-gard` (GARD breakpoint results), `*.fit.bf` (fitted model), `*_idmap.tsv` (sequence-name mapping, because HyPhy mangles long names). | `step10b_prep_and_trees.py`, `step16_gard_per_type.py` | `step17_selection_per_type.py` |
| `trees_v2/HSV{1,2}/` | IQ-TREE output per gene, 9 files each. **`.treefile` is the ML tree and the only one read downstream.** The rest are IQ-TREE's own auxiliary output, retained unmodified as the run record: `.contree` (consensus tree), `.bionj` (distance-based starting tree), `.iqtree` (full human-readable report), `.log` (run log), `.ckp.gz` (checkpoint), `.mldist` (ML distance matrix), `.model.gz` (model-selection results), `.splits.nex` (bootstrap split support, NEXUS). | `step10b_prep_and_trees.py` | `step16`, `step17` |
| `trees_v2/tree_summary_v2.tsv` | Per gene per type: chosen substitution model, log-likelihood, tree length, bootstrap support summary. | `step10b_prep_and_trees.py` | terminal — QC record |
| `alignment_qc_v2.tsv` | Per-gene alignment statistics: sequences, columns, gap fraction, parsimony-informative sites. | `step09b` | terminal — QC record |
| `masking_report_v2.tsv` | Per gene, how many columns were masked and why. | `step09c` | terminal — QC record |
| `outlier_sequences_v2.tsv` | Sequences flagged as alignment outliers, with the metric that flagged them. | `step09b` | terminal — QC record |
| `alignment_rebuild.log`, `masking.log`, `trees_v2/trees.log`, `trees_v2/stdout_trees.log` | Run logs; `stdout_trees.log` is IQ-TREE's raw stdout, captured verbatim. | respective scripts | terminal |

---

## 2. QC records

### `07_qc_reports/` — 10.4 MB

Numbered to match the stage that produced each. These are the audit trail: every
sequence that was dropped, and the reason. **No downstream script reads them**
except `step25_qc_reconciliation.py`, which cross-checks them against each other,
and `step24_build_supplements.py`, which rolls them into S1/S3.

| File | Contents |
|---|---|
| `04_genome_filtering_log.tsv` | Per genome: every filter applied, pass/fail, reason. |
| `04b_hsv2_exclusion_analysis.tsv` | Why each excluded HSV-2 genome was excluded. |
| `04c_lab_strain_audit.tsv` | Lab-passaged strains identified and their disposition. |
| `04e_species_verification.tsv`, `04e_species_disagreement_full_list.tsv` | Species assignment checked against sequence identity; the second lists every case where the NCBI label and the sequence disagree. |
| `04f_tier2_qa_report.tsv` | Second-tier QA on borderline genomes. |
| `05_metadata_completeness.tsv` | Per metadata field, how many records have it. |
| `06_cds_extraction_log.tsv` | Per gene per genome: extraction outcome. |
| `06b_translation_validation_per_gene.tsv` | Translation check results per gene. |
| `06c_per_gene_extraction_success.tsv` | Extraction success rate per gene. |
| `06d_rl2_rs1_nucleotide_diagnostic.tsv` | Diagnostic for RL2 and RS1, two genes with repeat-driven extraction problems (producer: `step06d_rl2_rs1_diagnostic.py`). |
| `07_cds_qc_log.tsv`, `07d_cds_qc_log_complete.tsv` | Per-sequence CDS QC outcomes; the `_complete` version is the authoritative one. |
| `07b_us7_failure_breakdown.tsv` | Why US7 fails QC disproportionately. |
| `07c_length_filter_audit.tsv` | Effect of the length filter, per gene. |
| `07e_qc_reconciliation.tsv` | Cross-check that the counts in all the above agree (producer: `step25_qc_reconciliation.py`). |
| `08_dedup_log.tsv` | Per sequence: whether it was removed as a duplicate, and of what. **A frozen input** — `step08`'s dedup half cannot be re-run independently. |
| `08b_dedup_cluster_sanity.tsv` | Sanity check on the dedup clusters. |
| `09_conservation_metric_diagnostic.tsv`, `09b_alignment_projection_validation.tsv` | Checks that conservation metrics and coordinate projection behave as intended. |
| `11_prereq2_trees_summary.tsv`, `11_prereq3_cleaned_masked_summary.tsv`, `11b_tier_b_masking_investigation.tsv`, `11c_tier_c_lc_annotation.tsv` | **Withdrawn-v1 records.** Small summaries from the superseded selection pipeline, retained because `CORRECTIONS_AND_PENDING_WORK.md` cites them. |
| `qc_reconciliation.log` | Run log. |

### `07_qc_reports/hyphy/` — 1.3 GB, withdrawn

Output of the **superseded v1 selection pipeline**. Nothing live reads it. It is
retained only as the record behind the withdrawn results; the current selection
results are in `14_selection/`.

| Path | Contents |
|---|---|
| `stage1/tier_{a,b,c}/` | 234 raw HyPhy JSON, 1.1 GB. Full model fits from the withdrawn stage-1 runs, in the three tiers the v1 pipeline used. |
| `stage2/` | 62 raw HyPhy JSON, 124 MB: `<GENE>_<TYPE>_masked_{FEL,MEME}.json`, where `<TYPE>` includes `combined` — the pooled-type runs that are exactly what the withdrawal was about. |
| `gard/` | 10 v1 GARD JSON, `<GENE>_combined_GARD.json`. Pooled-type GARD, superseded by the per-type GARD in `14_selection/gard/`. |
| `stage1_summary.tsv`, `stage2_summary.tsv`, `gard_summary.tsv` | Per-gene result summaries from the three withdrawn run stages. |
| `stage1_progress.log`, `stage2_progress.log`, `gard_progress.log`, `tier_c_masked_progress.log` | Run progress logs for the withdrawn stages. |
| `stage2_trigger_genes.tsv`, `stage2_trigger_genes_tightened.tsv`, `stage2_trigger_genes_final.tsv` | The genes selected to advance to v1 stage 2, under three successive criteria. The three versions are kept because the criteria changed mid-analysis, which is itself part of what was withdrawn. |
| `stage1_biological_validation.tsv`, `stage1_outlier_investigation.tsv`, `validation_gB_codons_59_61.tsv`, `retrospective_fdr_on_old_stage2.tsv` | The v1 validation and retrospective-FDR checks. |

Those 14 small tables and logs are **the part of the withdrawn arm worth
keeping** — they are what the corrections register cites.

---

## 3. Scripts

### `09_scripts/` — 37 live scripts

Numbered by stage. Four are not stage scripts:

| File | Role |
|---|---|
| `check_dependencies.py` | Verifies every external tool and Python package is present at the expected version. Run first. |
| `install_prerequisites.py` | Installs the external tools (IQ-TREE, HyPhy, MAFFT). |
| `hyphy_runner.py` | Shared HyPhy subprocess wrapper: builds command lines, enforces timeouts, parses JSON, caches completed runs via `.done` markers. **Imported by `step16` and `step17`** — not standalone. |
| `step30_provenance_audit.py` | Audits this dictionary against the files on disk: fails, naming the gaps, if any data file has no entry here. Read-only. Run as `snakemake provenance_audit` or directly. **Numbered as a stage script but produces no data** — it is a check, and it must be able to fail without blocking a pipeline run, which is why it is not in `rule all`. |

Stage scripts `step02`–`step29` correspond one-to-one with the numbered output
directories; each is named for what it produces, and the Producer columns
throughout this document are the authoritative mapping.

One stage script is **partially withdrawn**:
`step08_dedup_and_step09_align.py`. Its step-8 half (deduplication) produced
`07_qc_reports/08_dedup_log.tsv`, a frozen input that cannot be regenerated
independently of this script — which is the reason the script is retained. Its
step-9 half wrote the withdrawn `06_alignments/` and must not be run; the
script's own header says so.

### `09_scripts/deprecated_v1/` — 10 scripts

The withdrawn v1 selection pipeline: `step10_and_diagnostics_c1_c2_c3.py`,
`step11_*` (6 files), `step11b`, `step11d`, `step12_*` (2 files). Retained as
the evidence underpinning the corrections register — **do not run them.** Their
outputs are the withdrawn `07_qc_reports/hyphy/` tree.

### `08_documentation/`

| File | Contents |
|---|---|
| `methods_draft.md` | Manuscript methods, §1–§19. The authority on what was done and why. |
| `results_draft.md` | Manuscript results. |
| `DATA_DICTIONARY.md` | This file. |
| `MANUSCRIPT_CURATION.md` | Pre-submission assessment, 2026-11-21: which drafted claims reproduce from the result tables, the issues to fix before submission, and the figure plan. |
| `README.md` | Documentation index. |

### `bin/` — 75 MB, **not tracked in git**

External tool binaries (IQ-TREE 2.4.0 and its support files). Reinstalled by
`install_prerequisites.py`; excluded from git as platform-specific binaries.

---

## 4. Analysis stages (11–20)

### `11_conservation/` — 1.3 MB

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `per_residue_conservation.tsv` | Per gene per reference residue: Shannon entropy, identity fraction, gap fraction, within-type and cross-type conservation. **The residue-level conservation authority.** | `step13_conservation_layer.py` | `step21`, `step22`, `step27` |
| `gene_level_conservation.tsv` | The same metrics aggregated per gene. | `step13` | `step21`, `step27` |
| `conserved_windows.tsv` | Contiguous runs of conserved residues, the candidate regions for epitope selection. | `step13` | `step21` |
| `conservation.log` | Run log. | `step13` | terminal |

### `12_iedb/` — 21 MB

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `iedb_raw_hsv1.jsonl`, `iedb_raw_hsv2.jsonl` | Unmodified IEDB API responses. Kept so curation decisions can be re-derived. | `step14_iedb_curation.py` | `step14` (self) |
| `iedb_curated_epitopes.tsv` | Experimentally validated HSV epitopes passing curation: sequence, host, assay, MHC restriction, mapped reference coordinates. **The benchmark truth set.** | `step14` | `step19`, `step21` |
| `iedb_excluded_epitopes.tsv` | Epitopes excluded, each with the exclusion reason. | `step14` | terminal — audit |
| `iedb_residue_flags.tsv` | Per residue: whether it falls inside a validated epitope. | `step14` | `step21`, `step22` |
| `iedb_curation_summary.tsv` | Counts per curation decision. | `step14` | `step24` |
| `iedb.log` | Run log. | `step14` | terminal |

### `13_epitopes/` — 48 MB

Three prediction arms: CD8 (class I), CD4 (class II), and B-cell.

**A standing caveat applies to all three.** None of the three predictors'
published thresholds behaves as documented on this data — BepiPred flags 37.3% of
residues, DiscoTope flags none at its high-confidence cut, and
`classii_best_percentile_rank` falls below 1% for 76% of candidates because it is
a minimum over ~100 alleles. All three are therefore carried as **continuous
annotations, never as filters.** See `methods_draft.md` §16–§17.

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `class1_peptide_predictions.tsv` | Per 8–11-mer peptide per HLA class I allele: predicted affinity and percentile rank. | `step15_tcell_class1_prediction.py` | `step19`, `step21` |
| `class1_binders_by_allele.tsv` | Binders per allele at the strong/weak thresholds. | `step15` | `step23` |
| `class1_prediction_summary.tsv` | Per-gene prediction counts. | `step15` | `step24` |
| `classii_predictions.tsv.gz` | **610,848 rows.** Per 15-mer per HLA class II allele: `Score_EL`, `Rank_EL`, `Score_BA`, `nM` affinity, `Rank_BA`, binding `Core`, and `inverted` (1 where NetMHCIIpan models the peptide back-to-front in the groove — 13,362 rows, 4.9%, exclusively HLA-DP). Gzipped. | `step29_classii_predictions.py` | `step21`, `step23` |
| `classii_binders_by_allele.tsv` | Class II binders per allele. | `step29` | `step23` |
| `bcell_linear_scores.tsv` | Per residue: BepiPred-3.0 score (ESM-2 650M + 5-fold FFNN ensemble), plus the rolling-window-9 smoothed value. | `step28_bcell_predictors.py` | `step21`, `step22` |
| `bcell_conformational_scores.tsv` | Per residue of each solved structure: DiscoTope-3.0 score (ESM-IF1 + 100-model XGBoost ensemble). | `step28` | `step21`, `step22` |
| `bcell_predictor_provenance.tsv` | Exact model versions, weights, thresholds and command lines for both B-cell predictors. | `step28` | `step24` |
| `bcell_inputs/` | The antigen FASTA submitted to BepiPred-3.0. | `step28` | external predictor |
| `classii_inputs/` | The four submission FASTA — `classii_input_ALL40.fasta`, `classii_input_PILOT_gD.fasta`, `classii_input_HALF_A.fasta`, `classii_input_HALF_B.fasta` — and the two allele lists, `alleles_RUN1_DR_15.txt` (15 DR alleles) and `alleles_RUN2_DQ_DP_12.txt` (12 DQ/DP alleles). Split into two runs because the DTU server caps a submission at 15 alleles; the 27-allele IEDB reference set does not fit in one. | `step29` | DTU web server, run manually |
| `bcell_raw/bepipred3/` | Unmodified BepiPred-3.0 output: `raw_output.csv` (per-residue scores) and the three returned FASTA. | external predictor | `step28` |
| `bcell_raw/discotope3/` | Unmodified DiscoTope-3.0 output, one CSV per structure chain: `<PDB>_<CHAIN>_discotope3.csv`, 7 files across 3 structures. Per chain because DiscoTope was run in single-chain mode. | external predictor | `step28` |
| `classii_raw/` | Unmodified NetMHCIIpan-4.3 server output, gzipped: `run1_DR_15alleles.xls.gz`, `run2_DQ_DP_12alleles.xls.gz`. NetMHCIIpan is licence-gated with no public source, so these raw files **are** the provenance — the version string is parsed from the `#` command line inside them. | DTU web server, run manually | `step29`, `step24` |
| `classii_raw/pilot/` | `905068_NetMHCIIpan.xls.gz`, the single-antigen pilot submission (gD only) used to confirm the output format and register before committing to the full run, plus `README.md` recording why it exists. | DTU web server, run manually | terminal — provenance record |
| `class1.log`, `classii.log`, `bcell.log` | Run logs for `step15`, `step29`, `step28`. | respective scripts | terminal |
| `stdout_epitopes.log` | Raw stdout from the B-cell and class II adapter runs. Was sitting at the dataset root as `13_epitopes_stdout.log`; moved here in the 2026-11-21 cleanup. | `step28`, `step29` | terminal |

### `14_selection/` — 587 MB

**The live selection results.** Per virus type, never pooled.

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `gard/HSV{1,2}/` | GARD recombination-breakpoint JSON per gene, plus `.done` cache markers. | `step16_gard_per_type.py` | `step17` |
| `gard_summary_per_type.tsv` | Breakpoints detected per gene per type. | `step16` | `step17`, `step27` |
| `gard_provenance.json` | GARD version, command lines, run timestamps. | `step16` | `step24` |
| `partitions/HSV{1,2}/` | Per-partition alignments and trees where GARD found breakpoints: `.fasta`, `.bionj`, `.treefile.gz`, `.alnmd5` (alignment checksum), `.cache`. | `step16` | `step17`, `step26` |
| `hyphy/HSV{1,2}/` | **610 raw HyPhy JSON** (FEL, MEME, FUBAR, BUSTED, aBSREL) plus `.done` markers. The primary selection evidence. | `step17_selection_per_type.py` | `step18`, `step26` |
| `gene_level_results.tsv` | Per gene per type: gene-wide selection tests, ω estimates, p-values. | `step17` | `step18`, `step21`, `step27` |
| `site_level_results.tsv` | Per site per type: site-wise selection statistics from each method. | `step17` | `step18`, `step21`, `step26` |
| `branch_level_results.tsv` | Per branch: episodic-selection results. | `step17` | `step18`, `step27` |
| `selection_runs.tsv` | One row per HyPhy invocation: method, gene, type, exit status, runtime. The completeness ledger. | `step17` | `step24`, `step25` |
| `selected_sites_final.tsv` | **The FDR-corrected selected-site set.** Sites surviving multiple-testing correction, with the methods supporting each. | `step18_fdr_and_concordance.py` | `step21`, `step22`, `step26`, `step27` |
| `selection_concordance.tsv` | Agreement between methods, site by site. | `step18` | `step27` |
| `site_robustness.tsv` | Whether each selected site survives partition-boundary and model perturbation. | `step26_site_robustness.py` | `step21` |
| `selection.log`, `gard.log`, `fdr.log`, `site_robustness.log` | Run logs for `step17`, `step16`, `step18`, `step26`. | respective scripts | terminal |
| `stdout_selection.log`, `stdout_selection2.log`, `stdout_selection_final.log`, `stdout_gard.log` | Raw HyPhy stdout, captured verbatim across the three selection passes and the GARD pass. Retained because HyPhy reports convergence warnings on stdout that do not appear in the JSON. | `step16`, `step17` | terminal |

### `15_benchmark/` — 8.2 MB

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `labelled_peptides.tsv` | Peptides labelled positive/negative against the IEDB truth set, the benchmark input. | `step19_benchmark_predictors.py` | `step19` (self) |
| `benchmark_per_gene.tsv` | Predictor performance per gene. | `step19` | `step27` |
| `benchmark_summary.tsv` | Overall predictor performance: AUC, precision/recall at the published thresholds. **This is where the threshold problem is quantified.** | `step19` | `step21`, `step27` |
| `benchmark.log` | Run log. | `step19` | terminal |

### `16_shortlist/` — 6.5 MB

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `epitope_candidates_ranked.tsv` | **The primary integrated output.** Every candidate with all annotation layers joined: conservation, selection, class I/II binding, B-cell scores, IEDB validation, structure, and the composite rank. | `step21_integrate_and_rank.py` | `step22`, `step23`, `step24`, `step27` |
| `shortlist_primary.tsv` | The top-ranked candidates. | `step21` | `step24`, `step27` |
| `shortlist_cross_protective.tsv` | Candidates conserved across both virus types. | `step21` | `step24` |
| `shortlist_by_gene.tsv` | Best candidates per gene. | `step21` | `step24` |
| `integration.log` | Run log, including the diagnostic explaining the class II rank distribution. | `step21` | terminal |

### `17_structures/` — 0.26 MB + 2.7 MB PDB

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `pdb/*.pdb` | Three solved HSV glycoprotein structures, as downloaded from the PDB. | `step22_structure_mapping.py` | `step22`, `step28` |
| `chain_assignments.tsv` | Which PDB chain corresponds to which gene. **Needed because chains within one structure are different proteins** — e.g. 3M1C chain A is gH/UL22 and chain B is gL/UL1; treating them as replicates produces a meaningless correlation. | `step22` | `step28` |
| `structure_residue_annotation.tsv` | Per structural residue: solvent accessibility, secondary structure, mapped reference position, and the conservation/selection annotations joined on. | `step22` | `step21`, `step27` |
| `conserved_exposed_patches.tsv` | Surface patches that are both conserved and solvent-exposed. | `step22` | terminal — reported |
| `figures/*.pml` | PyMOL scripts regenerating the structure figures. | `step22` | terminal — run manually |
| `structures.log` | Run log. | `step22` | terminal |

### `18_population_coverage/` — 3.6 MB

Coverage computed by the Bui et al. (2006) method: per locus
`1 − (1 − Σf)²`, combined `1 − Π(1 − cov_locus)`.

| File | Contents | Producer | Consumer |
|---|---|---|---|
| `population_genotype_map.p` | Pickled IEDB allele-frequency reference data, 274 populations. | `step23_population_coverage.py` | `step23` (self) |
| `allele_frequency_provenance.json` | Source and version of the frequency data. | `step23` | `step24` |
| `coverage_by_population.tsv` | Class I coverage per population. | `step23` | `step27` |
| `coverage_summary.tsv` | Class I coverage summary. | `step23` | `step24`, `step27` |
| `coverage_by_population_classii.tsv` | Class II coverage per population, **bracketed**: an upper bound treating loci as independent, and a conservative bound computing `cov(DP) = cov(DPA1) × cov(DPB1)`. Carries `n_loci_typed`, `loci_typed`, and `scorable`. | `step23` | `step27` |
| `coverage_summary_classii.tsv` | Class II summary, over `scorable` populations only — the 64 of 274 typed at all five presenting loci. **Restricting to these matters: scoring untyped loci as zero coverage produced a 16-point artifact that read as a geographic finding.** | `step23` | `step24`, `step27` |
| `coverage_combined_classI_classII.tsv` | Combined class I + class II coverage. | `step23` | `step27` |
| `minimal_panel_greedy.tsv` | Greedy-selected minimal allele panel reaching target coverage, **low-complexity peptides barred from the search**. Columns: `alleles_bound` (alleles the peptide binds) and `new_alleles` (of those, not already covered); `new_alleles` cumulatively sums to `cumulative_alleles`. | `step23` | terminal — reported |
| `minimal_panel_greedy_unconstrained.tsv` | The same search without the low-complexity bar, kept so the cost of the exclusion is checkable rather than asserted. Not the panel the manuscript reports. | `step23` | terminal — comparison record |
| `population_coverage.log` | Run log. | `step23` | terminal |

### `19_supplements/` — 0.26 MB

Manuscript supplementary tables. All terminal.

| File | Contents |
|---|---|
| `S1_data_flow.tsv` | Sequence counts at every pipeline stage — the CONSORT-style flow table. |
| `S2_final_accessions.tsv` | **Every accession in the final analysis.** The reproducibility anchor for the untracked genome directories. |
| `S3_removed_sequences.tsv` | Every sequence removed, with the stage and reason. |
| `S4_per_gene_cohorts.tsv` | Per gene: how many sequences of each type entered the alignment. |
| `S5_software_versions.tsv` | Every tool and model version, including the predictor versions parsed from raw output. No "NOT RUN" rows remain. |
| `supplements.log` | Run log. |

Producer for all: `step24_build_supplements.py`.

### `20_figures/` — 1.2 MB

Six manuscript figures, each as `.pdf` and `.png`, plus `figures.log` (run log).
All terminal. Producer: `step27_figures.py`.

| Figure | Contents |
|---|---|
| `fig1_cohort_size` | Cohort composition and sequence counts. |
| `fig2_conservation` | Conservation across genes and residues. |
| `fig3_selection_gene_branch` | Gene- and branch-level selection results. |
| `fig4_site_variant_support` | Per-site selection support across methods. |
| `fig5_predictor_benchmark` | Predictor performance against the IEDB truth set. |
| `fig6_population_coverage` | Population coverage, class I and class II. |

---

## Frozen stages

Stages 00–08 are frozen: re-running them against today's NCBI would return a
different cohort, and `08_dedup_log.tsv` in particular cannot be regenerated
independently of `step08`'s alignment half. The `Snakefile` therefore reproduces
the pipeline from stage 13 onward and treats stages 00–08 as given inputs. The
stage 00–08 scripts are retained as the provenance record of how the frozen
inputs were made.

## Withdrawn material

Deleted in the 2026-11-21 cleanup, recoverable from git history at commit
`4f56731`:

- `06_alignments/` — v1 alignments, out of codon register because the two virus
  types were pooled before aligning. Superseded by `06_alignments_v2/`.
- `_archived_discarded/` — superseded audit scratch.

Deleted and **not** recoverable (never tracked):

- `07_qc_reports/hyphy/_tmp/` — 410 MB of HyPhy scratch and cache dumps.
- Two stray `test_UL27_HSV1_fubar*.json` probes.

Retained but withdrawn, read by nothing:

- `07_qc_reports/hyphy/stage1/`, `stage2/` — 1.3 GB of v1 HyPhy JSON.
- `09_scripts/deprecated_v1/` — the 10 v1 selection scripts.

## Data integrity

**One file was found corrupted on disk, 2026-11-21.**
`18_population_coverage/population_genotype_map.p` is the pickled IEDB
allele-frequency reference data. Git's binary heuristic did not detect it as
binary, so LF was converted to CRLF on checkout and 380,776 CRLF pairs were
written into it: 4,000,753 bytes on disk against 3,619,977 committed, a
difference of exactly the number of CRLF pairs. The pickle became unloadable —
its `GLOBAL` opcode read `collections\r` / `OrderedDict\r` — and `step23`
aborted on it.

The committed blob was clean, so the file was restored from `git show HEAD:…`
and verified to unpickle (2 top-level keys, `I` and `II`). All 136 other binary
files in the repository (`.gz`, `.png`, `.pdf`, `.zip`) were checked by
decompressing or reading their magic bytes and are intact; those formats carry
NUL bytes early enough for git to classify them correctly.

`.gitattributes` now pins `*.p`, `*.pkl`, `*.npy` and the other binary formats
so this cannot recur. If a pickle or archive in this repository ever fails to
load, check its size against `git cat-file -s HEAD:<path>` before assuming the
data is lost.
