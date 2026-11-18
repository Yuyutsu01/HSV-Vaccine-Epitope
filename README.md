# Herpes Simplex Virus (HSV-1 & HSV-2) Evolutionary Selection & Vaccine Epitope Profiling

[![Status](https://img.shields.io/badge/Status-100%25_Completed-success.svg)](#)
[![Pipeline](https://img.shields.io/badge/HyPhy-2.5.103_MP-blue.svg)](#)
[![Targets](https://img.shields.io/badge/Total_Runs-342%2F342_Complete-brightgreen.svg)](#)

A comprehensive, phylogenetically curated evolutionary selection and recombination detection pipeline evaluating clinical isolates of **Herpes Simplex Virus Type 1 (HSV-1)** and **Type 2 (HSV-2)** across all structural glycoproteins, replication machinery, and regulatory tegument targets.

---

## 🔬 Scientific Highlights & Vaccine Design Insights

1. **Strictly Clonal, Low-Plasticity Subunit Targets:**
   * **Glycoprotein D (`US6`):** Zero recombination breakpoints ($\Delta\text{AICc} = 0.00$) and high structural conservation (only codon 53 under episodic selection in HSV-1). Confirmed as an optimal monophyletic subunit candidate resistant to chimeric evasion.
   * **Glycoprotein L (`UL1`):** Zero recombination breakpoints ($\Delta\text{AICc} = 0.00$) with tight chaperone conservation, functioning as a stable co-immunogen alongside gH.
   * **Glycoproteins C (`UL44`), M (`UL10`), and K (`UL53`):** All exhibit strictly clonal evolutionary trajectories without mosaic chimeric breakpoints.

2. **Recombinant Modular Antigens:**
   * **Glycoprotein B (`UL27`):** Recombination breakpoint identified at **nucleotide position 735** ($\Delta\text{AICc} = 310.84$), physically delineating Domain I from Domain II ectodomains. Subunit immunogen design must accommodate domain-specific phylogenetic variation.
   * **Glycoprotein H (`UL22`):** Significant breakpoint at **nucleotide position 604** ($\Delta\text{AICc} = 468.23$), separating the N-terminal gL-binding chaperone region from the core fusion apparatus.
   * **DNA Polymerase (`UL30`):** Strong recombination breakpoint at **nucleotide position 2018** / codon 673 ($\Delta\text{AICc} = 6,619.11$) directly bisecting the exonuclease domain and the Region II acyclovir drug-resistance pocket.
   * **VP11/12 Tegument Protein (`UL46`):** Recombination breakpoint at **nucleotide position 1713** ($\Delta\text{AICc} = 2,385.71$), decoupling the N-terminal kinase-activating domain from the C-terminal structural tail.

---

## 📊 Pipeline Completion Architecture (342 / 342 Runs Complete)

```
========================================================================================
1. STAGE 1 — SELECTION PRE-SCREENING (270 / 270 Runs - 100% COMPLETE)
========================================================================================
  • Tier A (Standard Screening across 33 targets × 3 methods):          99 / 99  (100%)
  • Tier B (Expanded Cohorts across 33 targets × 3 methods):            99 / 99  (100%)
  • Tier C Masked (Repeat-masked US4, US7, US8, UL49 × 3 methods):       36 / 36  (100%)
  • Tier C Unmasked (Baseline comparison sets × 3 methods):              36 / 36  (100%)
  • Methods: FUBAR (Bayesian), BUSTED (Gene-wide LRT), aBSREL (Branch-site)
  • Biological Validation: Group 1 UL30 (PASS), Group 2 gB (PASS), 
                          Group 3 Conserved (PASS), Group 4 Masking (PASS)

========================================================================================
2. STAGE 2 — DEEP DIVERSIFYING SELECTION (62 / 62 Runs - 100% COMPLETE)
========================================================================================
  • Prioritized Targets: 31 high-priority vaccine antigens & regulatory genes
  • Methods: FEL (Fixed Effects Likelihood) & MEME (Mixed Effects Model of Evolution)
  • Codon-Level Mapping: Pervasive vs episodic selection at p < 0.10 and p < 0.05

========================================================================================
3. STAGE 3 — RECOMBINATION DETECTION (10 / 10 Runs - 100% COMPLETE)
========================================================================================
  • Method: HyPhy GARD (Genetic Algorithm for Recombination Detection)
  • Targets: Top-10 combined vaccine & regulatory genes (UL27, US6, UL22, UL1, 
            UL44, UL10, UL53, UL30, UL47, UL46)
  • Inferred: Breakpoints, ΔAICc improvements, partition trees, and clonal verdicts
========================================================================================
```

---

## 📁 Repository Directory Structure

```
HSV-Vaccine-Epitope/
├── README.md                                 # Master repository documentation
├── .gitignore                                # Excludes binary executables & cache
├── hsv_curated_dataset/
│   ├── 00_raw_downloads/                    # Raw NCBI GenBank genome dumps
│   ├── 01_reference_genomes/                # HSV-1 (strain 17) & HSV-2 (HG52) refs
│   ├── 02_filtered_genomes/                 # Quality-filtered non-chimeric genomes
│   ├── 03_metadata/                         # Isolate provenance & clinical metadata
│   ├── 04_cds_unaligned/                    # Extracted raw coding sequences
│   ├── 05_cds_qc_passed/                    # Stop/frameshift validated CDS
│   ├── 06_alignments/                       # In-frame codon alignments & IQ-TREEs
│   │   ├── HSV1/                            # HSV-1 cohort alignments & trees
│   │   ├── HSV2/                            # HSV-2 cohort alignments & trees
│   │   ├── combined/                        # Combined inter-clade alignments & trees
│   │   └── masked/                          # Repeat-masked alignments (US4, US7, etc.)
│   ├── 07_qc_reports/                       # Comprehensive HyPhy selection analyses
│   │   └── hyphy/
│   │       ├── stage1_summary.tsv           # Master Stage 1 screening summary
│   │       ├── stage1_biological_validation.tsv # 4-group biological validation audit
│   │       ├── stage1_outlier_investigation.tsv # Outlier resolution report
│   │       ├── stage2_trigger_genes_final.tsv   # 31 prioritized Stage 2 targets
│   │       ├── stage2_summary.tsv           # Master Stage 2 FEL/MEME results
│   │       ├── gard_summary.tsv             # Master 10-gene GARD recombination summary
│   │       ├── stage1/                      # 270 Stage 1 JSON output files
│   │       ├── stage2/                      # 62 Stage 2 JSON output files
│   │       └── gard/                        # 10 GARD recombination JSON output files
│   ├── 08_documentation/                    # Formal methods draft & pipeline logs
│   └── 09_scripts/                          # Reproducible Python execution pipelines
│       ├── step12_stage2_pipeline.py        # Automated Stage 2 FEL/MEME orchestrator
│       └── step12_gard_pipeline.py          # Automated GARD recombination orchestrator
└── _archived_discarded/                     # Pre-screening prototypes & archives
```

---

## 📈 GARD Recombination Screening Summary Table

| Gene | Target Antigen / Protein | Taxa | Length (nt) | Breakpoints | Location | $\Delta\text{AICc}$ | Recombination Status |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **`UL27`** | Glycoprotein B (gB) | 302 | 2,736 | 1 | nt 735 | 310.84 | **RECOMBINANT** |
| **`US6`** | Glycoprotein D (gD) | 145 | 1,182 | 0 | None | 0.00 | **STRICTLY CLONAL** |
| **`UL22`** | Glycoprotein H (gH) | 276 | 2,514 | 1 | nt 604 | 468.23 | **RECOMBINANT** |
| **`UL1`** | Glycoprotein L (gL) | 110 | 672 | 0 | None | 0.00 | **STRICTLY CLONAL** |
| **`UL44`** | Glycoprotein C (gC) | 213 | 1,533 | 0 | None | 0.00 | **STRICTLY CLONAL** |
| **`UL10`** | Glycoprotein M (gM) | 193 | 1,422 | 0 | None | 0.00 | **STRICTLY CLONAL** |
| **`UL53`** | Glycoprotein K (gK) | 142 | 1,014 | 0 | None | 0.00 | **STRICTLY CLONAL** |
| **`UL30`** | DNA Polymerase Catalytic Subunit | 312 | 3,723 | 1 | nt 2018 | 6,619.11 | **RECOMBINANT** |
| **`UL47`** | VP13/14 Major Tegument Protein | 227 | 2,106 | 0 | None | 0.00 | **STRICTLY CLONAL** |
| **`UL46`** | VP11/12 Tegument Phosphoprotein | 252 | 2,217 | 1 | nt 1713 | 2,385.71 | **RECOMBINANT** |

---

## 🛠️ Reproduction & Requirements

* **Python 3.10+**: `biopython`, `pandas`, `numpy`
* **HyPhy**: v2.5.103 (MP multi-threaded) via Docker: `docker pull hyphy/hyphy`
* **IQ-TREE**: v2.2+ for maximum likelihood phylogenies

To inspect or re-execute any stage with checkpoint caching:
```bash
# Run Stage 2 Selection Analysis (FEL & MEME)
python hsv_curated_dataset/09_scripts/step12_stage2_pipeline.py

# Run Top-10 GARD Recombination Screening
python hsv_curated_dataset/09_scripts/step12_gard_pipeline.py
```

---

## 📜 Authors & Citation
* Dataset curated and analyzed for rational HSV multi-antigen vaccine design.
* HyPhy models implemented according to Pond et al. (MBE 2006, 2012, 2015, 2020).
