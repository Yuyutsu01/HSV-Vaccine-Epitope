# Multi-method concordance is not a substitute for variant support: a five-method selection screen across 649 herpes simplex virus genomes

**Manuscript draft, 2026-11-21.**
Every quantity in this manuscript was recomputed from the named result table. Figure
callouts are marked `[FIGURE n NEAR HERE]` at the point in the text where each figure
should be placed. Section 11 restates the placement as a table for the production editor.

> **Revision note, 2026-11-21.** Assembling §3.12 exposed three defects in the
> conservation layer: reference residues were numbered on the masked alignment (8,072 of
> 22,955 numbers wrong), a peptide's between-type flags were set from however many of its
> residues had a conservation row, and between-type identity was judged on cohort
> consensus residues while candidate peptides are cut from reference proteins. All three
> are fixed, every affected stage was re-run, and the candidate numbers below are the
> corrected ones — most consequentially the cross-protective set, 2,120 → **1,729**
> peptides, and the minimal panel, four → **five**. The selection screen, the predictor
> benchmark and the conservation percentages were unaffected; §3.6 and Figures 2–5 are
> unchanged. `16_shortlist/excluded_not_in_other_type.tsv` records what the new check
> removed.

**Authors.** *(to be completed)*
**Affiliations.** *(to be completed)*
**Corresponding author.** *(to be completed)*

---

## Abstract

**Background.** Site-level tests for diversifying selection are routinely combined by
requiring two or more methods to agree, on the reasoning that concordance is conservative.
We tested that reasoning on herpes simplex virus, a system where within-type coding
sequences are 99.4–99.7% identical and where the assumption is therefore under maximum
strain.

**Methods.** We curated 1,305 NCBI Virus records into 649 analysis-ready genomes, extracted
and quality-controlled 12,054 coding sequences for 20 genes, and collapsed them to 4,163
unique sequences in 40 gene × type cohorts. HSV-1 and HSV-2 were analysed separately
throughout. Codon-aware alignments, occupancy masking, per-type maximum-likelihood trees
and GARD recombination screening produced 61 analysis units, on which five HyPhy methods
(FEL, MEME, FUBAR, BUSTED, aBSREL) were run in a single pass: 305 completed runs. Every
site called under selection was then annotated with the number of isolates actually
carrying the variant codon. Separately, class I, class II and B-cell epitope predictions
were generated, benchmarked where a curated truth set existed, and integrated into a
ranked candidate list.

**Results.** Gene-wide and branch-level selection is readily detected: BUSTED is
significant in 32 of 61 analysis units (24 of 40 cohorts, 18 of 20 genes) and aBSREL in 75
of 8,789 branches. At site level the picture inverts. FEL calls 154 sites at *q* < 0.05 and
every one has β < α — pervasive selection in these cohorts is uniformly purifying. MEME and
FUBAR together call 92 sites, of which exactly 2 satisfy a pre-registered two-method
concordance rule. Both fail on variant support: they rest on 4 of 293 and 2 of 276
isolates, their FUBAR posteriors sit a thousandth above the 0.95 cut-off, their MEME β̂
values are the unbounded estimates a near-invariant column produces, and their supporting
accessions include one sequence that fails CDS quality control and one pair submitted
consecutively in a single batch. Of the 92 candidate sites, 76 are robustly supported
(median 54 variant-carrying isolates) but almost all rest on FUBAR alone. **No site in this
study is both multi-method concordant and robustly supported.**

**Conclusions.** On low-divergence alignments, requiring methods to agree selects *for*
artifacts rather than against them, because the failure modes of MEME and FUBAR coincide at
near-invariant columns. Concordance and variant support are independent checks and both are
required. We report the negative site-level result, the positive gene- and branch-level
result that does not depend on single alignment columns, and a fully documented curated
dataset of 4,163 quality-controlled coding sequences with per-type alignments, coordinate
maps and a reproducible pipeline.

**Keywords.** herpes simplex virus; selection inference; HyPhy; multiple-testing
correction; negative results; comparative genomics.

---

## Impact statement

Combining site-level selection methods by requiring agreement is standard practice and is
widely described as conservative. We show that on highly conserved coding sequences it is
not conservative but actively misleading: the two methods most often combined fail in the
same place, at columns carrying a single non-synonymous change, so concordance concentrates
artifacts instead of filtering them. The correction is cheap — count how many sequences
carry the variant — and the consequence of omitting it is a published positive result that
rests on one or two sequences. The finding generalises beyond herpesviruses to any
low-divergence dataset, which increasingly means any within-species pathogen genomics
study.

---

## 1. Introduction

Detecting positive selection in viral coding sequences is a routine first step in antigen
prioritisation, and the standard toolkit — FEL, MEME, FUBAR, BUSTED and aBSREL as
implemented in HyPhy — is mature and well characterised. What is less well characterised is
how these methods behave together at the low end of the divergence range. Most published
benchmarks use datasets with substantially more variation than a within-species cohort of a
large double-stranded DNA virus provides.

Herpes simplex virus is an unusually demanding test case. HSV-1 and HSV-2 are distinct
viruses whose shared genes differ by 9–62% at the amino acid level, while within either
type the same genes are 99.4–99.7% identical. The two facts pull in opposite directions: a
pooled analysis has ample apparent variation but is a between-species comparison, and a
per-type analysis is biologically correct but leaves individual alignment columns almost
invariant. An earlier version of the present analysis pooled the two types and recovered a
recombination verdict for *US6*/gD that disappears entirely when each type is screened on
its own; that verdict was the species split being read as recombination, and the
pooled-cohort results were withdrawn rather than updated. This manuscript reports the
per-type analysis only.

Against that background we set out to do three things. First, to build a curated,
quality-controlled, per-type dataset for 20 HSV genes that others can reuse. Second, to run
a complete five-method selection screen on it with multiple-testing correction fixed in
advance. Third — and this became the principal result — to ask whether the conventional
practice of requiring two or more site-level methods to agree does what it is assumed to
do when the underlying alignments are this conserved.

It does not. The screen's two concordant sites are the two least credible sites it
produced. We report that as the main finding, report the negative site-level result it
implies, and report separately the gene- and branch-level signal, which is substantial and
does not depend on any single alignment column.

---

## 2. Methods

A complete methods document, with per-stage file provenance, is provided as
Supplementary Methods (`methods_draft.md`). What follows is sufficient to evaluate the
results; section numbers in parentheses refer to the full document.

### 2.1 Sequence retrieval and genome curation (§1)

Complete genome records for HSV-1 (NCBI taxid 10298) and HSV-2 (taxid 10310) were
retrieved from NCBI Virus with the NCBI `datasets` CLI, returning 1,305 records. Records
were excluded at genome level for falling outside the type-specific length window, carrying
more than 1% ambiguous bases, being annotated as laboratory-adapted, engineered,
recombinant or synthetic, or duplicating another record's content: 417 passed, 29 were
flagged but retained, and 859 were excluded, each decision with its reason in Supplementary
Table S3.

Virus type was verified computationally rather than taken from the database label. Each
genome was compared against both type references (HSV-1 NC_001806.2, HSV-2 NC_001798.2) and
assigned by average nucleotide identity; 80 records disagreed with their database label and
the computed assignment was used throughout.

A tier-2 set of 229 HSV-2 genomes — excluded at tier 1, overwhelmingly for 1–2% ambiguous
bases — was readmitted, because applying the ambiguity rule strictly leaves HSV-2 cohorts
too small for selection inference. 649 genomes entered per-gene extraction: 303 HSV-1, 117
HSV-2 tier 1 and 229 HSV-2 tier 2. **178 of the 286 HSV-2 accessions in the final cohorts
(62%) are tier-2 genomes and the HSV-1 cohorts contain none**, so HSV-2 inference rests on
assemblies of lower average quality. This is stated as a limitation throughout and is not
corrected for.

### 2.2 CDS extraction, quality control and deduplication (§2–§3)

Twenty genes were chosen a priori: the entry and surface glycoproteins (*US6*/gD,
*UL27*/gB, *UL22*/gH, *UL1*/gL, *US8*/gE, *US7*/gI, *US4*/gG, *UL53*/gK, *UL10*/gM,
*UL49A*/gN), the immediate-early and regulatory proteins (*RL2*/ICP0, *RS1*/ICP4,
*UL54*/ICP27), the tegument proteins (*UL48*/VP16, *UL46*, *UL47*, *UL49*/VP22) and the
enzymes (*UL23*/TK, *UL44*/gC, *UL30*/DNA polymerase).

Coding sequences were extracted by reference-guided annotation transfer and assessed
independently of extraction success. A CDS passed only if it began with ATG, ended in a
stop codon, had a length divisible by three, contained no internal stop codon and carried
ambiguous bases below the per-gene threshold, with a length-ratio window of 0.95–1.05
against the reference. **12,054 QC-passed CDS were carried forward.** Exact duplicates
within a gene × type cohort were collapsed to 4,163 unique sequences, with the number of
input isolates each representative stands for retained as a weight, so that every
population-frequency quantity downstream is computed over the 12,054 isolates represented
rather than over the unique haplotypes.

Two QC defects are reported rather than silently repaired. A rebuild step for the
repeat-resident genes *RL2* and *RS1* omitted the ambiguity and length-ratio tests, so 81
sequences carrying 0.5–16.8% ambiguous bases reached the final cohorts; occupancy masking
removes the repeat regions where that ambiguity concentrates, leaving ambiguous bases at
0.007–0.042% of alignment cells. Separately, 8 HSV-1 *US7*/gI sequences were retained
despite failing the length-ratio window, because the excess is a clean in-frame +84 nt
insertion with no internal stop — a real gI length polymorphism rather than a
mis-annotation. Both defects and their measured impact are documented in the Supplementary
Methods and in the project's corrections register.

### 2.3 Alignment, masking and phylogenetics (§4–§6)

Protein sequences were aligned with FAMSA (via `pyfamsa` 0.7.0) and the nucleotide
alignment produced by back-projecting each codon onto its aligned amino acid, so that every
alignment column is a codon column in one consistent frame. Each alignment was validated
against hard invariants — length divisible by three, no internal stop in any row, every
row's ungapped sequence identical to its input — and all 40 per-type alignments pass.

Columns were masked where fewer than 50% of clinical sequences carried a residue; 4,037
codon columns were dropped across the 40 cohorts, the median cohort losing none and the
maximum loss being 63.4% (*RS1*/HSV-1). A per-cohort coordinate map records, for every
column, the unmasked column, the masked codon column and the reference residue, so that all
residue-level claims are reported in reference numbering.

Maximum-likelihood trees were inferred per gene × type with IQ-TREE 2.4.0, model chosen by
ModelFinder, 1,000 ultrafast bootstrap replicates and identical sequences retained.
Combined HSV-1 + HSV-2 alignments were built but are used **only** for between-type
conservation comparison; they are never used for phylogenetics, recombination detection or
selection inference, for the reason given in §1.

### 2.4 Recombination screening (§7)

GARD was run per gene × type, on the same alignments used for selection and before the
selection analysis, so that any detected breakpoint partitions the alignment selection is
then run on. HyPhy 2.5.103(MP) was executed in a pinned Docker image. Sixteen of 40 cohorts
carry at least one breakpoint, giving 37 partitions; with the 24 unpartitioned cohorts this
yields **61 analysis units**. Breakpoint ΔAICc ranges from 3.70 (*UL44*/HSV-1, retained as
the conservative choice) to 1,769.78 (*RL2*/HSV-1). The *RS1*/HSV-1 search is flagged as
truncated: it was stopped by decision at its confirmed two-breakpoint model after
approximately 16 h and 5,215 models, so the two breakpoints reported are supported but a
third is not excluded.

### 2.5 Selection inference and multiple-testing correction (§8)

For each analysis unit, five HyPhy methods were run in a single pass with no gating of one
method on another's result: BUSTED, FEL, MEME, FUBAR and aBSREL. The significance policy
was fixed in advance:

- FEL and MEME site *p*-values are corrected by Benjamini–Hochberg **within each
  (gene, cohort, partition, method) group**; sites are reported at *q* < 0.05.
- BUSTED gene-level *p*-values are corrected by Benjamini–Hochberg across the family.
- aBSREL uses the method's own Holm–Bonferroni corrected *p*-values.
- FUBAR sites are reported at posterior probability ≥ 0.95, with ≥ 0.90 as a secondary,
  explicitly labelled tier.
- A site is called under selection in the primary analysis only if at least two methods
  agree, with the direction taken from β versus α.

The FEL and MEME grouping is consequential and is stated explicitly because a reader who
recomputes naively will get a different number: the 61 analysis units are separate
experiments on separate alignments, and pooling them would let a gene with many codons set
the threshold for a gene with few. Pooled Benjamini–Hochberg over all FEL rows at once
returns 95 sites at *q* < 0.05 rather than the 154 reported here.

### 2.6 Variant support (§8.1)

A significance threshold says nothing about how many sequences carry the variant that
produced the signal, and on alignments this conserved that question is decisive: a single
non-synonymous change at an otherwise invariant codon drives β/α towards infinity and MEME
will fire with an essentially unbounded β. Every reported site was therefore annotated with
its variant support, under a rule fixed before the outcome was examined: a site is
`NOT_ROBUST` if fewer than 5 of the isolates represented carry the variant, or if every
supporting sequence fails CDS quality control.

### 2.7 Epitope prediction, benchmarking and integration (§11–§17)

NetMHCpan-4.1 requires an academic licence that could not be obtained non-interactively, so
**MHCflurry 2.3.30** was used as the free, already-trained equivalent. Nothing was trained
or tuned in this work; every predictor is used as an instrument with its published models
and thresholds. All 9- and 10-mers were enumerated from the type reference proteins and
scored against each of the 27 alleles of the class I reference panel separately, giving
45,688 peptide predictions.

Validation used 688 IEDB epitopes curated to human-host records with a positive assay
outcome whose source antigen matched one of the 20 reference proteins exactly. Absence of
an IEDB record is never treated as a negative label and no negative epitope set is
constructed anywhere in this work; the benchmark therefore measures recovery of known
epitopes, not discrimination. Benchmarking is gene-grouped so that performance is not
inflated by one well-studied antigen.

**NetMHCIIpan-4.3** (licence-gated; run on the DTU web server and archived),
**BepiPred-3.0** and **DiscoTope-3.0** were run over the same reference proteins and
structures. None of the three is benchmarked here against an HSV-specific truth set, and
none of their published thresholds behaves on this data as documented, so **all three are
carried as continuous annotations and none is used as a filter anywhere in the pipeline**.

Candidates were assembled by joining predicted binding, per-residue conservation,
known-epitope flags, N-linked glycosylation sequons and the positive-selection filter. The
priority score is a fixed, pre-declared linear combination of percentile ranks
(0.45 × binding + 0.35 × conservation + 0.15 × conserved-in-both-types
− 0.05 × glycosylation-sequon overlap) and is an ordering device, not a trained model or a
predicted probability of immunogenicity.

Population coverage was computed from the IEDB allele-frequency data using the published
method of Bui *et al.* (2006): per locus, coverage = 1 − (1 − Σ*f*)², and combined
coverage = 1 − Π(1 − coverage_locus). The formula was reimplemented rather than the
downloaded tool executed, and the downloaded data file was read with a restricted
unpickler admitting only `collections.OrderedDict`, so no code from the download can
execute.

### 2.8 Reproducibility

Package versions are pinned (Python 3.11); IQ-TREE 2.4.0 and the NCBI `datasets` CLI are
bundled; HyPhy runs in a Docker image whose digest is recorded in the provenance JSON of
each stage. A Snakemake workflow declares every stage and its dependencies, and each stage
is also runnable as a standalone script. The analysed cohort is pinned by accession in
Supplementary Table S2. A provenance audit (`step30_provenance_audit.py`) checks that every
one of the 2,514 files in the repository has an entry in the data dictionary and fails if
one does not. All figures are regenerated from the result tables by a single script; no
value on any figure is typed in by hand.

---

## 3. Results

### 3.1 A per-type dataset of 4,163 coding sequences

After genome filtering, computational type verification, per-gene CDS extraction and QC,
and collapse of exact duplicates, **4,163 unique coding sequences** remain across **40 gene
× type cohorts** for 20 genes, representing **12,054 isolates**. The full pipeline, with
the count surviving each stage, is shown in Figure 1.

> **[FIGURE 1 NEAR HERE]** — the study design and data flow. This is the figure a reader
> needs before any number in this section, and it should sit immediately after this
> paragraph, before the discussion of cohort size. It also marks the withdrawn v1 branch
> and the two-track split at the selection stage, both of which the text returns to.

Cohort sizes are uneven, and the unevenness constrains what the HSV-2 arm can support.
HSV-1 cohorts run to 215 unique sequences (*RS1*); the smallest cohorts anywhere are
*UL49A* at 28 (HSV-2) and 32 (HSV-1) and *UL1* at 32 (HSV-2), and the two repeat-resident
genes are thin in HSV-2 as well (*RS1* 57, *RL2* 75). Per-gene cohort sizes are in
Figure S1.

> **[FIGURE S1 NEAR HERE]** — supplementary; cited here and not reproduced in the main
> text, because cohort size per gene is reference material rather than argument.

Two properties of the cohort must travel with every HSV-2 result. First, 62% of HSV-2
accessions come from the relaxed tier-2 set (§2.1). Second, patient identity is not recorded
in these records, so the cohorts are not samples of independent infections; serial isolates
from one patient are present and are weighted, not removed.

### 3.2 Recombination is type-specific, and the pooled analysis was wrong

GARD found at least one breakpoint in 16 of 40 cohorts, yielding 37 partitions and, with
the 24 unpartitioned cohorts, 61 analysis units.

The negative result here is as informative as the positive ones. **US6/gD shows no
breakpoint in either type.** The earlier pooled-cohort analysis had returned a recombination
verdict for *US6*; running GARD within each type removes it. That verdict was the
HSV-1/HSV-2 species split being read as recombination, and it is the direct evidence that
the pooled-cohort results were invalid — which is why they are withdrawn rather than
updated.

### 3.3 These are strongly conserved sequences, and that is the central constraint

Mean within-type identity is **99.38% (HSV-1)** and **99.73% (HSV-2)**; between-type
identity ranges from **37.5% (*US4*/gG)** to **90.8% (*UL30*/pol)**. 281 windows of ≥ 9
amino acids are conserved within both types and identical between them.

> **[FIGURE 2 NEAR HERE]** — within-type and between-type identity across the 20 genes.
> Place it with this paragraph. Panel (a) carries the within-type figures that frame the
> power limitation discussed below; panel (b) carries the between-type range and the
> *US4*/gG outlier the next paragraph explains.

*US4*/gG is confirmed as a poor cross-protective target, but for a reason that is
arithmetic rather than biological: HSV-1 gG is 238 aa against 699 aa in HSV-2 and occupies
only 238 of 698 alignment columns, so its apparent identity is dominated by non-overlap.

The high within-type identity is the central constraint on this study's power. A selection
analysis on sequences this conserved has limited ability to detect site-level diversifying
selection. That is a property of the virus, and it is stated as such rather than
compensated for by relaxing thresholds.

### 3.4 Gene-wide and branch-level selection is readily detected

All **305 HyPhy runs** (61 units × 5 methods) completed, each behind a completion marker.

BUSTED detects gene-wide episodic diversifying selection in **32 of 61 analysis units** at
*q* < 0.05 — equivalently 24 of the 40 gene × type cohorts and 18 of the 20 genes. Six genes
(*RS1*, *UL22*, *UL30*, *US4*, *US7*, *US8*) are significant in both types. aBSREL
identifies **75 of 8,789 branches** under episodic selection at Holm-corrected *p* < 0.05,
concentrated in *RS1*/HSV-1 (13 branches), *RL2*/HSV-1 (7), *UL44*/HSV-2 (6) and
*US7*/HSV-1 (6).

> **[FIGURE 3 NEAR HERE]** — gene-wide *q*-values for the significant units and
> branch-level counts. Place it directly after this paragraph. Panel (a) shows the 32
> significant units only; the 29 sub-threshold units are omitted for legibility and are
> listed in the result table.

These are lineage- and gene-level signals. They do not depend on any single alignment column
and they are unaffected by the site-level problems described next.

### 3.5 Pervasive site-level selection is uniformly purifying

FEL calls 154 sites at *q* < 0.05, and **every one has β < α**. No site-level positive
selection is contributed by FEL at all. For coding sequences at 99.4–99.7% identity this is
the expected result, and it frames everything that follows.

### 3.6 The two concordant sites are the two least credible sites in the screen

MEME and FUBAR together call **92 sites** as positively selected by at least one method (87
FUBAR, 7 MEME), and exactly **2 sites** satisfy the pre-registered ≥ 2-method concordance
rule: HSV-1 *RL2*/ICP0 residue 308 and HSV-1 *US7*/gI residue 155, both FUBAR + MEME.

Neither survives a check on how many sequences carry the variant (Table 2).

**Table 2. The two method-concordant sites, and why neither is reported as a finding.**

| | *RL2*/ICP0, HSV-1, residue 308 | *US7*/gI, HSV-1, residue 155 |
|---|---|---|
| Majority codon | GGC (Gly) × 183 | CCA (Pro) × 125 |
| Variant codons | CGC × 1, TCC × 2 | GAC × 1, TCA × 1 |
| Isolates carrying a variant | **4 of 293 (1.4%)** | **2 of 276 (0.7%)** |
| Largest codon change | 2 nt | **3 nt** |
| FUBAR posterior | 0.9560 | 0.9523 |
| MEME *q* (β̂) | 0.0246 (1,024) | 0.0051 (1,842) |
| FEL | *p* = 0.169, n.s. | *q* = 0.637, n.s. |
| Supporting accessions | HM585511, OP431901, **OR771687 — fails CDS QC** | **PP379031, PP379033 — consecutive accessions, one submission batch** |
| Verdict | `NOT_ROBUST` | `NOT_ROBUST` |

Every indicator points the same way. Support is two to four isolates out of roughly 290.
The FUBAR posteriors sit a thousandth above the 0.95 cut-off. The MEME β̂ values of 1,024
and 1,842 are the unbounded estimates a near-invariant column produces, not biological rate
estimates. FEL supports neither. A three-nucleotide codon change observed in a single
isolate is far more consistent with a sequencing or assembly artifact than with a point
substitution. One site's only two variants come from consecutive accessions in a single
submission; the other's include a sequence that fails the CDS ambiguity criterion — the
same QC defect described in §2.2, surfacing at exactly the level where it matters.

**Both sites are therefore reported as candidate artifacts and not as evidence of positive
selection.** The pre-registered rule's output is left unchanged in
`selected_sites_final.tsv` so that it stands on the record, and the robustness table
qualifies it.

The wider pattern inverts the concordance rule's suggestion. Of the 92 sites called by at
least one method, **76 are robust** — a median of 54 variant-carrying isolates against a
median cohort of 297 — while 7 are weak and 9 are not robust. The robust ones are almost
entirely FUBAR-only: well supported in the data, but resting on a single method.

> **No site in this study is both multi-method concordant and robustly supported.**

> **[FIGURE 4 NEAR HERE]** — the headline figure, and it must appear on the same spread as
> this statement. Each of the 92 candidate sites is plotted against the number of isolates
> carrying its variant; shape encodes the variant-support verdict and colour the virus
> type. Both concordant sites fall to the left of the five-isolate support threshold, and
> the upper panel shows that significance is unrelated to support across the whole screen.

This is the honest summary of the site-level screen, and it is also a methodological result
in its own right. Requiring two methods to agree — intended as the conservative choice —
*selected for* artifacts here, because MEME fires on singleton non-synonymous changes at
near-invariant columns while FUBAR's posterior grazes threshold at the same columns. The
two methods fail in the same place, so their agreement carries far less information than
agreement between independent checks would. Concordance and variant support are independent
checks, and both are needed.

### 3.7 Predictor validation is strong in HSV-1 and weak in HSV-2

Over 27 scorable gene × type groups spanning 17 genes, predicted affinity gives **mean AUC
0.809** (median 0.892, range 0.425–0.997), mean average precision 0.326, and **6.0-fold
enrichment** of known epitopes in the top decile over 428 positives.

The split by type is stark and is a limitation rather than a nuisance: **HSV-1 groups give
mean AUC 0.931 and HSV-2 groups 0.600**. HSV-2 class I records in IEDB are fewer and are
largely inherited from HSV-1 cross-reactivity studies, so predictions for HSV-2 are
correspondingly less well validated.

> **[FIGURE 5 NEAR HERE]** — AUC per gene × type, with marker area encoding the number of
> curated IEDB positives behind each estimate. Place it here. The size encoding is load
> bearing: the two lowest HSV-2 AUCs are computed over some of the largest positive sets in
> the panel, so the weak HSV-2 performance is not an artifact of thin truth sets.

The control that matters most: **conservation alone performs near chance** (mean AUC 0.566,
1.38-fold top-decile enrichment). Conservation is therefore not a proxy for immunogenicity,
which is the direct justification for using predicted binding and conservation as two
independent filters rather than either alone.

### 3.8 Candidate shortlist

Joining predicted binding, per-residue conservation, IEDB flags, glycosylation sequons and
the positive-selection filter gives 12,771 ranked candidates. The primary shortlist —
predicted binder, every residue ≥ 99.0% conserved within type, not positively selected — is
**9,488 rows / 8,010 distinct peptides**, of which **3,207 rows / 1,729 distinct peptides**
are identical across both types and form the cross-protective subset. A row is one peptide
in one virus type, so a peptide may appear twice; 1,478 of the 1,729 cross-protective
peptides have a row in each type and 251 in only one, because a peptide present and
identical in both can still clear the within-type conservation threshold in only one of
them. The distinct-peptide counts are the ones to quote as candidate epitopes.

Cross-protective status is verified directly against the two reference proteins rather
than inferred from the alignment: a peptide qualifies only if it occurs, contiguously, in
the other type's protein. That check is not redundant with per-residue identity, because
an indel inside the window makes every residue correspond while the contiguous peptide
exists in one type only — it removed 2 peptides here
(`16_shortlist/excluded_not_in_other_type.tsv`), and before the conservation layer was
corrected it would have removed 444.

Two quality flags travel with any use of this list. **1,460 of the 9,488 primary rows
(15.4%) are low-complexity**, drawn from the serine-, glycine-, alanine- and proline-rich
tandem repeats of ICP4 and ICP0 and the repeat region of gG — 1,307 distinct peptides,
concentrated in *RS1* (267 rows), *RL2* (137) and *US4* (113), such as `AAAAAADLLF` and
`AAAAPVVSA`. These are perfectly conserved and score as binders but are not credible
epitopes. Six appear in the top 50 by priority score, so the ranking suppresses them only
partly; they are flagged rather than dropped so that the shortlist size is not changed
silently by a filter added after the fact, and they are barred from the minimal panel in
§3.9, where it matters. Separately, **the positive-selection exclusion removes only 5 rows**,
and §3.6 shows the two sites driving it are probably artifacts: this filter is close to a
no-op on these data and should not be presented as having shaped the shortlist.

As a positive control, the cross-protective shortlist recovers **`RMLGDVMAV`** — a published
HLA-A\*02 gB epitope — at **rank 13 of 3,207** by priority score, and rank 7 once
low-complexity repeats and duplicate rows are set aside, correctly flagged as IEDB-known
(gB 559–567 in HSV-2 numbering, 562–570 in HSV-1).

The top of the ranking is dominated by HLA-A\*02:03 and the A2 supertype. That is a known
property of class I predictors on hydrophobic 9-mers, not a biological finding, which is
why coverage is assessed per allele rather than from the top of the list.

### 3.9 Population coverage, and what it does and does not measure

Across **239 class I populations**, the cross-protective shortlist reaches **99.71% World
coverage**, with a median of 97.5% across populations and **73.6% of populations at or
above 90%**. A greedy minimal panel reaches that same World coverage with five peptides
(Table 3).

**Table 3. Greedy minimal class I panel, low-complexity peptides excluded.**

| Step | Peptide | Alleles bound | New alleles | Cumulative alleles | World coverage (%) |
|---|---|---|---|---|---|
| 1 | `ILDNPGFVTF` | 13 | 13 | 13 | 90.58 |
| 2 | `EVIDKINAK` | 8 | 8 | 21 | 99.22 |
| 3 | `AMAAGRFGW` | 9 | 3 | 24 | 99.54 |
| 4 | `AQTGVRVAL` | 7 | 2 | 26 | 99.65 |
| 5 | `AAQYRNAVV` | 1 | 1 | 27 | 99.71 |

Low-complexity peptides are barred from the panel search: 1,307 are barred across the
candidate pool, leaving 1,559 of the 1,729 cross-protective candidates eligible. The
exclusion costs nothing measurable — the unconstrained search also needs five peptides and
also reaches 27 alleles and 99.71%, differing only in using the low-complexity `AAAVPLLSA`
— and the unconstrained panel is retained for comparison.

The panel is a coverage construct, not a recommendation. Two peptides from its earlier
four-peptide version, `AEYGGSFRF` and `AAFDGGSTA`, were removed by the corrections in the
revision note: both were flagged cross-protective but occur in HSV-2 only. That is the
clearest illustration of why the direct presence check in §3.8 is necessary — it was the
headline deliverable that the defect had reached.

Class II coverage over the 64 populations typed at all five presenting loci is **98.87%
conservative** (median 99.63%, minimum 90.75%, in the population the IEDB file records as
`Korea; South`) to **99.98% upper bound**, with World coverage of 99.64–100%. Combined with
class I across the 51 populations present in both frequency sets, the mean is 99.98%.

> **[FIGURE 6 NEAR HERE]** — three panels: class I coverage across the 239 populations,
> the class II conservative-to-upper-bound bracket across the 64 fully typed populations,
> and the greedy minimal panel. Place it with this subsection. Panel (b) is the one to
> point readers at, because the bracket and the omitted-population count are the honest
> content of the class II analysis.

Four caveats are reported with these figures, and the third is the largest single
correction made in this analysis.

1. **HLA-C is absent from the 27-allele panel**, so combined class I coverage is
   conservative.
2. **DP and DQ coverage is bracketed, not computed.** The IEDB file stores class II
   frequencies per chain, while a DP or DQ molecule is a heterodimer and NetMHCIIpan
   predicts for the αβ pair. The upper bound treats each chain locus as contributing
   independently and overestimates; the conservative figure requires both chains and
   underestimates, because DQA1 and DQB1 are in strong linkage disequilibrium. Both are
   reported and neither is presented as *the* class II coverage.
3. **Only populations typed at all five presenting loci are summarised.** Of 274
   populations, 64 have all five. A population with no HLA-DRB1 record is not a population
   with zero DR coverage; it is one that cannot be scored at DR. Treating an untyped locus
   as zero puts the mean conservative coverage at **82.3% instead of 98.9%** — a 16-point
   artifact that would have read as a geographic finding. HLA-DRB3, -DRB4 and -DRB5 have no
   frequency data in the file at all, so 4 of the 27 panel alleles cannot contribute and
   the DR term is DRB1-only.
4. **Coverage saturates, and it describes the allele panel rather than the shortlist.**
   Because the metric counts alleles with at least one binder, all shortlist tiers give
   identical coverage once every panel allele is hit. Coverage discriminates *panels*, not
   shortlist depth. Low-coverage populations are an artifact of sparse typing, not biology:
   Central America returns 7.76% because very few alleles are typed in that dataset.

### 3.10 Structural mapping

Candidates in the structurally characterised entry glycoproteins were mapped onto 2C36
(gD), 2GUM (gB) and 3M1C (gH–gL): **3,145 residues across 7 chains** assigned at
0.995–1.000 identity, of which 1,099 are solvent-exposed (RSA ≥ 0.25), and **32 conserved,
exposed surface patches** of ≥ 5 contiguous residues were identified.

The most promising is **gB 257–264 (`SRVEAFHR`)** in HSV-1 numbering, 252–259 in HSV-2,
where it occurs contiguously and identically: fully conserved within both types,
solvent-exposed, and annotated as inside a known IEDB epitope at all eight residues.

> **[FIGURE 7 NEAR HERE]** — gB with the `SRVEAFHR` patch highlighted, and gD coloured by
> DiscoTope score. Place it here. Panel (a) is the paper's strongest structural claim and
> panel (b) is the only structural view of the B-cell arm.

No positively selected site falls within any of the four structurally covered proteins, so
the structural annotation reports `positively_selected = False` for all 3,145 residues.
This is a limit of structural coverage, not evidence about those residues.

### 3.11 B-cell and class II predictions, and the behaviour of their thresholds

**BepiPred-3.0** was run on all 40 reference proteins (23,184 residues), with a hard
register check — every residue letter in the predictor output must equal the reference
residue at the same position — that all 23,184 residues passed. 8,647 residues (37.3%)
clear the published default threshold of 0.1512, a proportion high enough that the
threshold is reported but not used.

**DiscoTope-3.0** was run on the same three structures, mapping 3,145 of 3,145 residues onto
reference numbering with no unmapped rows and no residue-letter mismatches. Its published
high-confidence threshold of 0.90 is **never reached**: scores span 0.0005–0.659 (median
0.071), and only 9 residues reach 0.50.

> **[FIGURE S2 NEAR HERE]** — supplementary; cited here. It shows both score distributions
> and what their published thresholds do on this data, which is the evidence for carrying
> both predictors as annotations rather than filters.

Two internal checks are worth reporting. Crystallographic replicate chains of the same
protein, scored independently, agree at *r* = 0.956 (gD, *n* = 274) and *r* = 0.886–0.904
(gB, *n* = 570), which bounds how precise any single-structure score can be. Residues inside
a curated IEDB epitope score higher than the rest (median 0.089 vs 0.061, Mann–Whitney
one-sided *p* = 5.3 × 10⁻¹⁵) — genuinely out of sample, since DiscoTope has no access to
IEDB, but a small effect that is not a performance estimate. The conformational layer covers
4 of 20 genes, because only those have a solved structure; for the other 16 the result is
missing, not negative.

**NetMHCIIpan-4.3** produced 610,848 predictions (22,624 15-mers × 27 alleles), every one
verified to sit at its reported position in the reference protein. 4,861 peptide–allele
pairs (0.8%) are strong binders at the server's own 1% threshold and 25,855 (4.2%) are weak
or better; 0.8% against a 1% threshold is a sanity check on the rank column rather than a
finding. The summary statistic is not usable as a filter: the best percentile rank is a
minimum over 27 alleles × every overlapping 15-mer, so it falls below 1% for 76% of
candidates. Per-allele binding is the quantity to use.

---

### 3.12 The candidate target set

This section states what the resource contains, so that a reader deciding whether to
synthesise anything can see the whole basis for the selection in one place. Everything in
it is predicted. Nothing in it is measured, and §5 should be read alongside it.

**What survives every filter.** Of 12,771 ranked candidates, **1,729 distinct peptides** are
predicted binders, ≥ 99.0% conserved at every residue within their type, not positively
selected, and present contiguously and identically in both HSV-1 and HSV-2. Excluding the
low-complexity repeats leaves **1,559 usable cross-protective candidates**, and **480 of
those 1,559 (30.8%) lie inside a curated IEDB epitope**. That last figure is the arm's
strongest internal evidence, because IEDB played no part in the ranking: the priority score
is binding, conservation, cross-type identity and a glycosylation penalty, and the
known-epitope flag is carried as an annotation only.

**Table 4. Top of the cross-protective ranking, low-complexity repeats excluded.** One row
per distinct peptide; position is in the numbering of the type shown.

| # | Peptide | Gene | Position | Best allele | Inside an IEDB epitope |
|---|---|---|---|---|---|
| 1 | `FLAASALGV` | *UL22*/gH | HSV-2 805 | A\*02:03 | — |
| 2 | `MLFDQGTLV` | *RL2*/ICP0 | HSV-1 739 | A\*02:03 | yes |
| 3 | `FLSTLPSDV` | *UL48*/VP16 | HSV-2 102 | A\*02:03 | — |
| 4 | `FLGSCICFI` | *US7*/gI | HSV-2 274 | A\*02:03 | yes |
| 5 | `FLIAYQPLL` | *UL27*/gB | HSV-1 448 | A\*02:03 | yes |
| 6 | `MLIKGVDLV` | *UL30*/pol | HSV-1 950 | A\*02:03 | yes |
| **7** | **`RMLGDVMAV`** | *UL27*/gB | HSV-2 559 | A\*02:03 | **yes — positive control** |
| 8 | `ALGYQLAFV` | *UL22*/gH | HSV-2 471 | A\*02:03 | — |
| 9 | `FILDNPGFV` | *UL30*/pol | HSV-2 304 | A\*02:03 | — |
| 10 | `YLCDNFCPA` | *UL30*/pol | HSV-2 282 | A\*02:03 | — |

**Read the top of this table with care, and the whole list by allele.** Every one of the
ten binds HLA-A\*02:03, and that is a property of class I predictors on hydrophobic 9-mers
with aliphatic anchors, not a finding about HSV. The ranking is an ordering device over a
pre-declared formula, not a predicted probability of immunogenicity, so "rank 1" means
only "highest under that formula". The defensible use of the list is per allele and per
gene, which is why coverage (§3.9) is computed over the 27-allele reference panel rather
than from the head of the ranking.

**Where the candidates are.** The cross-protective set concentrates in *UL30*/pol (356
distinct peptides), *UL27*/gB (250) and *RS1*/ICP4 (204), then *UL47* (120), *UL48*/VP16
(92) and *UL46* (88). Part of that is real conservation and part is simply protein length —
these are among the longest products in the panel. At the other end, *US4*/gG yields 1 and
*UL1*/gL 6, consistent with gG being the most divergent product in the genome (§3.3).

**The three defensible claims.** Taken together the arm supports exactly three statements,
and it is worth separating them from the list as a whole:

1. **The ranking recovers known biology.** `RMLGDVMAV`, a published HLA-A\*02 gB epitope,
   appears at rank 7 of the de-duplicated cross-protective list without IEDB contributing
   to the score, and 30.8% of the usable set falls inside curated epitopes.
2. **A five-peptide panel covers the reference allele panel.** 99.71% World coverage, all
   27 alleles, with the caveats in §3.9 — HLA-C absent, coverage saturating, and the metric
   describing the panel rather than the shortlist.
3. **One structural target scores on every criterion.** gB `SRVEAFHR` (HSV-1 257–264,
   HSV-2 252–259): conserved within both types at every residue, contiguous and identical
   between them, solvent-exposed, and inside a curated IEDB epitope at all eight residues.
   It is the one candidate in this work that a structural immunologist could act on
   directly.

**What the set is for.** It is a conservation-filtered, cross-type, coverage-assessed,
IEDB-annotated starting pool from which 20–50 peptides could be chosen for synthesis and
assay. It is not a vaccine design, and the honest reading of §3.7 limits it further: the
affinity predictor ranks known epitopes well in HSV-1 (mean AUC 0.931) and barely better
than chance in HSV-2 (0.600), while conservation alone is near chance in both (0.566). Every
HSV-2 entry therefore rests on a predictor this study could not validate for HSV-2, and the
CD4 and B-cell layers (§3.11) are unvalidated throughout. The set is offered as a resource
with those limits attached, which is why this manuscript leads with §3.6 and not with this
table.

## 4. Discussion

### 4.1 Why concordance failed

The practice of requiring two or more site-level methods to agree rests on an implicit
assumption of independence: that the methods' false positives are uncorrelated, so that
agreement is evidence. On these data that assumption is false in a specific and predictable
way.

MEME tests for episodic selection and estimates β without an upper bound. At a column where
292 of 293 sequences carry one codon and one carries another non-synonymous codon, the
maximum-likelihood β̂ runs away — the values observed here are 1,024 and 1,842 — and the
likelihood-ratio statistic follows. FUBAR places a Bayesian posterior on the same column
and, with essentially one informative observation, returns a posterior that grazes any
threshold placed near the top of its range; both observed posteriors are within 0.006 of
the 0.95 cut-off. The two methods are not failing independently. They are failing on the
same feature of the data, from opposite directions, and concordance therefore concentrates
that failure rather than filtering it.

FEL, which tests for pervasive rather than episodic selection and is far less sensitive to
singletons, supports neither site. That is the useful signal in the concordance pattern, and
it is lost by a rule that counts methods rather than weighing them.

The remedy is not a higher threshold. Raising the FUBAR posterior cut-off to 0.96 would
remove both sites, but only by accident; the next dataset's singleton will sit at 0.97. The
remedy is to ask a question the significance test does not ask — how many sequences carry
the variant — and to treat that as an independent check rather than a post-hoc filter. In
this screen the check is decisive: it removes both concordant sites and retains 76 sites
with a median of 54 variant-carrying isolates, which are the sites a reader should look at
next even though each rests on one method.

### 4.2 Generality

Nothing in the mechanism is specific to herpesviruses. It requires only that alignment
columns be near-invariant, which is the normal condition for within-species cohorts of
slowly evolving organisms and increasingly the normal condition for pathogen genomics as
sequencing depth outruns diversity. We would expect the same behaviour in any dataset where
the typical column carries zero or one non-synonymous change, and the diagnostic is
cheap: plot significance against variant support, as in Figure 4, and look for points in
the top-left.

### 4.3 What the positive results do and do not show

The gene- and branch-level results are unaffected by this problem, because neither depends
on a single column. BUSTED's significance in 32 of 61 units and aBSREL's 75 branches are
real signals, concentrated in the immediate-early regulatory genes and the glycoproteins
with known antigenic variation. They are not, however, evidence about any particular
residue, and they should not be read as such; the mapping from a gene-level verdict to a
vaccine-relevant site is exactly the step this study shows cannot currently be made on
these data.

### 4.4 The resource

The curated dataset is offered as the third contribution: 4,163 quality-controlled coding
sequences in 40 per-type cohorts, codon-aware alignments with coordinate maps to reference
numbering, per-type trees, GARD partitions, the complete five-method HyPhy output, and a
documented, audited pipeline. Two QC defects are documented with their measured impact
rather than silently repaired, and the analysed cohort is pinned by accession so that the
dataset is the one the results were computed on.

---

## 5. Limitations

These are stated as constraints on interpretation, not as caveats to be discounted.

- **No site-level positive selection is established.** The screen's two concordant sites
  fail a minimum-variant-support check, and no site is both concordant and robust. This is
  the result most likely to be mis-cited.
- **The HSV-2 arm rests on lower-quality assemblies.** 62% of HSV-2 accessions come from the
  relaxed tier-2 set; the HSV-1 cohorts contain none.
- **The CD4 arm is prediction, not evidence.** NetMHCIIpan-4.3 was run over the 27-allele
  class II reference set, but no curated set of experimentally validated HSV class II
  epitopes was assembled to benchmark against. No peptide is claimed to be a CD4 epitope.
- **The B-cell arm is prediction, not evidence.** Neither predictor is benchmarked against
  experimentally mapped HSV B-cell epitopes here, and neither threshold behaves as
  published on this data — BepiPred's flags 37% of all residues, DiscoTope's flags none.
- **No immunogenicity measurement.** Every candidate is predicted. The IEDB benchmark shows
  the predictor ranks known epitopes well in HSV-1 and poorly in HSV-2; it does not show
  that any shortlisted peptide is immunogenic.
- **No epidemiological inference.** The cohorts are convenience samples with unknown
  geographic composition — country is unknown for 82.45% of records — and unknown patient
  structure.
- **The *RS1*/HSV-1 recombination search is truncated**, and a third breakpoint is not
  excluded.
- **Structural coverage is 4 of 20 genes.** For the other 16 the conformational result is
  missing, not negative.

---

## 6. Conclusions

A five-method selection screen over 61 analysis units and 305 HyPhy runs, on herpes simplex
virus coding sequences at 99.4–99.7% within-type identity, detects substantial gene-wide and
branch-level selection and no credible site-level positive selection. The two sites that
satisfy a conventional two-method concordance rule are the two least credible sites the
screen produced, because the two methods that agree fail on the same feature of the data.
Concordance is not a substitute for variant support, and on low-divergence alignments it is
not even conservative. We recommend that site-level selection scans on such data report
variant support alongside significance as a matter of course, and we provide the curated
dataset, the pipeline and the figures that make that check reproducible.

---

## 7. Data and code availability

All result tables, alignments, trees, coordinate maps and figures are in the project
repository. Supplementary tables S1–S5 give the sequence counts surviving each filter (S1),
the final accession list with isolate weights (S2), every removed sequence with its reason
(S3), per-gene cohort sizes (S4) and software versions (S5).

**Two items must be completed before submission.** First, the 587 MB of raw HyPhy JSON
output should be deposited in Zenodo and the versioned DOI cited here rather than a
repository URL, which journals increasingly require and which protects against the
repository changing after review. Second, the repository must be confirmed public and the
commit referenced here must be the one pushed.

---

## 8. Main figures

All figures are produced by `09_scripts/step27_figures.py` from the result tables; no value
on any figure is typed in by hand. Each is written as 300 dpi PNG and as vector PDF.

**Figure 1. Study design and data flow.**
`20_figures/fig1_study_design.*`
The pipeline from 1,305 NCBI Virus records to the final candidate shortlist, with the count
surviving each stage read from the result tables. The split below the HyPhy stage is the
paper's argument in one picture: the same 305 runs give a positive gene- and branch-level
result (left, blue) and a negative site-level result (right, orange). The withdrawn v1
analysis, which pooled HSV-1 and HSV-2 into single alignments, is marked at the foot rather
than omitted.

**Figure 2. Sequence conservation within and between types.**
`20_figures/fig2_conservation.*`
(a) Mean within-type nucleotide identity for each of the 20 genes, HSV-1 and HSV-2 plotted
on the same row and joined by a connector. (b) Mean between-type identity, directly
labelled. *US4*/gG is the lowest at 37.5%, for the length-mismatch reason given in §3.3.

**Figure 3. Selection is detected gene-wide and on branches.**
`20_figures/fig3_selection_gene_branch.*`
(a) BUSTED −log₁₀(*q*) for the 32 of 61 analysis units significant at *q* < 0.05; the 29
sub-threshold units are omitted for legibility and listed in `gene_level_results.tsv`.
(b) Number of branches under episodic selection per gene × type at Holm-corrected
*p* < 0.05 (aBSREL), 75 of 8,789 branches in total. Gene labels carry the virus type, since
eight genes appear once per type.

**Figure 4. Significance is unrelated to how many isolates carry the variant.**
`20_figures/fig4_site_variant_support.*`
Each of the 92 candidate sites is plotted against the number of isolates carrying its
variant codon (log scale). Upper panel: the sites with a MEME *q* value. Lower panel: the
FUBAR-only sites, which have no *q*; their vertical spread is jitter and encodes nothing.
Marker shape is the variant-support verdict and colour is virus type. The dashed line is the
five-isolate support threshold. Both method-concordant sites are labelled and both fall to
the left of it. **This is the paper's central figure.**

**Figure 5. Predictor validation is strong in HSV-1 and weak in HSV-2.**
`20_figures/fig5_predictor_benchmark.*`
AUC of predicted affinity against curated IEDB epitopes, per gene × type group, with marker
area proportional to the number of curated positives behind each estimate (5 to 75 per
group). Mean AUC is 0.931 for HSV-1 groups and 0.600 for HSV-2.

**Figure 6. Population coverage of the cross-protective shortlist.**
`20_figures/fig6_population_coverage.*`
(a) Class I coverage across 239 populations, sorted. (b) Class II coverage across the 64
populations typed at all five presenting loci, shown as the band between the conservative
and upper-bound estimates; the 210 populations with an untyped locus are omitted, because
scoring an untyped locus as zero produces a false zero. (c) The greedy minimal class I
panel, which reaches 99.71% World coverage and all 27 panel alleles with five peptides.

**Figure 7. Shortlisted surface on the deposited structures.**
`20_figures/fig7_structure.*`
(a) gB (PDB 2GUM, chain A) with the conserved exposed patch `SRVEAFHR` (gB 257–264)
highlighted; this patch is conserved in both types at every residue, solvent-exposed, and
inside a curated IEDB epitope at every residue. (b) gD (PDB 2C36, chain A) coloured by
DiscoTope-3.0 conformational B-cell score. Van der Waals spheres in a principal-axis
projection, painted far to near; the panels are at different magnifications and the scale
bars are what makes them comparable.

## 9. Supplementary figures

**Figure S1. Per-gene cohort size.**
`20_figures/figS1_cohort_size.*`
Unique coding sequences per gene × type after deduplication; 4,163 in 40 cohorts.

**Figure S2. B-cell predictor scores and the behaviour of their thresholds.**
`20_figures/figS2_bcell_scores.*`
(a) BepiPred-3.0 linear B-cell score across 23,184 residues, split by whether the residue
was called an epitope. The two distributions overlap across the published 0.1512 threshold:
residues called epitopes run down to 0.054 and residues not called run up to 0.257. (b)
DiscoTope-3.0 conformational score per chain, with the published high-confidence cut-off of
0.90 marked. No residue in the three structures reaches it; the highest score is 0.659.
Together these panels are the evidence for carrying both predictors as continuous
annotations and never as filters.

## 10. Tables

**Table 1. Selection screen: what was tested and what was called.**

| Level | Method | Tested | Called under selection |
|---|---|---|---|
| Gene-wide | BUSTED | 61 analysis units | **32** at *q* < 0.05 |
| Branch | aBSREL | 8,789 branches | **75** at Holm-corrected *p* < 0.05 |
| Site, pervasive | FEL | 23,155 codons | **154** at *q* < 0.05, **all with β < α** |
| Site, episodic / Bayesian | MEME, FUBAR | 23,155 codons | **92** by ≥ 1 method (87 FUBAR, 7 MEME) |
| Site, primary analysis | ≥ 2 methods concordant | 92 candidates | **2**, both `NOT_ROBUST` |
| Site, concordant **and** robust | — | 92 candidates | **0** |

**Table 2.** The two method-concordant sites — reproduced in §3.6.

**Table 3.** Greedy minimal class I panel — reproduced in §3.9.

**Table 4.** Top of the cross-protective ranking — reproduced in §3.12.

**Table 5. Supplementary tables.**

| | Contents |
|---|---|
| S1 | Sequence counts surviving each filter stage |
| S2 | Final accession list with isolate weight per unique sequence |
| S3 | Every removed sequence with the reason for removal |
| S4 | Per-gene cohort sizes |
| S5 | Software versions, generated from the running environment |

## 11. Figure placement summary

For the production editor. Placement is given relative to the section it supports.

| Figure | Place after | Why there |
|---|---|---|
| 1 | §3.1, first paragraph | The reader needs the pipeline before any count in the Results |
| 2 | §3.3, first paragraph | Carries the conservation figures that frame the power limitation |
| 3 | §3.4, second paragraph | Immediately after the BUSTED and aBSREL counts |
| **4** | §3.6, after the boxed statement | **The central figure; must share a spread with that sentence** |
| 5 | §3.7, second paragraph | After the HSV-1/HSV-2 split is stated, since the size encoding defends it |
| 6 | §3.9, after Table 3 | Supports both the class I panel and the class II bracket |
| 7 | §3.10, second paragraph | Immediately after the `SRVEAFHR` claim |
| S1 | cited in §3.1 | Supplementary; reference material, not argument |
| S2 | cited in §3.11 | Supplementary; supports the decision not to use either threshold |

---

## 12. Target journals

In order of fit, with the reason:

1. **Virus Evolution** — selection-inference methodology on a viral dataset; the concordance
   finding is squarely in scope and the journal publishes negative methodological results.
2. **Microbial Genomics** — comparable fit, receptive to negative results, and its Impact
   Statement format suits the methodological framing.
3. **BMC Genomics / PLOS Computational Biology** — broader audience, slower review.
4. **Scientific Data** — only for the curated dataset, as a companion resource paper.

Vaccine-focused venues (*Vaccine*, *npj Vaccines*) should be avoided without wet-lab
validation: the epitope predictions are unvalidated and reviewers there will say so. The
shortlist is offered here as a resource, and the manuscript leads with the methodological
result rather than the candidate list for exactly that reason.
