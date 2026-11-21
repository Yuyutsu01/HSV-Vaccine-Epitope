# Results

**Status:** draft for the manuscript, 2026-11-21. Every number here is read from a named
result table; none is transcribed by hand. Figures are in `20_figures/` and are regenerated
by `step27_figures.py` from the same tables. Methods are in `methods_draft.md`; the audit
trail, including nine documentation corrections and three substantive findings from the
verification pass, is in `CORRECTIONS_AND_PENDING_WORK.md`.

**One framing note before the numbers.** The headline of this study is a *negative* result
at site level, and it is stated as such. An earlier version of this analysis reported two
positively selected sites; §4.3 shows that both are artifacts of threshold and singleton
variants, and they are withdrawn as findings. The gene-level and branch-level results do
not depend on single alignment columns and stand.

---

## 1. Cohort

1,305 genome records were retrieved for HSV-1 (taxid 10298) and HSV-2 (taxid 10310). After
genome-level filtering, type verification by average nucleotide identity, per-gene CDS
extraction and QC, and collapse of exact duplicates, **4,163 unique coding sequences**
remain across **40 gene × type cohorts** for 20 genes, representing **12,054 isolates**
(`19_supplements/S1_data_flow.tsv`, `S4_per_gene_cohorts.tsv`; **Figure 1** for the
pipeline as a whole, **Figure S1** for cohort size per gene).

Cohort sizes are uneven and this limits what the HSV-2 arm can support. HSV-1 cohorts run
to 215 unique sequences (RS1); the smallest cohorts anywhere are UL49A at 28 (HSV-2) and 32
(HSV-1) and UL1 at 32 (HSV-2), and the two repeat-resident genes are thin in HSV-2 as well
(RS1 57, RL2 75). Two facts about the HSV-2 cohort must travel with every HSV-2 result:

- **178 of the 286 HSV-2 accessions (62%) come from a relaxed "tier 2" set** — genomes
  excluded at tier 1, overwhelmingly for 1–2% ambiguous bases, and readmitted because
  applying the ambiguity rule strictly leaves cohorts too small to analyse. The HSV-1
  cohorts contain none. HSV-2 inference therefore rests on assemblies of lower average
  quality than HSV-1 inference.
- Patient identity is not recorded in these records, so the cohorts are not samples of
  independent infections; serial isolates from one patient are present and are weighted,
  not removed (methods §9).

## 2. Alignment, masking and recombination

Rebuilding the alignments codon-aware was the single largest correction to the earlier
pipeline: within-type nucleotide identity for the indel-bearing genes rose from 64–82% to
97.6–99.5%, and masked lengths now land on the true reference protein lengths. Occupancy
masking dropped 4,037 codon columns across the 40 cohorts; the median cohort lost none and
the maximum loss was 63.4% (RS1/HSV-1), concentrated in the repeat-resident genes.

GARD found at least one breakpoint in **16 of 40 cohorts**, yielding **37 partitions**;
with the 24 unpartitioned cohorts this gives **61 analysis units**
(`14_selection/gard_summary_per_type.tsv`). Breakpoint ΔAICc ranges from 3.70 (UL44/HSV-1,
retained as the conservative choice) to 1,769.78 (RL2/HSV-1).

Importantly, **US6/gD shows no breakpoint in either type**. The earlier pooled-cohort
analysis had returned a recombination verdict for US6; running GARD within each type
removes it, confirming that the earlier verdict was the HSV-1/HSV-2 species split being
read as recombination. This is the direct evidence that the pooled-cohort results were
invalid and it is why they are withdrawn rather than updated.

One caveat is recorded rather than smoothed over: the RS1/HSV-1 GARD search is
**truncated**. It is the heaviest unit in the screen and was stopped by decision at its
confirmed two-breakpoint model after ~16 h and 5,215 models. The two breakpoints reported
are supported; a third is not excluded.

## 3. Conservation

These are strongly conserved coding sequences, and that is the central constraint on the
study's power. Mean within-type identity is **99.38% (HSV-1)** and **99.73% (HSV-2)**;
between-type identity ranges from **37.5% (US4/gG)** to **90.8% (UL30/pol)**
(`11_conservation/gene_level_conservation.tsv`; **Figure 2**). 281 windows of ≥9 aa are
conserved within both types and identical between them
(`11_conservation/conserved_windows.tsv`).

gG/US4 is confirmed as a poor cross-protective target: HSV-1 gG is 238 aa against 699 aa in
HSV-2 and occupies only 238 of 698 alignment columns, so its apparent "identity" is
dominated by non-overlap.

## 4. Selection

All **305 HyPhy runs** (61 units × 5 methods) completed, each behind a completion marker
(`14_selection/selection_runs.tsv`).

### 4.1 Gene-wide and branch-level signal is substantial

BUSTED detects gene-wide episodic diversifying selection in **32 of 61 analysis units** at
q < 0.05 after Benjamini–Hochberg correction across the family — equivalently **24 of the
40 gene × type cohorts** and **18 of the 20 genes**. Six genes (RS1, UL22, UL30, US4, US7,
US8) are significant in both types. aBSREL identifies **75 of 8,789 branches** under
episodic selection at Holm-corrected p < 0.05, concentrated in RS1/HSV-1 (13 branches),
RL2/HSV-1 (7), UL44/HSV-2 (6) and US7/HSV-1 (6) (**Figure 3**).

These are lineage- and gene-level signals. They do not depend on any single alignment
column and they are not affected by the site-level issues below.

### 4.2 FEL finds only purifying selection

FEL calls 154 sites at q < 0.05, and **every one has β < α**: pervasive selection in these
cohorts is uniformly purifying. No site-level positive selection is contributed by FEL at
all. For coding sequences at 99.4–99.7% identity this is the expected result and it frames
everything that follows.

### 4.3 The site-level result is negative

MEME and FUBAR together call **92 sites** as positively selected by at least one method
(87 FUBAR, 7 MEME), and **2 sites** satisfy the pre-registered ≥2-method concordance rule:
HSV-1 RL2/ICP0 residue 308 and HSV-1 US7/gI residue 155, both FUBAR+MEME.

**Neither survives a check on how many sequences actually carry the variant**
(`14_selection/site_robustness.tsv`; **Figure 4**):

| | RL2/ICP0, HSV-1, res. 308 | US7/gI, HSV-1, res. 155 |
|---|---|---|
| Majority codon | GGC (Gly) × 183 | CCA (Pro) × 125 |
| Variant codons | CGC × 1, TCC × 2 | GAC × 1, TCA × 1 |
| Isolates with a variant | **4 of 293 (1.4%)** | **2 of 276 (0.7%)** |
| Largest codon change | 2 nt | **3 nt** |
| FUBAR posterior | 0.9560 | 0.9523 |
| MEME q (β̂) | 0.0246 (1,024) | 0.0051 (1,842) |
| FEL | p = 0.169, n.s. | q = 0.637, n.s. |
| Supporting accessions | HM585511, OP431901, **OR771687 — fails CDS QC** | **PP379031, PP379033 — adjacent, one submission batch** |

Every indicator points the same way. Support is two to four isolates out of ~290. The FUBAR
posteriors sit a thousandth above the 0.95 cut-off. The MEME β̂ values of 1,024 and 1,842
are the unbounded estimates a near-invariant column produces, not biological rate
estimates. FEL supports neither. A three-nucleotide codon change observed in a single
isolate is far more consistent with a sequencing or assembly artifact than with a point
substitution. One site's only two variants come from consecutive accessions in a single
submission; the other's include a sequence that fails the CDS ambiguity criterion.

**Both sites are therefore reported as candidate artifacts and not as evidence of positive
selection.** `selected_sites_final.tsv` is left unchanged so the pre-registered rule's
output stands on the record, and `site_robustness.tsv` qualifies it.

The wider pattern inverts the concordance rule's suggestion. Of the 92 sites called by at
least one method, **76 are robust** (median 54 variant-carrying isolates against a median
cohort of 297), 7 are
weak and 9 are not robust. The robust ones are almost entirely FUBAR-only: well supported in
the data, but resting on a single method.

> **No site in this study is both multi-method concordant and robustly supported.**

This is the honest summary of the site-level screen. It is also a methodological result
worth reporting in its own right: requiring two methods to agree, intended as the
conservative choice, *selected for* artifacts here, because MEME fires on singleton
non-synonymous changes at invariant columns and FUBAR's posterior grazes threshold at the
same columns. Concordance and variant support are independent checks and both are needed.

## 5. Epitope prediction and its validation

Because NetMHCpan-4.1 requires an academic licence that could not be obtained
non-interactively, MHCflurry 2.3.30 was used as the free, already-trained equivalent.
Nothing was trained or tuned in this work. All 9- and 10-mers were enumerated from the type
reference proteins and scored against each of the 27 alleles of the class I reference panel
separately, giving **45,688 peptide predictions**
(`13_epitopes/class1_peptide_predictions.tsv`).

Validation is gene-grouped against 688 curated IEDB epitopes, so performance is not
inflated by one well-studied antigen. Over **27 scorable gene × type groups** spanning 17
genes, predicted affinity gives **mean AUC 0.809** (median 0.892, range 0.425–0.997), mean
average precision 0.326, and **6.0-fold enrichment** of known epitopes in the top decile
over 428 positives (`15_benchmark/benchmark_summary.tsv`; **Figure 5**).

The split by type is stark and is a limitation, not a nuisance: **HSV-1 groups give mean AUC
0.931, HSV-2 groups 0.600**. HSV-2 class I records in IEDB are fewer and are largely
inherited from HSV-1 cross-reactivity studies. Predictions for HSV-2 are correspondingly
less well validated.

The control that matters most: **conservation alone performs near chance** (mean AUC 0.566,
1.38-fold top-decile enrichment). Conservation is therefore not a proxy for immunogenicity,
which is the direct justification for using predicted binding and conservation as two
independent filters rather than either alone.

## 6. Candidate shortlist

Joining predicted binding, per-residue conservation, IEDB flags, glycosylation sequons and
the positive-selection filter gives 12,771 ranked candidates. The primary shortlist —
predicted binder, every residue ≥99.0% conserved within type, not positively selected — is
**9,488 rows / 8,010 distinct peptides**, of which **3,207 rows / 1,729 distinct peptides**
are identical across both types and form the cross-protective subset (`16_shortlist/`).
1,478 of those 1,729 have a row in each type and 251 in only one, because a peptide present
and identical in both types can still clear the within-type conservation threshold in only
one of them.

Cross-protective status is verified against the two reference proteins directly — the
peptide must occur contiguously in the other type's protein — rather than inferred from
per-residue identity through the alignment. The two are not the same claim: an indel inside
the window lets every residue correspond while the contiguous peptide exists in one type
only. The check removed 2 peptides (`16_shortlist/excluded_not_in_other_type.tsv`); before
the 2026-11-21 conservation-layer corrections it would have removed 444.

Two quality flags should travel with any use of this list:

- **1,460 of the 9,488 primary rows (15.4%) are low-complexity**, 1,307 distinct peptides
  drawn from the alanine-, serine-, glycine- and proline-rich tandem repeats of ICP4 and
  ICP0 and the gG repeat region — concentrated in RS1 (267 rows), RL2 (137) and US4 (113),
  peptides such as `AAAAAADLLF` and `AAAAPVVSA`. These are perfectly conserved and score as
  binders but are not credible epitopes. 6 appear in the top 50 by priority score, so the
  ranking suppresses them only partly; they are flagged rather than dropped, and barred
  from the minimal panel in §7, where it matters.
- The positive-selection exclusion removes only **5 rows**, and §4.3 shows the two sites
  driving it are probably artifacts. This filter is close to a no-op on these data and
  should not be presented as having shaped the shortlist.

As a positive control, the cross-protective shortlist recovers **`RMLGDVMAV`** — a published
HLA-A\*02 gB epitope — at **rank 13 of 3,207** by priority score, and rank 7 once
low-complexity repeats and duplicate rows are set aside, correctly flagged as IEDB-known
(gB 559–567 in HSV-2 numbering, 562–570 in HSV-1).

The top of the ranking is dominated by HLA-A\*02:03 and the A2 supertype. That is a known
property of class I predictors on hydrophobic 9-mers, not a biological finding, which is why
coverage is assessed per allele rather than from the top of the list.

## 7. Population coverage

Across **239 class I populations**, the cross-protective shortlist reaches **99.71% World
coverage**, with a median of 97.5% across populations and **73.6% of populations at or above
90%** (`18_population_coverage/`; **Figure 6**). A greedy minimal panel reaches 90.6% World
coverage with one peptide (`ILDNPGFVTF`), 99.2% with two, and **99.71% with five, covering
all 27 panel alleles** (`ILDNPGFVTF`, `EVIDKINAK`, `AMAAGRFGW`, `AQTGVRVAL`, `AAQYRNAVV`).

Low-complexity peptides are barred from the panel search — 1,307 are barred across the
candidate pool, leaving **1,559 of the 1,729** cross-protective candidates eligible —
because §6 treats them as non-credible epitopes and a panel offered as a design must not
rest on one. The exclusion costs nothing measurable: the unconstrained search also needs
five peptides and also reaches 27 alleles and 99.71%, differing only in using the
low-complexity `AAAVPLLSA`. The unconstrained panel is retained in
`minimal_panel_greedy_unconstrained.tsv` so the comparison is checkable.

The panel's earlier four-peptide version included `AEYGGSFRF` and `AAFDGGSTA`, both of
which the 2026-11-21 corrections removed: both were flagged cross-protective but occur in
HSV-2 only. The defect had reached the headline deliverable, which is the reason the
direct presence check above now exists.

Three caveats are reported with these figures:

1. **HLA-C is absent from the 27-allele panel**, so combined class I coverage is
   conservative.
2. **Coverage saturates.** Because the metric counts alleles with at least one binder, all
   three shortlist tiers give identical coverage once every panel allele is hit. Coverage
   discriminates *panels*, not shortlist depth; the per-population average hit counts are
   the quantity that distinguishes them.
3. **Low-coverage populations are an artifact of sparse typing**, not biology. Central
   America returns 7.76% because very few alleles are typed in that dataset. These are
   reported unadjusted and flagged.

## 8. Structural mapping

Candidates in the structurally characterised entry glycoproteins were mapped onto 2C36 (gD),
2GUM (gB) and 3M1C (gH–gL): **3,145 residues across 7 chains** assigned at 0.995–1.000
identity, of which 1,099 are solvent-exposed (RSA ≥ 0.25) and **32 conserved, exposed
surface patches** of ≥5 contiguous residues were identified (`17_structures/`).

The most promising is **gB 257–264 (`SRVEAFHR`)**, shown on the structure in
**Figure 7**: fully conserved within both types,
solvent-exposed, and annotated as inside a known IEDB epitope at all eight residues.

No positively selected site falls within any of the four structurally covered proteins, so
`structure_residue_annotation.tsv` reports `positively_selected = False` for all 3,145
residues. This is a limit of structural coverage, not evidence about those residues.

## 8a. B-cell epitope prediction

Score distributions for both B-cell predictors, and the behaviour of their published
thresholds on this data, are in **Figure S2**.

**BepiPred-3.0** was run on all 40 reference proteins (23,184 residues), and the per-residue
scores pass a strict register check against the reference translations before being used.
8,647 residues (37.3%) clear the published default threshold of 0.1512 — a proportion high
enough that the threshold is reported but not used; the continuous score enters the
candidate table as `bcell_linear_mean_score` (median 0.120 over 12,771 rows, range
0.018–0.323) as an annotation only. The highest-scoring genes are UL49 (VP22), UL1 (gL) and
UL54 (ICP27), all rich in low-complexity and disordered sequence, which is a known bias
direction for contact-trained B-cell predictors rather than a result.

**DiscoTope-3.0** was run on the same three structures as §8, scoring 7 chains and mapping
**3,145 of 3,145 residues** onto reference numbering with no unmapped rows and no
residue-letter mismatches. Its published high-confidence threshold of 0.90 is **never
reached**: scores span 0.0005–0.659 (median 0.071), and only 9 residues reach 0.50. Two
checks are worth reporting. Crystallographic replicate chains of the same protein, scored
independently, agree at *r* = 0.956 (gD, n = 274) and *r* = 0.886–0.904 (gB, n = 570),
which bounds how precise any single-structure score can be. And residues inside a curated
IEDB epitope score higher than the rest (median 0.089 vs 0.061, Mann–Whitney one-sided
*p* = 5.3 × 10⁻¹⁵) — genuinely out of sample, since DiscoTope has no access to IEDB, but a
small effect that is not a performance estimate.

The conformational layer covers 4 of 20 genes, because only those have a solved structure.
For the other 16 the result is missing, not negative.

## 8b. CD4 / class II binding and population coverage

**NetMHCIIpan-4.3** was run over the 40 reference proteins and the 27-allele IEDB class II
reference set: **610,848 predictions**, every one of which was verified to sit at its
reported position in the reference protein. 4,861 peptide–allele pairs (0.8%) are strong
binders at the server's own 1% threshold and 25,855 (4.2%) are weak or better — 0.8% being
close to the 1% a correctly calibrated %Rank should yield, which is a sanity check on the
rank column rather than a finding.

As with the B-cell arm, the summary statistic is not usable as a filter:
`classii_best_percentile_rank` is a minimum over 27 alleles × every overlapping 15-mer, so
it falls below 1% for **76%** of candidates. It is carried as an annotation; per-allele
binding is the quantity to use.

**Class II population coverage** over the 64 populations typed at all five presenting loci:
mean **98.87%** conservative (median 99.63%, minimum 90.75%, in the population the
IEDB file records as `Korea; South`) to 99.98% upper bound, and World **99.64–100%**. Combined with class I across the 51 populations in both
frequency sets, mean **99.98%**.

Three things qualify those numbers, and the second is the largest single correction made in
this analysis. (i) DP and DQ are heterodimers while the IEDB file stores chain frequencies,
so molecule coverage is **bracketed between a conservative and an upper-bound estimate**,
not computed. (ii) Populations differ in which loci were typed — only 64 of 274 have all
five — and treating an untyped locus as zero coverage puts the mean at **82.3% instead of
98.9%**, a 16-point artifact that would have read as a geographic finding. (iii) HLA-DRB3,
-DRB4 and -DRB5 have no frequency data in the file, so 4 of the 27 alleles cannot
contribute and the DR term is DRB1-only.

And as for class I, the figure describes the **allele panel**, not the shortlist: every
panel allele has at least one predicted binder among 22,624 peptides, so the weak and
strong sets return identical coverage. It is at the panel's ceiling by construction.

## 9. What this study does not establish

- **The CD4 arm is prediction, not evidence.** NetMHCIIpan-4.3 has now been run over the
  27-allele class II reference set, but no class II predictor is benchmarked here: unlike
  class I, no curated set of experimentally validated HSV class II epitopes was assembled
  to benchmark against. No peptide is claimed to be a CD4 epitope. *(Earlier drafts stated
  that this is a CD8-only study and that all three of NetMHCIIpan, BepiPred-3.0 and
  DiscoTope-3.0 were licence-blocked. Both statements are withdrawn: the two B-cell tools
  are free for academic use, and the class II predictions were obtained from the DTU web
  server.)*
- **The B-cell arm is prediction, not evidence.** Both predictors are run on the reference
  antigen, neither is benchmarked against experimentally mapped HSV B-cell epitopes here,
  and neither threshold behaves as published on this data — BepiPred's flags 37% of all
  residues, DiscoTope's flags none. The scores are a prioritisation signal and nothing more.
  No B-cell epitope is claimed.
- **No immunogenicity measurement.** Every candidate is predicted. The IEDB benchmark shows
  the predictor ranks known epitopes well in HSV-1 and poorly in HSV-2; it does not show
  that any shortlisted peptide is immunogenic. Absence of an IEDB record is never treated
  as a negative label, and no negative epitope set is constructed anywhere in this work.
- **No epidemiological inference.** The cohorts are convenience samples with unknown
  geographic composition and unknown patient structure; sampling-bias tables are reported
  (`03_metadata/sampling_bias_tables.tsv`) and no claim is made that rests on cohort
  frequencies being population frequencies.
- **No site-level positive selection.** Stated once more because it is the result most
  likely to be mis-cited: the screen's two concordant sites fail a minimum-variant-support
  check, and no site is both concordant and robust.

## Figures

All figures are produced by `09_scripts/step27_figures.py` from the result tables; no
value on any of them is typed in by hand. Each is written as 300 dpi PNG and as vector PDF.

| Figure | File | Shows |
|---|---|---|
| 1 | `20_figures/fig1_study_design.*` | Study design and data flow, including the two-track split at the selection stage and the withdrawn v1 branch |
| 2 | `20_figures/fig2_conservation.*` | Within-type and between-type identity across 20 genes |
| 3 | `20_figures/fig3_selection_gene_branch.*` | BUSTED gene-wide q-values for the 32 significant units, and aBSREL branch counts |
| 4 | `20_figures/fig4_site_variant_support.*` | Variant support behind every site called under selection |
| 5 | `20_figures/fig5_predictor_benchmark.*` | MHCflurry AUC per gene × type against curated IEDB epitopes, with marker area encoding the number of positives |
| 6 | `20_figures/fig6_population_coverage.*` | Class I coverage across 239 populations, the class II conservative-to-upper-bound bracket, and the greedy minimal panel |
| 7 | `20_figures/fig7_structure.*` | gB with the conserved exposed patch `SRVEAFHR`, and gD coloured by DiscoTope score |
| S1 | `20_figures/figS1_cohort_size.*` | Unique sequences per gene × type after deduplication |
| S2 | `20_figures/figS2_bcell_scores.*` | BepiPred and DiscoTope score distributions, and how their published thresholds behave here |

Figure 7 is drawn directly from the deposited coordinates by `step27_figures.py` as van der
Waals spheres in a principal-axis projection; it needs no external renderer. The PyMOL
session scripts in `17_structures/figures/` remain available for interactive inspection and
need PyMOL or ChimeraX.
