# Manuscript curation for publication

Assessment of `methods_draft.md` and `results_draft.md` against the result tables,
2026-11-21. Every number below was recomputed from the named file; nothing is taken
from the drafts on trust.

---

## 1. Verdict

**The science is in better shape than most submitted manuscripts.** The drafts are
internally consistent, every headline number reproduces exactly, and the
limitations section is unusually candid. The central result — that no site in the
study is both multi-method concordant and robustly supported — is stated as the
headline rather than buried, and the withdrawn earlier finding is documented
rather than quietly dropped.

What is missing is not analysis. It is **framing, one substantive fix, and
figures**. Specifically:

- one scientific inconsistency that a reviewer will find (§3.1);
- no study-design figure, which for this kind of paper is the figure readers look
  for first (§5.1);
- a CD4/class II arm with no figure at all (§5.2);
- a structural arm that exists only as unrendered PyMOL scripts (§5.3).

---

## 2. What was verified

All recomputed from the tables, independently of the drafts:

| Claim | Draft | Recomputed | |
|---|---|---|---|
| Records retrieved | 1,305 | 1,305 | OK |
| Unique CDS / cohorts / isolates | 4,163 / 40 / 12,054 | 4,163 / 40 / 12,054 | OK |
| HyPhy runs completed | 305 | 305 (61 units × 5 methods) | OK |
| BUSTED significant units | 32 of 61 | 32 of 61 | OK |
| FEL sites q<0.05 | 154 | 154 | OK |
| FEL sites with β<α | all | 154 of 154 | OK |
| FUBAR posterior >0.95 | 87 | 87 | OK |
| MEME q<0.05 | 7 | 7 | OK |
| Sites called by ≥1 method | 92 | 92 | OK |
| Concordant sites | 2 | 2 | OK |
| Mean AUC (affinity) | 0.809 | 0.8088 | OK |
| HSV-1 / HSV-2 AUC split | 0.931 / 0.600 | 0.931 / 0.600 | OK |
| Conservation-alone AUC | 0.565 | 0.5651 | OK |
| Within-type identity | 99.38 / 99.73% | 99.38 / 99.73% | OK |
| Between-type range | 37.5–90.8% | 37.5 (US4) – 90.8 (UL30) | OK |
| Conserved windows ≥9 aa | 283 | 283 | OK |
| Low-complexity shortlist rows | 1,527 / 9,450 (16.2%) | 1,527 / 9,450 (16.2%) | OK |
| Class I coverage / populations | 99.71% / 239 | 99.71% / 239 | OK |

**Note on the FDR scope.** BH correction is applied per
`(gene, cohort, partition, method)` group, not pooled across the family. Pooling
instead gives 95 FEL sites rather than 154. The per-group scope is what
`step18` does and what `fdr.log` records, and it is the defensible choice, but
the methods text should state the grouping explicitly — a reader who recomputes
naively will get a different number and assume an error.

---

## 3. Must fix before submission

### 3.1 The minimal panel contains a peptide the paper calls not credible

This is the one genuine scientific inconsistency.

Results §6 states that low-complexity peptides from the ICP0/ICP4 repeats "are
not credible epitopes." The greedy minimal panel in
`18_population_coverage/minimal_panel_greedy.tsv` is:

| Step | Peptide | Gene | Cumulative alleles | World coverage | `low_complexity` |
|---|---|---|---|---|---|
| 1 | `ILDNPGFVTF` | UL30 | 13/27 | 90.58% | False |
| 2 | `EVIDKINAK` | UL27 | 21/27 | 99.22% | False |
| 3 | `AEYGGSFRF` | UL27 | 25/27 | 99.62% | False |
| 4 | **`AAAAPLPSA`** | **RS1** | **27/27** | **99.71%** | **True** |

The headline 99.71% coverage depends on a poly-alanine peptide from ICP4 that the
manuscript elsewhere says should not be believed. Figure 6b shows it on the x-axis.

**FIXED, and it cost nothing.** `step23` now bars low-complexity peptides from
the greedy search (1,404 barred; 1,817 of 2,120 cross-protective candidates
remain eligible). The search simply finds a different fourth peptide:

| Step | Peptide | Gene | Cumulative alleles | World coverage |
|---|---|---|---|---|
| 1 | `ILDNPGFVTF` | UL30 | 13/27 | 90.58% |
| 2 | `EVIDKINAK` | UL27 | 21/27 | 99.22% |
| 3 | `AEYGGSFRF` | UL27 | 25/27 | 99.62% |
| 4 | **`AAFDGGSTA`** | **UL47** | **27/27** | **99.71%** |

`AAAAPLPSA` is replaced by `AAFDGGSTA` (UL47, `low_complexity=False`) at the
**same 27 alleles and the same 99.71% World coverage**. The cost of the
exclusion is 0.0 percentage points and 0 alleles.

An earlier version of this document predicted the fix would cost 0.09 points and
drop the panel to 3 peptides. That was a guess that the greedy search would
simply stop; it does not — it substitutes. The manuscript's headline coverage
number is unchanged, only the identity of the fourth peptide.

The unconstrained panel is still written to
`minimal_panel_greedy_unconstrained.tsv` so the comparison is on the record
rather than asserted.

### 3.2 ~~Supplementary table S1 has an unexplained 40-sequence gap~~ — WITHDRAWN

**This finding was wrong and is retracted.** An earlier version of this document
reported that S1's 4,203 → 4,163 step was unexplained and that methods §1
misstated the reference-anchor count. Both halves were incorrect:

- `S1_data_flow.tsv` **does** explain it, in the `detail` column: step 8 reads
  "40 alignments, includes one reference anchor each" and step 9 reads "40 gene ×
  type cohorts, reference anchor removed." The gap is documented precisely where
  it should be. The original check printed only `stage`, `level` and `n`, and
  never looked at `detail`.
- Methods §1's "one reference anchor per type" is about the **649 genomes
  entering per-gene extraction**, where the two reference genomes are indeed one
  per type. The 40 per-cohort anchors in the alignments are a different quantity
  and are documented separately.

The pipeline, S1 and methods §1 are consistent and correct. Verified directly:
`hyphy_ready/HSV1/RL2.fasta` carries 187 sequences and no `REF_` row, matching
S4's 187, and `step10b_prep_and_trees.py` drops the anchor before HyPhy with the
comment "these cohorts are clinical isolates."

Kept here rather than deleted because a retracted finding in a document about
checking findings should be visible.

### 3.3 `alleles_added` did not sum to `cumulative_alleles` — FIXED

In the published `minimal_panel_greedy.tsv`, step 3 added 7 to a cumulative 21
and reached 25, not 28, because the column counted alleles the peptide *binds*,
not alleles newly covered. Now split into `alleles_bound` and `new_alleles`;
`new_alleles` cumulatively sums to `cumulative_alleles` exactly.

---

## 4. Should fix

- **State the FDR grouping** in methods §8 (see §2 above).
- **Results §1 should carry the qualifier methods §1 has.** The results text says
  "178 of the 286 HSV-2 accessions (62%)" while S1 reports 346 analysis-ready
  HSV-2 genomes. Methods resolves this — 286 is the count *in the per-gene
  cohorts*, after CDS QC — but the results draft omits "in the per-gene cohorts"
  and the two numbers appear to contradict.
- **The 16-point class II coverage artifact deserves its own short methods
  subsection, not a parenthesis.** Treating untyped loci as zero coverage moved
  the mean from 82.3% to 98.9%. That is a reusable methodological warning for
  anyone using the IEDB population-coverage data, and it is currently buried in
  results §8b(ii). Reviewers reward this kind of thing when it is foregrounded.
- **Give the two withdrawn sites accession-level provenance in a supplement.**
  Results §4.3 names HM585511, OP431901, OR771687, PP379031, PP379033 in a table
  cell. Promote to a supplementary table with the codon, the QC status, and the
  submission batch, so the artifact argument is checkable.

---

## 5. Figures

### 5.1 The missing figure that matters most: study design

**There is no Figure 1 showing the pipeline.** For a multi-stage computational
study with a negative headline, this is the figure that determines whether a
reviewer understands the paper. It should show, as a single flow:

```
1,305 NCBI records
  → genome QC + ANI type verification         → 649 genomes (303 HSV-1, 117 + 229 HSV-2 tier-2)
  → per-gene CDS extraction + QC              → 12,054 CDS
  → exact deduplication                        → 4,163 unique, 40 gene × type cohorts
  → codon-aware MSA, masking, per-type trees
  → GARD                                       → 61 analysis units (37 partitions + 24 whole)
  → 5 HyPhy methods × 61 units                 → 305 runs
       ├─ gene/branch level  → 32 of 61 units, 75 of 8,789 branches     [POSITIVE]
       └─ site level         → 92 sites → 2 concordant → 0 robust       [NEGATIVE]
  → epitope prediction (class I, class II, B-cell) + IEDB benchmark
  → integration → 12,771 candidates → 9,450 primary → 3,202 cross-protective
```

Mark the tier-2 readmission and the withdrawn v1 branch explicitly. Make the
two-track split at the selection stage visually obvious — it is the paper's
argument in one picture.

### 5.2 Class II coverage has no figure

Results §8b reports 610,848 class II predictions and coverage of 98.87–99.98%,
and **no figure shows any of it.** Figure 6 is class I only. Either:

- extend Figure 6 to three panels (a: class I by population, b: class II
  bracketed conservative-to-upper-bound, c: minimal panel); or
- add a separate class II figure showing the conservative/upper-bound bracket per
  population, with the 64 scorable populations distinguished from the 210 that
  are not.

The bracket is the honest part of that analysis and deserves to be seen.

### 5.3 The structural arm is not rendered

`17_structures/figures/` contains four `.pml` scripts and no images. A
vaccine-epitope paper is expected to show the candidates on a structure. Render
at minimum one panel: gB with the conserved exposed patch **gB 257–264
(`SRVEAFHR`)** highlighted, which is the paper's strongest structural claim
(fully conserved in both types, solvent-exposed, inside a known IEDB epitope at
all eight residues). Add gD/2C36 with the DiscoTope score as a surface colour
ramp if a second panel is wanted.

### 5.4 Defects in the existing six

I viewed all six rendered PNGs.

| Fig | Verdict | Defects |
|---|---|---|
| 1 cohort size | Fine, but supplementary material | Legend box overlaps the "4,163 sequences…" caption |
| 2 conservation | Good | "gG/US4 is the poorest…" annotation collides with the `59.0` bar label. (I checked the apparently-missing RL2 HSV-1 point: HSV-1 99.30 and HSV-2 99.28 simply overlap. Not a defect.) |
| 3 selection | Good, dense | Panel (a) is 61 rows — at the edge of legibility at print width. Gene labels in (b) don't distinguish type except by colour (RS1 appears twice) |
| 4 site support | **Strongest figure in the set** | Lower panel has no y-axis label; a reviewer will ask what y encodes |
| 5 benchmark | Good | Legend overlaps the "mean AUC…" caption, obscuring "HSV-2 0.600". Duplicate gene labels (UL47, UL23, UL46, UL48, UL49, RS1, UL54, UL27 each appear twice) with no type suffix. No `n` per group — an AUC over 3 positives is drawn identically to one over 80 |
| 6 coverage | Good | Class I only (§5.2); panel b includes the low-complexity peptide (§3.1) |

**The legend-over-caption collision appears in Figures 1, 2 and 5** — it is a
systematic bug in `step27_figures.py`, not three coincidences. Fix once:
place the summary annotation with `ax.text` in figure coordinates below the axes,
or move the legend outside the axes with `bbox_to_anchor`.

### 5.5 Proposed final figure set

| # | Figure | Status |
|---|---|---|
| 1 | Study design and data flow | **To build** (§5.1) |
| 2 | Conservation within and between types | Exists, fix annotation collision |
| 3 | Selection: gene-wide and branch-level | Exists, trim panel (a) to significant units |
| 4 | Site-level support vs significance — **the headline** | Exists, add y-axis label |
| 5 | Predictor benchmark vs IEDB | Exists, disambiguate labels, encode n |
| 6 | Population coverage, class I **and class II** | Exists, extend (§5.2) |
| 7 | Structure: gB conserved exposed patch | **To build** (§5.3) |
| S1 | Per-gene cohort size (current Figure 1) | Demote to supplementary |
| S2 | B-cell score distributions, both predictors | **To build** — currently no figure for §8a |

That is 7 main figures, which is normal for this scope. Current Figure 1 moves to
supplementary because cohort size per gene is reference material, not argument.

---

## 6. Framing and target

The paper has three publishable contributions, and the order matters:

1. **A methodological result.** Requiring two methods to agree — the conservative
   choice — *selected for* artifacts, because MEME fires on singleton
   non-synonymous changes at near-invariant columns while FUBAR's posterior
   grazes threshold at the same columns. Concordance and variant support are
   independent checks and both are needed. This generalises well beyond HSV and
   is the most citable thing in the manuscript.
2. **A negative result, properly supported.** No site is both concordant and
   robust, in 305 HyPhy runs over 61 analysis units, at 99.4–99.7% sequence
   identity. Gene- and branch-level signal is real and unaffected.
3. **A curated resource.** 4,163 QC'd coding sequences, codon-aware per-type
   alignments, coordinate maps, and a fully documented pipeline.

**Lead with 1, support with 2, offer 3.** A draft that leads with the epitope
shortlist invites the reviewer question the paper cannot answer — is any of it
immunogenic? — whereas leading with the methodological result makes the shortlist
a resource rather than a claim.

Suggested venues, in order of fit:

- **Virus Evolution** — selection-inference methodology on a viral dataset; the
  concordance-rule finding is squarely in scope.
- **Microbial Genomics** — comparable fit, receptive to negative results.
- **BMC Genomics / PLOS Computational Biology** — broader, slower.
- **Scientific Data** — only for contribution 3, as a companion resource paper.

Avoid vaccine-focused venues (*Vaccine*, *npj Vaccines*) without wet-lab work;
the predictions are unvalidated and reviewers there will say so.

---

## 7. Before submission

- [x] Regenerate the minimal panel excluding low-complexity peptides (§3.1) —
      done; `AAFDGGSTA` replaces `AAAAPLPSA` at the same 27 alleles and 99.71%
- [x] ~~Fix the S1 reference-anchor gap~~ — withdrawn, the finding was wrong (§3.2)
- [x] Rename `alleles_added` → `alleles_bound` + `new_alleles` (§3.3)
- [x] Restore `population_genotype_map.p`, corrupted by CRLF conversion, and add
      `.gitattributes` so it cannot recur (see DATA_DICTIONARY §"Data integrity")
- [ ] State the per-group FDR scope in methods §8 (§2)
- [ ] Build Figures 1 and 7; extend Figure 6; build S2 (§5)
- [ ] Fix the legend/caption collision in `step27_figures.py` (§5.4)
- [ ] Add the y-axis label to Figure 4's lower panel (§5.4)
- [ ] Deposit the 587 MB of HyPhy JSON in Zenodo and cite the DOI in methods §19,
      rather than carrying it in the git repository
- [ ] Confirm the GitHub repository is public and the commit referenced in §18 is
      the one pushed

**On the last two.** Methods §19 should point at an archived, versioned DOI, not
a GitHub URL — journals increasingly require it, and a repository can change after
review. Depositing the bulk HyPhy output in Zenodo also resolves the practical
problem that the 637 MB commit currently fails to push to GitHub with an HTTP 408
timeout: the code and derived tables belong in git, the raw model fits belong in a
data archive with a DOI.
