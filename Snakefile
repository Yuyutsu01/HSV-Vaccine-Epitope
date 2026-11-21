# Snakemake workflow for the HSV-1 / HSV-2 conservation, selection and epitope pipeline.
# Correction A10: the pipeline previously had no workflow definition, so its dependency
# order existed only in prose and could not be re-executed or checked.
#
# Each rule wraps one stage script. The scripts own their per-gene loops, so the rules
# are stage-level and use the one file per stage that cannot exist unless the stage
# finished, rather than sentinel touch-files.
#
#   snakemake -n                        show the plan without running anything
#   snakemake --cores 4 all             run everything that is out of date
#   snakemake --cores 4 conservation    stop after the conservation layer
#   snakemake --dag | dot -Tsvg > dag.svg
#
# Stages 1-7 (download, filter, verify, extract, QC, dedup) are retained from the
# original run and are declared here as inputs, not rules: re-downloading from NCBI
# would change the cohort and invalidate every published number. Their scripts remain
# in 09_scripts/ and the accession list in 19_supplements/S2_final_accessions.tsv
# pins exactly what was analysed.
#
# HyPhy stages need Docker Desktop running. GARD and the trees take hours.

import os

D = "hsv_curated_dataset"
S = f"{D}/09_scripts"
PY = "python"

# Stages whose inputs are frozen (see note above).
FROZEN = [
    f"{D}/07_qc_reports/08_dedup_log.tsv",
    f"{D}/05_cds_qc_passed",
]


rule all:
    input:
        f"{D}/18_population_coverage/coverage_summary.tsv",
        f"{D}/18_population_coverage/coverage_summary_classii.tsv",
        f"{D}/17_structures/structure_residue_annotation.tsv",
        f"{D}/16_shortlist/epitope_candidates_ranked.tsv",
        f"{D}/13_epitopes/bcell_conformational_scores.tsv",
        f"{D}/14_selection/selected_sites_final.tsv",
        f"{D}/14_selection/site_robustness.tsv",
        f"{D}/07_qc_reports/07d_cds_qc_log_complete.tsv",
        f"{D}/20_figures/fig4_site_variant_support.pdf",
        f"{D}/19_supplements/S1_data_flow.tsv",
        f"{D}/03_metadata/sampling_bias_tables.tsv",


# ---------------------------------------------------------------- alignment layer

rule codon_msa:
    """Correction A1/A2: true codon-aware profile MSA (FAMSA protein MSA + codon
    back-projection), replacing the reference-anchored pairwise star alignment."""
    input: FROZEN
    output:
        qc = f"{D}/06_alignments_v2/alignment_qc_v2.tsv",
    shell: f"{PY} {S}/step09b_rebuild_codon_msa.py"


rule mask_and_map:
    """Correction A3/A4: 50% occupancy masking over clinical sequences only, whole
    codons, plus alignment-column to reference-residue coordinate maps."""
    input: rules.codon_msa.output.qc
    output:
        report = f"{D}/06_alignments_v2/masking_report_v2.tsv",
    shell: f"{PY} {S}/step09c_mask_and_coordinate_map.py"


rule trees:
    """Per-type ML trees with 1000 ultrafast bootstraps. No combined-cohort trees:
    a pooled HSV-1+HSV-2 tree is a species tree, not a within-population genealogy."""
    input: rules.mask_and_map.output.report
    output:
        log = f"{D}/06_alignments_v2/trees_v2/trees.log",
    threads: 4
    shell: f"{PY} {S}/step10b_prep_and_trees.py"


# ------------------------------------------------------------- conservation layer

rule conservation:
    """Per-residue Shannon entropy and percent identity, isolate-weighted by the
    dedup cluster sizes, within and between virus types."""
    input: rules.mask_and_map.output.report
    output:
        per_residue = f"{D}/11_conservation/per_residue_conservation.tsv",
        gene_level  = f"{D}/11_conservation/gene_level_conservation.tsv",
        windows     = f"{D}/11_conservation/conserved_windows.tsv",
    shell: f"{PY} {S}/step13_conservation_layer.py"


# ----------------------------------------------------------------- selection arm

rule gard:
    """Correction A6/A7: GARD runs PER TYPE and BEFORE selection. The earlier pooled
    HSV-1+HSV-2 runs detected the species split, not recombination."""
    input: rules.trees.output.log
    output:
        summary = f"{D}/14_selection/gard_summary_per_type.tsv",
    threads: 4
    shell: f"{PY} {S}/step16_gard_per_type.py"


rule selection:
    """BUSTED, FEL, MEME, FUBAR, aBSREL per type, partitioned at any GARD breakpoint."""
    input: rules.gard.output.summary
    output:
        done = f"{D}/14_selection/selection_runs.tsv",
    threads: 4
    shell: f"{PY} {S}/step17_selection_per_type.py"


rule fdr:
    """Correction A5: Benjamini-Hochberg correction within gene x cohort x partition
    x method, cross-method concordance, and sites mapped to reference residues."""
    input:
        sel  = rules.selection.output.done,
        maps = rules.mask_and_map.output.report,
    output:
        final = f"{D}/14_selection/selected_sites_final.tsv",
    shell: f"{PY} {S}/step18_fdr_and_concordance.py"


rule site_robustness:
    """Verification pass V12: variant support behind every reported site. A q-value says
    nothing about how many sequences carry the variant, and on alignments this conserved a
    single non-synonymous singleton drives beta/alpha to infinity. Annotates, never
    re-filters: selected_sites_final.tsv is left as the pre-registered rule produced it."""
    input:
        final = rules.fdr.output.final,
        qc    = f"{D}/07_qc_reports/07d_cds_qc_log_complete.tsv",
    output:
        robust = f"{D}/14_selection/site_robustness.tsv",
    shell: f"{PY} {S}/step26_site_robustness.py"


# ------------------------------------------------------------------- epitope arm

rule iedb:
    """Curated experimentally validated epitopes. Records excluded by the curation
    rules are NOT negatives; absence of evidence is never a negative label."""
    output:
        curated = f"{D}/12_iedb/iedb_curated_epitopes.tsv",
        flags   = f"{D}/12_iedb/iedb_residue_flags.tsv",
    shell: f"{PY} {S}/step14_iedb_curation.py"


rule class1:
    """HLA class I prediction (MHCflurry 2.3) with conservation joined. Predicts per
    allele, not best-allele-only, so population coverage can be computed downstream."""
    input: rules.conservation.output.per_residue
    output:
        peptides = f"{D}/13_epitopes/class1_peptide_predictions.tsv",
        by_allele = f"{D}/13_epitopes/class1_binders_by_allele.tsv",
    shell: f"{PY} {S}/step15_tcell_class1_prediction.py"


rule benchmark:
    """Gene-grouped benchmark of the predictors against curated IEDB epitopes.
    Validation only; nothing here is trained or tuned."""
    input:
        pred = rules.class1.output.peptides,
        iedb = rules.iedb.output.curated,
    output:
        summary = f"{D}/15_benchmark/benchmark_summary.tsv",
    shell: f"{PY} {S}/step19_benchmark_predictors.py"


rule population_coverage:
    """IEDB population coverage (Bui et al. 2006), reimplemented from the published
    formula over the downloaded allele frequency data. Class I, plus class II bracketed
    between a conservative and an upper-bound estimate, because the frequency file stores
    class II per chain while DP and DQ are heterodimers."""
    input:
        rules.class1.output.by_allele,
        # Literal path, not rules.classii.output: `rule classii` is defined further down
        # the file and a rules.* reference only resolves to a rule already defined.
        # Snakemake matches producers by output path, so ordering does not matter here.
        f"{D}/13_epitopes/classii_binders_by_allele.tsv",
    output:
        summary = f"{D}/18_population_coverage/coverage_summary.tsv",
        classii = f"{D}/18_population_coverage/coverage_summary_classii.tsv",
    shell: f"{PY} {S}/step23_population_coverage.py"


# ------------------------------------------------------------ bias and integration

rule sampling_bias:
    """Patient assignment from serial-isolate naming and BioSample, one-per-patient
    cohorts, and the feasibility statement for the stratified comparisons."""
    input: FROZEN
    output:
        tables = f"{D}/03_metadata/sampling_bias_tables.tsv",
    shell: f"{PY} {S}/step20_patient_and_sampling_bias.py"


rule classii:
    """NetMHCIIpan-4.3 class II predictions, adapted to reference-residue coordinates.

    NOT part of `all`, and deliberately so: NetMHCIIpan is licence-gated and its output
    cannot be regenerated from this repository. The rule exists so that when the raw output
    is dropped into 13_epitopes/classii_raw/, the adapter and everything downstream of it
    run as a normal part of the workflow rather than as a manual step.

    Run `snakemake --cores 1 classii` after placing the files. step29 --inspect reports the
    detected layout without writing, which is what to run first on any new output file."""
    input:
        raw = f"{D}/13_epitopes/classii_raw",
    output:
        pred    = f"{D}/13_epitopes/classii_predictions.tsv",
        binders = f"{D}/13_epitopes/classii_binders_by_allele.tsv",
    shell: f'{PY} {S}/step29_classii_predictions.py --xls "{{input.raw}}/*.xls.gz"'


rule bcell_linear:
    """BepiPred-3.0 linear B-cell epitope scores, adapted to reference-residue coordinates.

    The predictor itself is NOT a rule. It is third-party code (bp3 0.0.12.7 + fair-esm,
    ESM-2 650M) that needs its own environment and a GPU, and its raw output is archived
    in 13_epitopes/bcell_raw/bepipred3/ so this rule is reproducible from the repository
    without re-running a 2.6 GB language model. step28 re-derives nothing: it checks that
    every scored residue letter matches the reference protein at the same position and
    aborts if the two are out of register."""
    input:
        raw = f"{D}/13_epitopes/bcell_raw/bepipred3/raw_output.csv",
    output:
        linear = f"{D}/13_epitopes/bcell_linear_scores.tsv",
        prov   = f"{D}/13_epitopes/bcell_predictor_provenance.tsv",
    shell: f"{PY} {S}/step28_bcell_predictors.py --bepipred {{input.raw}}"


rule integrate:
    """Candidate shortlist: conservation x predicted binding x known-epitope flags
    x glycosylation sequons. The selection filter is recorded as NOT_YET_APPLIED
    rather than silently passed if the selection arm has not produced results."""
    input:
        pred = rules.class1.output.peptides,
        cons = rules.conservation.output.per_residue,
        iedb = rules.iedb.output.flags,
        sel  = rules.fdr.output.final,
        bcl  = rules.bcell_linear.output.linear,
    output:
        ranked = f"{D}/16_shortlist/epitope_candidates_ranked.tsv",
    shell: f"{PY} {S}/step21_integrate_and_rank.py"


rule structures:
    """Surface exposure (Shrake-Rupley RSA) on 2C36 gD, 2GUM gB, 3M1C gH-gL, with
    chain assignment by sequence identity and PyMOL session scripts."""
    input: rules.integrate.output.ranked
    output:
        mapping = f"{D}/17_structures/structure_residue_annotation.tsv",
    shell: f"{PY} {S}/step22_structure_mapping.py"


rule bcell_conformational:
    """DiscoTope-3.0 conformational B-cell epitope propensity on the three solved
    structures, joined to reference residues through step22's coordinate map. Runs AFTER
    the structural layer because step22 owns the pdb_resnum -> ref_residue mapping; the
    linear scores run before integration because step21 consumes them.

    As with bcell_linear, the predictor is not a rule: its raw per-chain CSVs are archived
    under 13_epitopes/bcell_raw/discotope3/ and this rule only re-does the join, which it
    verifies residue-letter by residue-letter against structure_residue_annotation.tsv."""
    input:
        raw = f"{D}/13_epitopes/bcell_raw/discotope3",
        map = rules.structures.output.mapping,
    output:
        conf = f"{D}/13_epitopes/bcell_conformational_scores.tsv",
    shell: f"{PY} {S}/step28_bcell_predictors.py --discotope {{input.raw}}"


rule supplements:
    """Correction A10: data-flow counts, final accession list, and every removed
    sequence with its reason, generated from the QC logs rather than transcribed."""
    input:
        rules.mask_and_map.output.report,
        rules.sampling_bias.output.tables,
    output:
        flow = f"{D}/19_supplements/S1_data_flow.tsv",
    shell: f"{PY} {S}/step24_build_supplements.py"


rule figures:
    """Main-text figures, drawn from the result tables and never from transcribed values.
    Palette is the validated categorical set; the robustness verdict is encoded by marker
    shape because the natural status colours put green against red, which fails
    colour-vision separation."""
    input:
        rules.site_robustness.output.robust,
        f"{D}/16_shortlist/epitope_candidates_ranked.tsv",
        f"{D}/18_population_coverage/coverage_summary.tsv",
        f"{D}/11_conservation/gene_level_conservation.tsv",
    output:
        fig4 = f"{D}/20_figures/fig4_site_variant_support.pdf",
    shell: f"{PY} {S}/step27_figures.py"


rule qc_reconciliation:
    """Verification pass V2/V3/V10: one complete QC record over the CDS actually carried
    forward. 07_cds_qc_log.tsv covers 11,504 of the 12,054 sequences in 05_cds_qc_passed/,
    so before this rule no single file described the QC state of the analysed data.
    Read-only with respect to the cohort, which is pinned in S2_final_accessions.tsv."""
    output:
        complete = f"{D}/07_qc_reports/07d_cds_qc_log_complete.tsv",
        recon    = f"{D}/07_qc_reports/07e_qc_reconciliation.tsv",
    shell: f"{PY} {S}/step25_qc_reconciliation.py"


rule provenance_audit:
    """Completeness check, not a correctness check: every data file on disk must have
    an entry in 08_documentation/DATA_DICTIONARY.md. Exits non-zero and names the gaps
    if any file is undocumented, so a new output cannot enter the repository unexplained.
    Read-only -- writes nothing. Not in `rule all` because it audits the tree rather
    than producing a result, and because it must be able to fail without blocking a run."""
    shell: f"{PY} {S}/step30_provenance_audit.py"


# --------------------------------------------------------------- partial targets

rule alignments:
    input: rules.mask_and_map.output.report

rule conservation_only:
    input: rules.conservation.output.per_residue

rule epitopes:
    input:
        rules.benchmark.output.summary,
        rules.population_coverage.output.summary
