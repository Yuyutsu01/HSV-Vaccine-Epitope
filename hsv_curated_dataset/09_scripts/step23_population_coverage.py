"""
Step 23 (PENDING WORK B6, coverage arm): HLA population coverage of the epitope shortlist.

Proposal pipeline step 8 ("add HLA population coverage using the IEDB population coverage
tool"). The official IEDB Population Coverage package (v3.0.2) ships its allele frequency data
as a single pickle. Rather than execute the downloaded tool, this script reads ONLY that data
file, through a restricted unpickler that refuses every global except collections.OrderedDict,
and implements the published calculation (Bui et al., BMC Bioinformatics 2006) directly:

    coverage at a locus  = 1 - (1 - sum of frequencies of covered alleles)^2      [diploid, HWE]
    combined coverage    = 1 - product over loci of (1 - coverage at that locus)

Verified at import time: the data file contains no code objects, only nested OrderedDicts of
(allele, frequency) pairs.

Coverage is reported for four epitope sets, and a greedy search reports the smallest panel that
reaches a target coverage -- the practical question for a multi-epitope construct.

Note on the allele panel: step15 used the 27-allele HLA-A/HLA-B reference panel, so HLA-C
contributes nothing here and combined coverage is an underestimate of what a full A/B/C panel
would give. This is stated in the output.

Inputs:
  - 13_epitopes/class1_binders_by_allele.tsv   (step15)
  - 16_shortlist/shortlist_primary.tsv, shortlist_cross_protective.tsv (step21)

Outputs:
  - 18_population_coverage/allele_frequency_provenance.json
  - 18_population_coverage/coverage_by_population.tsv
  - 18_population_coverage/coverage_summary.tsv
  - 18_population_coverage/minimal_panel_greedy.tsv
  - 18_population_coverage/population_coverage.log
"""

import os
import sys
import json
import pickle
import hashlib
import collections
from datetime import datetime

import numpy as np
import pandas as pd

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_EPI = os.path.join(DIR_BASE, "13_epitopes")
DIR_SHORT = os.path.join(DIR_BASE, "16_shortlist")
DIR_OUT = os.path.join(DIR_BASE, "18_population_coverage")

SOURCE_URL = "https://downloads.iedb.org/tools/retrieve_package.php?tool=population"
SOURCE_TOOL = "IEDB Population Coverage 3.0.2"
DATA_RELPATH = ("population_coverage/population_coverage_pickle/population_genotype_map.p")

# Populations highlighted in the log; the full table covers all of them.
FOCUS = ["World", "Europe", "East Asia", "South Asia", "Sub-Saharan Africa",
         "North Africa", "North America", "South America", "Oceania",
         "Southeast Asia", "Southwest Asia", "Central America", "East Africa",
         "West Africa", "Central Africa"]

COVERAGE_TARGETS = (0.80, 0.90, 0.95)
GREEDY_MAX = 40

LOG_PATH = os.path.join(DIR_OUT, "population_coverage.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


class DataOnlyUnpickler(pickle.Unpickler):
    """Refuses every global except collections.OrderedDict, so no code can run on load."""

    def find_class(self, module, name):
        if (module, name) == ("collections", "OrderedDict"):
            return collections.OrderedDict
        raise pickle.UnpicklingError(
            f"refused to load {module}.{name} from the downloaded data file")


def load_frequencies(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    digest = hashlib.sha256(raw).hexdigest()
    with open(path, "rb") as fh:
        data = DataOnlyUnpickler(fh).load()
    return data, digest


def locus_of(allele):
    return allele.split("*")[0]


def coverage(covered_alleles, pop_data):
    """
    (combined coverage, {locus: coverage}) for one population.
    Only loci present in both the population data and the covered set contribute.
    """
    per_locus, not_covered = {}, 1.0
    for locus, pairs in pop_data.items():
        hit = sum(f for a, f in pairs if a in covered_alleles)
        hit = min(hit, 1.0)
        cov = 1.0 - (1.0 - hit) ** 2
        per_locus[locus] = cov
        not_covered *= (1.0 - cov)
    return 1.0 - not_covered, per_locus


def average_hits(peptide_alleles, pop_data):
    """Expected number of epitope-allele hits per individual (IEDB's average_hit)."""
    total = 0.0
    for alleles in peptide_alleles:
        for locus, pairs in pop_data.items():
            total += 2.0 * sum(f for a, f in pairs if a in alleles)
    return total


# --------------------------------------------------------------- class II coverage
#
# Class II needs its own treatment, and the reason is a property of the data rather than a
# choice. The IEDB frequency file stores class II as INDIVIDUAL CHAINS -- HLA-DRB1,
# HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1 -- but a DP or DQ molecule is a heterodimer, and
# NetMHCIIpan predicts for the PAIR. There is no frequency for a pair in this file, so a
# pair's coverage has to be derived from its two chains, and that can be done two ways:
#
#   upper bound   treat each chain locus as contributing independently, exactly as the
#                 class I calculation does. This OVERESTIMATES, because carrying a covered
#                 DQA1 chain does not mean the individual forms a covered DQ molecule.
#   conservative  require BOTH chains for the molecule: cov(DP) = cov(DPA1) * cov(DPB1).
#                 This UNDERESTIMATES, because it assumes the two chains are inherited
#                 independently, and DQA1-DQB1 in particular are in strong linkage
#                 disequilibrium, so real haplotypes pair them far more often than chance.
#
# Both are reported. The true value lies between them, and neither is presented as "the"
# class II coverage.
#
# Second data limitation: this file carries frequencies for HLA-DRB1 only. HLA-DRB3,
# -DRB4 and -DRB5 have no entry, so 4 of the 27 panel alleles cannot contribute at all and
# the DR component is DRB1-only. That makes both figures conservative in that respect, and
# it is why the numbers here will not reproduce IEDB's published >99% for this panel.
#
# HLA-DM is present in the file but is a peptide-loading chaperone, not a presenting
# molecule, so it is excluded.

CLASSII_PRESENTING_LOCI = ("HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1")
CLASSII_NO_FREQ_DATA = ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5")


def split_classii_chains(alleles):
    """IEDB-form class II molecule names -> {locus: set(chain alleles)}.

    'HLA-DRB1*01:01'              -> {'HLA-DRB1': {'HLA-DRB1*01:01'}}
    'HLA-DQA1*05:01/DQB1*02:01'   -> {'HLA-DQA1': {...}, 'HLA-DQB1': {...}}
    """
    chains = {}
    for a in alleles:
        parts = str(a).split("/")
        for i, part in enumerate(parts):
            part = part.strip()
            if i > 0 and not part.startswith("HLA-"):
                part = "HLA-" + part           # the beta chain is written without the prefix
            locus = part.split("*")[0]
            chains.setdefault(locus, set()).add(part)
    return chains


def coverage_class_ii(chains, pop_data):
    """(upper bound, conservative, {locus: coverage}, loci_present) for one population.

    Populations differ in which loci were typed. A population with no HLA-DRB1 entry is
    NOT a population with zero DR coverage -- it is a population that cannot be scored at
    DR, and treating the two alike produces false zeros. Several populations in this file
    carry only one class II locus; scoring them as 0% would drag the mean down by about
    eleven points and the error would look like a real geographic finding. The loci that
    were actually present are therefore returned, and the summary statistics are computed
    over the populations that can support the calculation.
    """
    per_locus = {}
    for locus in CLASSII_PRESENTING_LOCI:
        pairs = pop_data.get(locus)
        if not pairs:
            continue
        hit = min(sum(f for a, f in pairs if a in chains.get(locus, ())), 1.0)
        per_locus[locus] = 1.0 - (1.0 - hit) ** 2
    loci_present = tuple(sorted(per_locus))

    dr = per_locus.get("HLA-DRB1", 0.0)
    dp = per_locus.get("HLA-DPA1", 0.0) * per_locus.get("HLA-DPB1", 0.0)
    dq = per_locus.get("HLA-DQA1", 0.0) * per_locus.get("HLA-DQB1", 0.0)
    conservative = 1.0 - (1.0 - dr) * (1.0 - dp) * (1.0 - dq)

    not_cov = 1.0
    for locus in CLASSII_PRESENTING_LOCI:
        not_cov *= (1.0 - per_locus.get(locus, 0.0))
    return 1.0 - not_cov, conservative, per_locus, loci_present


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 23: HLA population coverage of the epitope shortlist")

    # Locate the downloaded data file.
    candidates = [
        os.path.join(DIR_OUT, "population_genotype_map.p"),
        os.path.join(os.environ.get("TEMP", ""), "claude"),
    ]
    data_path = None
    if os.path.exists(candidates[0]):
        data_path = candidates[0]
    else:
        for root, _dirs, files in os.walk(candidates[1]):
            if "population_genotype_map.p" in files:
                data_path = os.path.join(root, "population_genotype_map.p")
                break
    if not data_path:
        log("ABORT: population_genotype_map.p not found. Download the IEDB Population "
            f"Coverage package from {SOURCE_URL}, extract it, and place "
            f"{DATA_RELPATH} at {os.path.join(DIR_OUT, 'population_genotype_map.p')}")
        return 1

    freqs, digest = load_frequencies(data_path)
    log(f"Loaded allele frequency data (restricted unpickler, data only): "
        f"{len(freqs.get('I', {}))} class I populations, sha256 {digest[:16]}...")

    # Keep a copy inside the project so the analysis is reproducible.
    local = os.path.join(DIR_OUT, "population_genotype_map.p")
    if not os.path.exists(local):
        with open(data_path, "rb") as src, open(local, "wb") as dst:
            dst.write(src.read())
        log(f"Copied the data file into {os.path.relpath(local, DIR_BASE)}")

    with open(os.path.join(DIR_OUT, "allele_frequency_provenance.json"), "w",
              encoding="utf-8") as fh:
        json.dump(dict(tool=SOURCE_TOOL, source_url=SOURCE_URL,
                       data_file=DATA_RELPATH, sha256=digest,
                       retrieved=datetime.now().strftime("%Y-%m-%d"),
                       method="Bui et al., BMC Bioinformatics 2006; "
                              "reimplemented, tool code not executed",
                       note="class I only; step15 panel covers HLA-A and HLA-B, not HLA-C"),
                  fh, indent=2)

    class_i = freqs["I"]

    binders_path = os.path.join(DIR_EPI, "class1_binders_by_allele.tsv")
    if not os.path.exists(binders_path):
        log(f"ABORT: {binders_path} missing -- run step15 first")
        return 1
    binders = pd.read_csv(binders_path, sep="\t")
    log(f"Loaded {len(binders)} peptide-allele binder rows")

    pep_alleles = (binders.groupby("peptide").allele
                          .apply(lambda s: set(s.dropna().astype(str))).to_dict())

    sets = {"all_predicted_binders": set(pep_alleles)}
    low_complexity = set()
    for name, fname in (("shortlist_primary", "shortlist_primary.tsv"),
                        ("shortlist_cross_protective", "shortlist_cross_protective.tsv")):
        path = os.path.join(DIR_SHORT, fname)
        if os.path.exists(path):
            d = pd.read_csv(path, sep="\t")
            peps = set(d.peptide.unique())
            sets[name] = peps & set(pep_alleles)
            if "low_complexity" in d.columns:
                low_complexity |= set(d.loc[d.low_complexity.astype(bool), "peptide"].unique())
            log(f"{name}: {len(sets[name])} peptides with allele-level binder data")
        else:
            log(f"{name}: {fname} not found, skipped")

    # Low-complexity peptides are reported in the coverage sets unchanged, but they are
    # barred from the minimal panel. They come from the serine/glycine/alanine/proline
    # tandem repeats of ICP0 and ICP4: perfectly conserved, scored as binders, and not
    # credible epitopes -- and they cannot even be assigned a unique position. A panel
    # that a reader is meant to take seriously as a design must not rest on one.
    log(f"Low-complexity peptides barred from the minimal panel: {len(low_complexity)}")

    rows, summary = [], []
    for name, peps in sets.items():
        if not peps:
            continue
        covered = set().union(*(pep_alleles[p] for p in peps))
        per_pop = []
        for pop, pop_data in class_i.items():
            comb, per_locus = coverage(covered, pop_data)
            hits = average_hits([pep_alleles[p] for p in peps], pop_data)
            rows.append(dict(epitope_set=name, population=pop,
                             coverage=round(100 * comb, 2),
                             average_hits=round(hits, 2),
                             **{f"coverage_{l}": round(100 * c, 2)
                                for l, c in per_locus.items()}))
            per_pop.append(comb)
        arr = np.array(per_pop)
        world = next((r["coverage"] for r in rows
                      if r["epitope_set"] == name and r["population"] == "World"), None)
        summary.append(dict(epitope_set=name, n_peptides=len(peps),
                            n_alleles_covered=len(covered),
                            populations=len(arr),
                            world_coverage=world,
                            mean_coverage=round(100 * arr.mean(), 2),
                            median_coverage=round(100 * np.median(arr), 2),
                            min_coverage=round(100 * arr.min(), 2),
                            max_coverage=round(100 * arr.max(), 2),
                            pct_populations_above_90=round(
                                100 * float((arr >= 0.90).mean()), 1)))

    cov = pd.DataFrame(rows)
    cov.to_csv(os.path.join(DIR_OUT, "coverage_by_population.tsv"), sep="\t", index=False)
    summ = pd.DataFrame(summary)
    summ.to_csv(os.path.join(DIR_OUT, "coverage_summary.tsv"), sep="\t", index=False)
    log(f"Wrote coverage_by_population.tsv ({len(cov)} rows) and coverage_summary.tsv")

    for r in summ.to_dict("records"):
        log(f"  {r['epitope_set']:<28} {r['n_peptides']:>5} peptides, "
            f"{r['n_alleles_covered']:>2} alleles -> World {r['world_coverage']}%, "
            f"mean {r['mean_coverage']}% across {r['populations']} populations "
            f"(median {r['median_coverage']}%, min {r['min_coverage']}%, "
            f"{r['pct_populations_above_90']}% of populations >=90%)")

    log("Coverage in focus populations, cross-protective shortlist:")
    target_set = ("shortlist_cross_protective" if "shortlist_cross_protective" in sets
                  else "shortlist_primary")
    sub = cov[cov.epitope_set == target_set].set_index("population")
    for pop in FOCUS:
        if pop in sub.index:
            log(f"    {pop:<22} {sub.loc[pop, 'coverage']:>6.2f}%  "
                f"(avg hits {sub.loc[pop, 'average_hits']:.1f})")

    # Greedy minimal panel against World frequencies.
    if target_set in sets and "World" in class_i:
        world_data = class_i["World"]

        def greedy_panel(pool):
            """Greedy maximum-coverage search. Returns one row per peptide added.

            alleles_bound is how many panel alleles the peptide binds; new_alleles is
            how many of those were not already covered. They differ, and an earlier
            version reported only the first under the name 'alleles_added', which made
            the table look as though it failed to add up.
            """
            chosen, covered, out = [], set(), []
            best_cov = 0.0
            for step in range(1, min(GREEDY_MAX, len(pool)) + 1):
                best = None
                for p in pool:
                    if p in chosen:
                        continue
                    c, _ = coverage(covered | pep_alleles[p], world_data)
                    if best is None or c > best[1]:
                        best = (p, c)
                if best is None or best[1] <= best_cov + 1e-9:
                    break
                new = len(pep_alleles[best[0]] - covered)
                chosen.append(best[0])
                covered |= pep_alleles[best[0]]
                best_cov = best[1]
                out.append(dict(step=step, peptide=best[0],
                                alleles_bound=len(pep_alleles[best[0]]),
                                new_alleles=new,
                                cumulative_alleles=len(covered),
                                world_coverage=round(100 * best_cov, 2)))
            return out

        full_pool = sorted(sets[target_set])
        clean_pool = [p for p in full_pool if p not in low_complexity]

        greedy_rows = greedy_panel(clean_pool)
        pd.DataFrame(greedy_rows).to_csv(
            os.path.join(DIR_OUT, "minimal_panel_greedy.tsv"), sep="\t", index=False)
        log(f"Greedy minimal panel from the {target_set} set, low-complexity peptides "
            f"excluded ({len(clean_pool)} of {len(full_pool)} candidates eligible):")
        for r in greedy_rows[:12]:
            log(f"    +{r['peptide']:<12} -> {r['cumulative_alleles']:>2} alleles, "
                f"World coverage {r['world_coverage']}%")
        for t in COVERAGE_TARGETS:
            hit = next((r for r in greedy_rows if r["world_coverage"] >= 100 * t), None)
            log(f"    {int(100 * t)}% World coverage: "
                + (f"{hit['step']} peptides" if hit
                   else f"not reached with {len(greedy_rows)} peptides from this set"))

        # The unconstrained panel is kept for the record so the cost of the exclusion is
        # visible, not asserted. It is not the panel the manuscript reports.
        unc_rows = greedy_panel(full_pool)
        pd.DataFrame(unc_rows).to_csv(
            os.path.join(DIR_OUT, "minimal_panel_greedy_unconstrained.tsv"),
            sep="\t", index=False)
        if greedy_rows and unc_rows:
            lc_used = [r["peptide"] for r in unc_rows if r["peptide"] in low_complexity]
            log(f"Unconstrained panel for comparison: {len(unc_rows)} peptides, "
                f"{unc_rows[-1]['cumulative_alleles']} alleles, "
                f"World {unc_rows[-1]['world_coverage']}% "
                f"(low-complexity peptides used: {lc_used or 'none'})")
            log(f"Cost of excluding low-complexity peptides: "
                f"{round(unc_rows[-1]['world_coverage'] - greedy_rows[-1]['world_coverage'], 2)} "
                f"percentage points of World coverage, "
                f"{unc_rows[-1]['cumulative_alleles'] - greedy_rows[-1]['cumulative_alleles']} "
                f"panel alleles.")

    # ------------------------------------------------------------ class II coverage
    c2_path = os.path.join(DIR_EPI, "classii_binders_by_allele.tsv")
    if not os.path.exists(c2_path):
        log("Class II: classii_binders_by_allele.tsv not present; run step29. "
            "Class II coverage not computed.")
    elif "II" not in freqs:
        log("Class II: the frequency file has no class II block; coverage not computed.")
    else:
        log("")
        log("=" * 78)
        log("CLASS II POPULATION COVERAGE")
        c2 = pd.read_csv(c2_path, sep="\t")
        class_ii = freqs["II"]
        log(f"Loaded {len(c2)} class II peptide-allele binder rows, "
            f"{c2.allele.nunique()} alleles, {len(class_ii)} populations")

        c2_sets = {
            "classii_binders_weak_or_better": c2,
            "classii_strong_binders": c2[c2.percentile_rank <= 1.0],
        }

        world_chains = class_ii.get("World", {})
        all_alleles = sorted(c2.allele.unique())
        chains_all = split_classii_chains(all_alleles)
        log("panel decomposed into chains:")
        for locus in CLASSII_PRESENTING_LOCI + CLASSII_NO_FREQ_DATA:
            ours = chains_all.get(locus, set())
            have = set(a for a, _ in world_chains.get(locus, []))
            matched = ours & have
            note = ("  <-- NO FREQUENCY DATA for this locus in the IEDB file; "
                    "these alleles cannot contribute"
                    if locus in CLASSII_NO_FREQ_DATA else "")
            log(f"    {locus:10} panel {len(ours):>2}, matched to frequency data "
                f"{len(matched):>2}{note}")
            if ours and not note and (ours - have):
                log(f"        not found in the frequency data: {sorted(ours - have)}")

        c2_rows, c2_summary = [], []
        for name, sub in c2_sets.items():
            alleles = sorted(sub.allele.unique())
            chains = split_classii_chains(alleles)
            ub, cons = [], []
            for pop, pop_data in class_ii.items():
                u, c, per_locus, present = coverage_class_ii(chains, pop_data)
                c2_rows.append(dict(
                    epitope_set=name, population=pop,
                    coverage_upper_bound=round(100 * u, 2),
                    coverage_conservative=round(100 * c, 2),
                    n_loci_typed=len(present),
                    loci_typed=";".join(l.replace("HLA-", "") for l in present),
                    scorable=len(present) == len(CLASSII_PRESENTING_LOCI),
                    **{f"coverage_{l}": round(100 * v, 2)
                       for l, v in per_locus.items()}))
                ub.append(u); cons.append(c)

            # Summaries over the FULLY TYPED populations only. Including partially typed
            # ones would report their missing loci as absent coverage.
            full = [r for r in c2_rows
                    if r["epitope_set"] == name and r["scorable"]]
            fc = np.array([r["coverage_conservative"] for r in full]) / 100.0
            fu = np.array([r["coverage_upper_bound"] for r in full]) / 100.0
            w = next((r for r in c2_rows
                      if r["epitope_set"] == name and r["population"] == "World"), {})
            c2_summary.append(dict(
                epitope_set=name, n_peptides=int(sub.peptide.nunique()),
                n_alleles=len(alleles),
                populations_total=len(ub), populations_fully_typed=len(full),
                world_upper_bound=w.get("coverage_upper_bound"),
                world_conservative=w.get("coverage_conservative"),
                mean_upper_bound=round(100 * fu.mean(), 2),
                mean_conservative=round(100 * fc.mean(), 2),
                median_conservative=round(100 * float(np.median(fc)), 2),
                min_conservative=round(100 * fc.min(), 2),
                max_conservative=round(100 * fc.max(), 2),
                pct_fully_typed_pops_above_90=round(
                    100 * float((fc >= 0.90).mean()), 1)))

        pd.DataFrame(c2_rows).to_csv(
            os.path.join(DIR_OUT, "coverage_by_population_classii.tsv"),
            sep="\t", index=False)
        c2s = pd.DataFrame(c2_summary)
        c2s.to_csv(os.path.join(DIR_OUT, "coverage_summary_classii.tsv"),
                   sep="\t", index=False)
        log("")
        log("Wrote coverage_by_population_classii.tsv and coverage_summary_classii.tsv")
        for r in c2s.to_dict("records"):
            log(f"  {r['epitope_set']:<32} {r['n_peptides']:>6} peptides, "
                f"{r['n_alleles']} alleles")
            log(f"      World: {r['world_conservative']}% to {r['world_upper_bound']}% "
                f"(conservative to upper bound)")
            log(f"      across the {r['populations_fully_typed']} of "
                f"{r['populations_total']} populations typed at all 5 presenting loci: "
                f"mean {r['mean_conservative']}% to {r['mean_upper_bound']}%, "
                f"median {r['median_conservative']}%, min {r['min_conservative']}%, "
                f"{r['pct_fully_typed_pops_above_90']}% of them >=90%")

        log("")
        log("Class II coverage in focus populations (strong binders, "
            "conservative to upper bound):")
        c2cov = pd.DataFrame(c2_rows)
        c2cov = c2cov[c2cov.scorable]
        sub = c2cov[c2cov.epitope_set == "classii_strong_binders"].set_index("population")
        for pop in FOCUS:
            if pop in sub.index:
                log(f"    {pop:<22} {sub.loc[pop, 'coverage_conservative']:>6.2f}% to "
                    f"{sub.loc[pop, 'coverage_upper_bound']:>6.2f}%")

        # ---- combined class I + class II, per population, where both are available
        ci = cov[cov.epitope_set == target_set].set_index("population").coverage / 100.0
        comb_rows = []
        for pop in sub.index:
            if pop not in ci.index:
                continue
            a = float(ci.loc[pop])
            for label, col in (("conservative", "coverage_conservative"),
                               ("upper_bound", "coverage_upper_bound")):
                b = float(sub.loc[pop, col]) / 100.0
                comb_rows.append(dict(population=pop, class_ii_bound=label,
                                      class_i_coverage=round(100 * a, 2),
                                      class_ii_coverage=round(100 * b, 2),
                                      combined=round(100 * (1 - (1 - a) * (1 - b)), 2)))
        combined = pd.DataFrame(comb_rows)
        combined.to_csv(os.path.join(DIR_OUT, "coverage_combined_classI_classII.tsv"),
                        sep="\t", index=False)
        log("")
        log(f"Wrote coverage_combined_classI_classII.tsv "
            f"({combined.population.nunique()} populations present in BOTH the class I and "
            f"class II frequency sets)")
        for label in ("conservative", "upper_bound"):
            g = combined[combined.class_ii_bound == label]
            w = g[g.population == "World"]
            log(f"    class II {label:12}: World combined "
                f"{float(w.combined.iloc[0]) if len(w) else float('nan'):.2f}%, "
                f"mean {g.combined.mean():.2f}% across {len(g)} populations")

        log("")
        log("Class II caveats, in order of how much they matter:")
        log("  1. DP and DQ are heterodimers but the IEDB file stores chain frequencies, "
            "so molecule coverage is bracketed, not computed. Both bounds are reported.")
        log("  2. Populations differ in which class II loci were typed. A population with "
            "no HLA-DRB1 entry cannot be scored at DR, and is NOT a population with zero "
            "DR coverage; summary statistics therefore cover only the populations typed "
            "at all five presenting loci, and every row carries loci_typed and scorable.")
        log("  3. HLA-DRB3/4/5 have no frequency data in this file, so 4 of the 27 panel "
            "alleles cannot contribute and the DR term is DRB1-only.")
        log("  4. Coverage reflects the ALLELE PANEL, not epitope quality: with 22,624 "
            "peptides every allele has some predicted binder, so these figures are close "
            "to the panel's maximum and are not a property of the shortlist.")
        log("  5. No class II predictor is benchmarked in this study.")

    log("")
    log("Caveat: HLA-C is absent from the step15 panel, so combined class I coverage is "
        "conservative.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
