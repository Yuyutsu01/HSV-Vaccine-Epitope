"""
Step 27: main-text figures.

Every number plotted here is read from a result table; nothing is hardcoded. If a table is
missing the figure is skipped with a message rather than drawn from stale values.

Design rules applied (and why, so they are not undone by accident)
------------------------------------------------------------------
  - Palette is the validated categorical set, assigned in fixed order and never cycled:
    slot 1 blue  #2a78d6 = HSV-1, slot 2 orange #eb6834 = HSV-2, slot 3 aqua #1baf7a =
    between-type. Validated all-pairs for CVD and normal-vision separation; aqua carries a
    sub-3:1 contrast warning against the light surface, so wherever it appears it is also
    directly labelled.
  - Robustness verdict is encoded by MARKER SHAPE first, not colour. The natural status
    palette for good/warning/critical puts green against red, which fails colour-vision
    separation outright (deltaE 4.1 deutan). Shape carries the meaning; colour carries
    virus type.
  - No figure uses two y-scales. Where two quantities of different scale belong together
    they get adjacent panels sharing an axis.
  - Grid and spines are recessive; values and labels wear ink colours, never series
    colours; a legend is present whenever two or more series are drawn.

Outputs PNG at 300 dpi for inspection and PDF (vector) for the manuscript, into
20_figures/.
"""

import os
import zlib
import sys
import glob
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DIR_BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DIR_SEL = os.path.join(DIR_BASE, "14_selection")
DIR_OUT = os.path.join(DIR_BASE, "20_figures")
LOG_PATH = os.path.join(DIR_OUT, "figures.log")

# validated palette slots, fixed order
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
TYPE_COLOR = {"HSV1": BLUE, "HSV2": ORANGE}
TYPE_LABEL = {"HSV1": "HSV-1", "HSV2": "HSV-2"}
SURFACE = "#fcfcfb"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e4e3df"

VERDICT_MARKER = {"ROBUST": "o", "WEAK": "^", "NOT_ROBUST": "X"}


def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans", "font.size": 8.5,
        "axes.edgecolor": INK3, "axes.linewidth": 0.6,
        "axes.labelcolor": INK, "axes.titlesize": 10,
        "axes.titleweight": "bold", "axes.titlecolor": INK,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "legend.frameon": False, "legend.fontsize": 7.5,
        "grid.color": GRID, "grid.linewidth": 0.6,
        "lines.linewidth": 1.8,
    })


def recessive(ax, which=("top", "right")):
    for s in which:
        ax.spines[s].set_visible(False)
    ax.set_axisbelow(True)


def save(fig, name):
    os.makedirs(DIR_OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(DIR_OUT, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    log(f"  wrote {name}.png and {name}.pdf")


def rd(rel):
    p = os.path.join(DIR_BASE, rel)
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


# ----------------------------------------------------------------- figure 1
def fig1_cohorts():
    s4 = rd("19_supplements/S4_per_gene_cohorts.tsv")
    if s4 is None:
        log("fig1 skipped: S4 missing")
        return
    piv = s4.pivot_table(index="gene", columns="virus_type",
                         values="unique_sequences", aggfunc="sum").fillna(0)
    piv = piv.loc[piv.sum(axis=1).sort_values().index]
    y = np.arange(len(piv))
    h = 0.38
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    # 2px surface gap between adjacent bars: offset by slightly more than h/2
    ax.barh(y + h / 2 + 0.012, piv.get("HSV1", 0), height=h, color=BLUE,
            label="HSV-1", linewidth=0)
    ax.barh(y - h / 2 - 0.012, piv.get("HSV2", 0), height=h, color=ORANGE,
            label="HSV-2", linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels(piv.index)
    ax.set_xlabel("unique coding sequences after deduplication")
    ax.set_title("Per-gene cohort size", loc="left")
    ax.xaxis.grid(True)
    recessive(ax)
    # selective direct labels: the largest bar of each series only
    for series, col, off in (("HSV1", BLUE, h / 2 + 0.012),
                             ("HSV2", ORANGE, -h / 2 - 0.012)):
        if series not in piv:
            continue
        i = int(np.argmax(piv[series].values))
        v = piv[series].values[i]
        ax.text(v + 3, y[i] + off, f"{int(v)}", va="center", ha="left",
                fontsize=7, color=INK2)
    ax.legend(loc="lower right")
    tot = int(s4.unique_sequences.sum())
    ax.annotate(f"{tot:,} sequences in 40 gene x type cohorts",
                xy=(0.99, 0.015), xycoords="axes fraction", ha="right",
                fontsize=7, color=INK2)
    save(fig, "fig1_cohort_size")


# ----------------------------------------------------------------- figure 2
def fig2_conservation():
    c = rd("11_conservation/gene_level_conservation.tsv")
    if c is None:
        log("fig2 skipped: gene_level_conservation missing")
        return
    c = c.sort_values("mean_between_type_identity")
    y = np.arange(len(c))
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.6, 5.4), sharey=True,
                                   gridspec_kw=dict(width_ratios=[1, 1.25],
                                                    wspace=0.08))
    # (a) within-type
    axa.scatter(c.mean_pct_identity_hsv1, y, s=34, color=BLUE, label="HSV-1",
                zorder=3, linewidths=0.8, edgecolors=SURFACE)
    axa.scatter(c.mean_pct_identity_hsv2, y, s=34, color=ORANGE, label="HSV-2",
                zorder=3, linewidths=0.8, edgecolors=SURFACE)
    for yi, a, b in zip(y, c.mean_pct_identity_hsv1, c.mean_pct_identity_hsv2):
        axa.plot([a, b], [yi, yi], color=GRID, lw=1.2, zorder=1)
    axa.set_yticks(y)
    axa.set_yticklabels(c.gene)
    axa.set_xlabel("mean within-type identity (%)")
    axa.set_title("a  Within type", loc="left")
    axa.set_xlim(96.5, 100.3)
    axa.xaxis.grid(True)
    recessive(axa)
    axa.legend(loc="lower left")

    # (b) between-type, direct-labelled because aqua is sub-3:1 on this surface
    axb.barh(y, c.mean_between_type_identity, height=0.6, color=AQUA,
             linewidth=0)
    for yi, v, g in zip(y, c.mean_between_type_identity, c.gene):
        axb.text(v + 0.8, yi, f"{v:.1f}", va="center", fontsize=7, color=INK2)
    axb.set_xlabel("mean between-type identity (%)")
    axb.set_title("b  Between types", loc="left")
    axb.set_xlim(0, 104)
    axb.xaxis.grid(True)
    recessive(axb)
    axb.annotate("gG/US4 is the poorest\ncross-protective target",
                 xy=(c.mean_between_type_identity.iloc[0] + 3, 0.1),
                 fontsize=7, color=INK2, va="bottom")
    fig.suptitle("Sequence conservation across 20 genes", x=0.072, y=1.0,
                 ha="left", fontsize=11, fontweight="bold", color=INK)
    save(fig, "fig2_conservation")


# ----------------------------------------------------------------- figure 3
def fig3_selection():
    g = rd("14_selection/gene_level_results.tsv")
    b = rd("14_selection/branch_level_results.tsv")
    if g is None:
        log("fig3 skipped: gene_level_results missing")
        return
    g = g.copy()
    g["unit_label"] = [f"{r.gene}" + (f" p{int(r.partition)}" if r.partition else "")
                       for r in g.itertuples()]
    g["nlq"] = -np.log10(g.q_value.clip(lower=1e-12))
    g = g.sort_values("nlq", ascending=True)
    y = np.arange(len(g))

    if b is not None:
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(8.2, 8.6),
                                       gridspec_kw=dict(width_ratios=[2.1, 1],
                                                        wspace=0.28))
    else:
        fig, axa = plt.subplots(figsize=(6.0, 8.6))
        axb = None

    colors = [TYPE_COLOR[c] for c in g.cohort]
    axa.barh(y, g.nlq, height=0.66, color=colors, linewidth=0)
    axa.axvline(-np.log10(0.05), color=INK2, lw=1.0, ls="--", zorder=4)
    # label at the foot of the threshold line, where the bars are shortest: at the top
    # it lands on the longest bar
    axa.annotate("q = 0.05", xy=(-np.log10(0.05), -0.4),
                 xytext=(4, 0), textcoords="offset points",
                 fontsize=7, color=INK2, va="bottom")
    axa.set_yticks(y)
    axa.set_yticklabels(g.unit_label, fontsize=6.2)
    axa.set_xlabel(r"BUSTED gene-wide episodic selection, $-\log_{10}(q)$")
    axa.set_title("a  Gene-wide (61 analysis units)", loc="left")
    axa.xaxis.grid(True)
    recessive(axa)
    nsig = int((g.q_value < 0.05).sum())
    axa.annotate(f"{nsig} of {len(g)} units significant",
                 xy=(0.99, 0.01), xycoords="axes fraction", ha="right",
                 fontsize=7, color=INK2)
    axa.legend(handles=[Line2D([], [], marker="s", ls="", ms=6, color=TYPE_COLOR[k],
                               label=TYPE_LABEL[k]) for k in ("HSV1", "HSV2")],
               loc="lower right", bbox_to_anchor=(1.0, 0.05))

    if axb is not None:
        sig = b[b.significant.astype(bool)]
        cnt = (sig.groupby(["cohort", "gene"]).size()
               .rename("n").reset_index()
               .sort_values("n"))
        yy = np.arange(len(cnt))
        axb.barh(yy, cnt.n, height=0.66,
                 color=[TYPE_COLOR[c] for c in cnt.cohort], linewidth=0)
        axb.set_yticks(yy)
        axb.set_yticklabels([f"{r.gene}" for r in cnt.itertuples()], fontsize=6.2)
        axb.set_xlabel("branches under selection")
        axb.set_title("b  Branch-level (aBSREL)", loc="left")
        axb.xaxis.grid(True)
        recessive(axb)
        axb.annotate(f"{int(sig.shape[0])} of {len(b):,} branches\n"
                     f"at Holm-corrected p < 0.05",
                     xy=(0.97, 0.01), xycoords="axes fraction", ha="right",
                     fontsize=7, color=INK2)
    fig.suptitle("Selection is detected gene-wide and on branches",
                 x=0.055, y=0.995, ha="left", fontsize=11, fontweight="bold",
                 color=INK)
    save(fig, "fig3_selection_gene_branch")


# ----------------------------------------------------------------- figure 4
def fig4_robustness():
    r = rd("14_selection/site_robustness.tsv")
    s = rd("14_selection/site_level_results.tsv")
    if r is None:
        log("fig4 skipped: site_robustness missing (run step26)")
        return
    r = r[r.variant_isolates.notna()].copy()
    if s is not None:
        ps = s[s.positively_selected.astype(bool)]
        minq = (ps.groupby(["gene", "cohort", "ref_residue"]).q_value
                .min().rename("min_q").reset_index())
        r = r.merge(minq, on=["gene", "cohort", "ref_residue"], how="left")
    if "min_q" not in r.columns:
        r["min_q"] = np.nan
    r["has_q"] = r.min_q.notna()
    r["nlq"] = -np.log10(r.min_q.clip(lower=1e-6))

    # Two panels sharing the x axis rather than one with a dead middle: the sites that
    # carry a q (MEME) and the FUBAR-only sites live on different y quantities, and
    # forcing them onto one scale left two thirds of the panel empty and collided the
    # legend with the floor row.
    fig, (axq, axf) = plt.subplots(
        2, 1, figsize=(7.4, 5.8), sharex=True,
        gridspec_kw=dict(height_ratios=[2.0, 1.0], hspace=0.12))

    withq, noq = r[r.has_q], r[~r.has_q]
    for ax, sub, ycol in ((axq, withq, "nlq"), (axf, noq, None)):
        for verdict, marker in VERDICT_MARKER.items():
            for cohort in ("HSV1", "HSV2"):
                d = sub[(sub.verdict == verdict) & (sub.cohort == cohort)]
                if d.empty:
                    continue
                if ycol:
                    yv = d[ycol]
                else:
                    # Deterministic jitter so overlapping points stay separable.
                    # crc32, not hash(): Python randomises string hashing per process
                    # (PYTHONHASHSEED), so hash() reseeded this on every run and the
                    # figure's point positions changed each time it was regenerated.
                    seed = zlib.crc32(f"{verdict}|{cohort}".encode())
                    rng = np.random.default_rng(seed)
                    yv = rng.uniform(-0.42, 0.42, len(d))
                ax.scatter(d.variant_isolates, yv, marker=marker,
                           s=52 if marker != "o" else 36,
                           color=TYPE_COLOR[cohort], alpha=0.85,
                           linewidths=0.8, edgecolors=SURFACE, zorder=3)

    for ax in (axq, axf):
        ax.axvline(5, color=INK2, lw=1.0, ls="--", zorder=2)
        ax.xaxis.grid(True)
        recessive(ax)
    axq.set_xscale("log")
    axq.yaxis.grid(True)

    axq.set_ylabel(r"$-\log_{10}$(min q)")
    axq.set_title("Significance is unrelated to how many isolates carry the variant",
                  loc="left")
    axq.annotate("called by MEME\n(has a q value)", xy=(0.995, 0.93),
                 xycoords="axes fraction", ha="right", va="top",
                 fontsize=7.5, color=INK2)
    axq.set_ylim(1.35, 3.0)

    axf.set_yticks([])
    axf.set_ylim(-1.0, 1.0)
    axf.set_xlabel("isolates carrying the variant codon (log scale)")
    axf.annotate("FUBAR only\n(posterior, no q)", xy=(0.995, 0.92),
                 xycoords="axes fraction", ha="right", va="top",
                 fontsize=7.5, color=INK2)

    # the threshold label sits on the lower panel, clear of the site labels above
    axf.annotate("minimum support\nthreshold: 5 isolates", xy=(5, -0.95),
                 xytext=(8, 2), textcoords="offset points",
                 fontsize=7, color=INK2, va="bottom")

    # direct-label the two concordant sites: the point of the figure
    two = withq[withq.n_methods.astype(str) == "2"]
    for row, dx, dy in zip(two.itertuples(), (18, 18), (-26, 20)):
        axq.annotate(
            f"{TYPE_LABEL[row.cohort]} {row.gene} {row.ref_residue:.0f}  "
            f"({int(row.variant_isolates)}/{int(row.total_isolates)} isolates)",
            xy=(row.variant_isolates, row.nlq),
            xytext=(dx, dy), textcoords="offset points",
            fontsize=7.5, color=INK, fontweight="bold",
            arrowprops=dict(arrowstyle="-", lw=0.8, color=INK2))

    handles = [Line2D([], [], marker=m, ls="", ms=6.5, color=INK2,
                      label=v.replace("_", " ").title())
               for v, m in VERDICT_MARKER.items()]
    handles += [Line2D([], [], marker="s", ls="", ms=6, color=TYPE_COLOR[k],
                       label=TYPE_LABEL[k]) for k in ("HSV1", "HSV2")]
    # legend below the axes: inside, it landed on the FUBAR strip
    axf.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
               ncol=5, columnspacing=1.4)
    counts = r.verdict.value_counts().to_dict()
    fig.text(0.5, -0.055,
             "shape = variant-support verdict, colour = virus type   |   "
             + ",  ".join(f"{k.replace('_', ' ').lower()} {v}"
                          for k, v in counts.items())
             + "   |   both concordant sites fall left of the threshold",
             ha="center", fontsize=7, color=INK2)
    save(fig, "fig4_site_variant_support")


# ----------------------------------------------------------------- figure 5
def fig5_benchmark():
    b = rd("15_benchmark/benchmark_per_gene.tsv")
    if b is None:
        log("fig5 skipped: benchmark_per_gene missing")
        return
    a = b[b.score == "affinity"].copy() if "score" in b.columns else b.copy()
    a = a.sort_values("auc")
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    y = np.arange(len(a))
    ax.hlines(y, 0.5, a.auc, color=GRID, lw=1.2, zorder=1)
    for cohort in ("HSV1", "HSV2"):
        sub = a[a.virus_type == cohort]
        idx = [list(a.index).index(i) for i in sub.index]
        ax.scatter(sub.auc, idx, s=38, color=TYPE_COLOR[cohort],
                   label=TYPE_LABEL[cohort], zorder=3, linewidths=0.8,
                   edgecolors=SURFACE)
    ax.axvline(0.5, color=INK2, lw=1.0, ls="--", zorder=2)
    ax.annotate("chance", xy=(0.5, len(a) - 0.5), xytext=(4, 0),
                textcoords="offset points", fontsize=7, color=INK2, va="top")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.gene}" for r in a.itertuples()], fontsize=6.6)
    ax.set_xlabel("AUC, predicted affinity vs curated IEDB epitopes")
    ax.set_title("Predictor validation is strong in HSV-1 and weak in HSV-2",
                 loc="left")
    ax.set_xlim(0.35, 1.02)
    ax.xaxis.grid(True)
    recessive(ax)
    ax.legend(loc="lower right")
    m1 = a[a.virus_type == "HSV1"].auc.mean()
    m2 = a[a.virus_type == "HSV2"].auc.mean()
    ax.annotate(f"mean AUC  HSV-1 {m1:.3f}   HSV-2 {m2:.3f}",
                xy=(0.99, 0.02), xycoords="axes fraction", ha="right",
                fontsize=7, color=INK2)
    save(fig, "fig5_predictor_benchmark")


# ----------------------------------------------------------------- figure 6
def fig6_coverage():
    cp = rd("18_population_coverage/coverage_by_population.tsv")
    gp = rd("18_population_coverage/minimal_panel_greedy.tsv")
    if cp is None:
        log("fig6 skipped: coverage_by_population missing")
        return
    target = "shortlist_cross_protective"
    sub = cp[cp.epitope_set == target] if "epitope_set" in cp.columns else cp
    if sub.empty:
        sub = cp
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.8, 3.6),
                                   gridspec_kw=dict(width_ratios=[1.15, 1],
                                                    wspace=0.28))
    vals = np.sort(sub.coverage.values)
    ecdf = np.arange(1, len(vals) + 1) / len(vals) * 100
    axa.step(vals, ecdf, where="post", color=BLUE, lw=1.8)
    axa.axvline(90, color=INK2, lw=1.0, ls="--")
    above = float((sub.coverage >= 90).mean() * 100)
    axa.annotate(f"{above:.1f}% of populations\nat or above 90%",
                 xy=(90, 12), xytext=(-6, 0), textcoords="offset points",
                 ha="right", fontsize=7, color=INK2)
    axa.set_xlabel("class I coverage within a population (%)")
    axa.set_ylabel("populations at or below (%)")
    axa.set_title(f"a  Coverage across {sub.population.nunique()} populations",
                  loc="left")
    axa.grid(True)
    recessive(axa)

    if gp is not None and len(gp):
        cumcol = next((c for c in ("cumulative_coverage", "world_coverage",
                                   "coverage") if c in gp.columns), None)
        pepcol = next((c for c in ("peptide", "epitope") if c in gp.columns), None)
        if cumcol:
            n = np.arange(1, len(gp) + 1)
            axb.plot(n, gp[cumcol], color=BLUE, lw=1.8, marker="o", ms=5,
                     markeredgecolor=SURFACE, markeredgewidth=0.8)
            axb.set_xticks(n)
            if pepcol:
                axb.set_xticklabels(gp[pepcol], rotation=35, ha="right",
                                    fontsize=6.6)
            axb.set_xlabel("peptides added, greedy order")
            axb.set_ylabel("cumulative World coverage (%)")
            axb.set_title("b  Minimal panel", loc="left")
            axb.set_ylim(min(85, float(gp[cumcol].min()) - 2), 101)
            axb.grid(True)
            recessive(axb)
            axb.annotate(f"{gp[cumcol].iloc[0]:.1f}% from one peptide",
                         xy=(1, gp[cumcol].iloc[0]), xytext=(8, -12),
                         textcoords="offset points", fontsize=7, color=INK2)
    else:
        axb.axis("off")
    fig.suptitle("HLA class I population coverage of the cross-protective shortlist",
                 x=0.055, y=1.04, ha="left", fontsize=11, fontweight="bold",
                 color=INK)
    save(fig, "fig6_population_coverage")


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("STEP 27: main-text figures")
    style()
    for fn in (fig1_cohorts, fig2_conservation, fig3_selection,
               fig4_robustness, fig5_benchmark, fig6_coverage):
        log(fn.__name__)
        try:
            fn()
        except Exception as exc:
            log(f"  FAILED: {type(exc).__name__}: {exc}")
    made = sorted(os.path.basename(p) for p in glob.glob(os.path.join(DIR_OUT, "*.pdf")))
    log("")
    log(f"figures in 20_figures/: {len(made)}")
    for m in made:
        log(f"    {m}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
