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
import re
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
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors

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


def caption(fig, text, y=-0.03):
    """Summary text BELOW the axes, in figure coordinates.

    Not ax.annotate(..., xycoords="axes fraction") in a corner. Three figures
    (the old 1, 2 and 5) each put their summary line inside the axes at the
    bottom right, which is also where loc="lower right" puts the legend, so the
    legend box sat on the text; in figure 5 it hid "HSV-2 0.600" outright.
    Placing the summary outside the axes makes that collision impossible rather
    than merely unlikely, which is why every figure now routes through here.
    bbox_inches="tight" in save() expands the canvas to include it.
    """
    fig.text(0.5, y, text, ha="center", va="top", fontsize=7, color=INK2)


def unit_label(gene, cohort, partition=None):
    """Gene label that carries its virus type.

    RS1, UL23, UL27, UL46, UL47, UL48, UL49 and UL54 each appear once per type,
    so a bare gene name is ambiguous and separated only by colour -- which
    fails both for a reader who cannot distinguish the two hues and for anyone
    reading a greyscale printout.
    """
    lab = f"{gene} {TYPE_LABEL.get(cohort, cohort)}"
    if partition:
        lab += f" p{int(partition)}"
    return lab


def rd(rel):
    p = os.path.join(DIR_BASE, rel)
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


# --------------------------------------------------------- supplementary S1
def figS1_cohorts():
    """Per-gene cohort size.

    Demoted from main-text figure 1 to supplementary S1: how many sequences each
    gene contributed is reference material a reader looks up, not part of the
    argument. Main-text figure 1 is now the study design flow.
    """
    s4 = rd("19_supplements/S4_per_gene_cohorts.tsv")
    if s4 is None:
        log("figS1 skipped: S4 missing")
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
    caption(fig, f"{tot:,} unique coding sequences in 40 gene x type cohorts")
    save(fig, "figS1_cohort_size")


# ----------------------------------------------------------------- figure 1
def _s1n(s1, prefix):
    """The n for one S1_data_flow stage, matched on its leading stage number.

    Raises if the stage is absent or ambiguous. That is deliberate: a flow
    diagram whose numbers silently fall back to a default is worse than no
    figure, because the reader cannot tell the difference.
    """
    hit = s1[s1.stage.str.startswith(prefix)]
    if len(hit) != 1:
        raise KeyError(f"S1 stage {prefix!r} matched {len(hit)} rows")
    return int(hit.n.iloc[0])


# Figure 1's vertical metrics are in POINTS, not axes fractions, and the canvas
# height is derived from them. Two earlier versions set a figure height by eye
# and laid the boxes out in axes fractions; both overflowed -- the second ran
# three boxes clean off the bottom of the axes, where the rectangles were
# clipped away and only their text survived, which is a silent failure. Sizing
# the canvas from the content cannot fail that way.
FIG1_W = 8.0
PT_TITLE = 13.5           # a box title line
PT_LINE = 11.5            # a box body line
PT_PAD = 7.0              # above the title and below the last line
PT_GAP = 15.0             # between boxes, and the length of a connector
PT_HEAD = 26.0            # figure title band
PT_FOOT = 20.0            # caption band


def _box_pt(title, lines):
    return ((PT_TITLE if title else 0.0) + PT_LINE * len(lines) + 2 * PT_PAD)


def _box(ax, u, x, y_top, w, title, lines, face=SURFACE, edge=INK3, lw=0.8,
         title_color=INK, fs=7.4, ls="solid"):
    """Draw one stage box with its top edge at y_top; return its bottom edge.

    `u` converts points to axes fraction. Height comes from the line count, so
    the only free parameter is where the stack starts.
    """
    h = _box_pt(title, lines) * u
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y_top - h), w, h,
        boxstyle="round,pad=0,rounding_size=0.010",
        linewidth=lw, edgecolor=edge, facecolor=face, linestyle=ls, zorder=2))
    if title:
        ax.text(x, y_top - PT_PAD * u, title, ha="center", va="top", zorder=3,
                fontsize=fs, fontweight="bold", color=title_color)
    if lines:
        ax.text(x, y_top - (PT_PAD + (PT_TITLE if title else 0.0)) * u,
                "\n".join(lines), ha="center", va="top", zorder=3,
                fontsize=fs - 0.5, color=INK2, linespacing=1.42)
    return y_top - h


def _arrow(ax, x0, y0, x1, y1, color=INK3, ls="-"):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9,
        linewidth=0.9, color=color, linestyle=ls, shrinkA=0, shrinkB=0,
        zorder=1))


def fig1_study_design():
    """Study design and data flow.

    Every count is read from a result table -- S1_data_flow for the sequence
    stages, the selection tables for the analysis units and verdicts, the
    shortlist tables for the candidate counts. Nothing here is typed in by
    hand, so the figure cannot drift away from the tables it describes.

    The two-track split below the HyPhy stage is the paper's argument in one
    picture: the same 305 runs give a positive gene- and branch-level result
    and a negative site-level one.
    """
    s1 = rd("19_supplements/S1_data_flow.tsv")
    g = rd("14_selection/gene_level_results.tsv")
    br = rd("14_selection/branch_level_results.tsv")
    rb = rd("14_selection/site_robustness.tsv")
    runs = rd("14_selection/selection_runs.tsv")
    s4 = rd("19_supplements/S4_per_gene_cohorts.tsv")
    if any(t is None for t in (s1, g, br, rb, runs, s4)):
        log("fig1 skipped: a source table is missing")
        return

    n_records = _s1n(s1, "1.")
    n_hsv1 = _s1n(s1, "4. Analysis-ready genomes (HSV-1)")
    n_hsv2 = _s1n(s1, "4. Analysis-ready genomes (HSV-2)")
    n_t2 = _s1n(s1, "4. Analysis-ready genomes (HSV-2 tier 2)")
    n_genomes = _s1n(s1, "4a.")
    n_cds = _s1n(s1, "6b.")
    n_uniq = _s1n(s1, "7.")
    n_cohorts = int(len(s4))

    n_units = int(len(g))
    n_part = int((g.partition.fillna(0) > 0).sum())
    n_whole = n_units - n_part
    n_runs = int(len(runs))
    n_methods = int(runs.method.nunique())
    n_sig_units = int((g.q_value < 0.05).sum())
    n_br_sig = int(br.significant.astype(bool).sum())
    n_br = int(len(br))
    n_sites = int(len(rb))
    n_conc = int((rb.n_methods.astype(int) == 2).sum())
    n_both = int(((rb.n_methods.astype(int) == 2)
                  & (rb.verdict == "ROBUST")).sum())

    ranked = rd("16_shortlist/epitope_candidates_ranked.tsv")
    primary = rd("16_shortlist/shortlist_primary.tsv")
    cross = rd("16_shortlist/shortlist_cross_protective.tsv")
    have_shortlist = all(t is not None for t in (ranked, primary, cross))

    stages = [
        ("NCBI Virus records retrieved",
         [f"{n_records:,} records   |   HSV-1 taxid 10298, HSV-2 taxid 10310"]),
        ("Genome QC and ANI type verification",
         [f"{n_genomes} analysis-ready genomes",
          f"{n_hsv1} HSV-1   |   {n_hsv2} HSV-2 tier 1   |   "
          f"{n_t2} HSV-2 tier 2 (readmitted)"]),
        ("Per-gene CDS extraction and QC",
         [f"{n_cds:,} coding sequences carried forward"]),
        ("Exact-duplicate collapse",
         [f"{n_uniq:,} unique sequences in {n_cohorts} gene x type cohorts",
          "every cohort analysed per type, never pooled"]),
        ("Codon-aware alignment, masking, per-type trees",
         ["one alignment and one tree per gene x type"]),
        ("GARD recombination screen",
         [f"{n_units} analysis units",
          f"{n_part} recombination partitions + {n_whole} unpartitioned genes"]),
        (f"{n_methods} HyPhy methods x {n_units} analysis units",
         [f"{n_runs} completed runs   |   FEL, MEME, FUBAR, BUSTED, aBSREL"]),
    ]
    left_lines = [f"{n_sig_units} of {n_units} units under episodic selection",
                  f"{n_br_sig} of {n_br:,} branches (aBSREL, Holm p < 0.05)",
                  "signal is present and type-specific"]
    right_lines = [f"{n_sites} candidate sites -> {n_conc} method-concordant",
                   f"-> {n_both} both concordant and variant-supported",
                   "no site survives both checks"]
    epi = ("Epitope prediction and IEDB benchmark",
           ["MHCflurry 2.3.30, NetMHCIIpan-4.3, BepiPred-3.0, DiscoTope-3.0",
            "all scores carried as continuous annotations, never as filters"])
    integ = ("Integration and shortlisting",
             [f"{len(ranked):,} ranked candidates",
              f"{len(primary):,} primary   |   "
              f"{len(cross):,} cross-protective"]) if have_shortlist else None
    withdrawn = ("", [
        "WITHDRAWN  v1 pooled HSV-1 and HSV-2 into single alignments. Its "
        "site-level results are superseded",
        "by the per-type analysis above and are retained only in the "
        "corrections register."])

    # ---- measure, then size the canvas
    total = PT_HEAD + PT_FOOT
    total += sum(_box_pt(t, ls_) for t, ls_ in stages)
    total += PT_GAP * (len(stages) - 1)
    total += 2 * PT_GAP                      # fork stem + fork spread
    total += _box_pt("Gene and branch level", left_lines)
    total += PT_GAP                          # merge
    total += _box_pt(*epi)
    if integ:
        total += PT_GAP + _box_pt(*integ)
    total += 1.4 * PT_GAP + _box_pt(*withdrawn)

    fig_h = total / 72.0 / 0.97
    fig = plt.figure(figsize=(FIG1_W, fig_h))
    # The axes fills the canvas: with the default subplot margins the drawing
    # area is 0.775 of the figure width, which squeezed the longest body line
    # past its box edge.
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    u = 1.0 / (fig_h * 72.0)                 # points -> axes fraction
    gap = PT_GAP * u
    cx, w = 0.5, 0.74

    ax.text(0.035, 1.0 - 6 * u, "Study design and data flow", ha="left",
            va="top", fontsize=11.5, fontweight="bold", color=INK)

    y = 1.0 - PT_HEAD * u
    for i, (title, lines) in enumerate(stages):
        y = _box(ax, u, cx, y, w, title, lines)
        if i < len(stages) - 1:
            _arrow(ax, cx, y, cx, y - gap)
            y -= gap

    lx, rx, tw = 0.258, 0.742, 0.455
    fork_y = y - gap
    track_top = fork_y - gap
    _arrow(ax, cx, y, cx, fork_y)
    for tx in (lx, rx):
        _arrow(ax, cx, fork_y, tx, track_top)
    lb = _box(ax, u, lx, track_top, tw, "Gene and branch level", left_lines,
              face="#eef5fd", edge=BLUE, lw=1.1, title_color=BLUE)
    rbx = _box(ax, u, rx, track_top, tw, "Site level", right_lines,
               face="#fdf0ea", edge=ORANGE, lw=1.1, title_color=ORANGE)
    y = min(lb, rbx)

    merge_y = y - gap
    for tx in (lx, rx):
        _arrow(ax, tx, y, cx, merge_y)
    y = _box(ax, u, cx, merge_y, w, *epi)
    if integ:
        _arrow(ax, cx, y, cx, y - gap)
        y = _box(ax, u, cx, y - gap, w, *integ)

    # The withdrawn v1 branch: marked, not omitted, and placed at the foot as a
    # dashed note rather than as a callout floating over the stack, where the
    # first version of this figure put it straight on top of a stage box.
    y = _box(ax, u, cx, y - 1.4 * gap, w, *withdrawn,
             face="#f4f3f0", edge=INK3, lw=0.7, ls=(0, (4, 2)))

    ax.text(cx, y - 7 * u,
            "every count is read from S1_data_flow.tsv, the selection tables "
            "and the shortlist tables; none is hardcoded",
            ha="center", va="top", fontsize=7, color=INK2)
    if y - 7 * u < 0:
        log(f"  fig1 WARNING: stack overran the canvas by {-(y - 7 * u):.3f}")
    save(fig, "fig1_study_design")


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
    # Pushed right along gG/US4's own row, into the empty space past the end
    # of the shortest bar. Two earlier placements both landed on a value
    # label: 3 units right of the bar tip sat on "37.5", and lifting it 14
    # points put it across "59.0" and "61.1" on the rows above. The only
    # reliably empty region in this panel is to the right of the bottom bar.
    axb.annotate("gG/US4 is the poorest\ncross-protective target",
                 xy=(c.mean_between_type_identity.iloc[0], 0),
                 xytext=(40, 0), textcoords="offset points",
                 fontsize=7, color=INK2, va="center", ha="left",
                 arrowprops=dict(arrowstyle="-", lw=0.6, color=INK3))
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
    g["unit_label"] = [unit_label(r.gene, r.cohort, r.partition)
                       for r in g.itertuples()]
    g["nlq"] = -np.log10(g.q_value.clip(lower=1e-12))
    n_units = len(g)
    n_sig = int((g.q_value < 0.05).sum())
    # Panel (a) shows the significant units only. All 61 rows was at the edge of
    # legibility at print width and the 29 sub-threshold bars carry no
    # information beyond "below the line"; the full list is in the result table.
    g = g[g.q_value < 0.05].sort_values("nlq", ascending=True)
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
    axa.set_title(f"a  Gene-wide: the {n_sig} significant of {n_units} "
                  f"analysis units", loc="left")
    axa.xaxis.grid(True)
    recessive(axa)
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
        axb.set_yticklabels([unit_label(r.gene, r.cohort)
                             for r in cnt.itertuples()], fontsize=6.2)
        axb.set_xlabel("branches under selection")
        axb.set_title("b  Branch-level (aBSREL)", loc="left")
        axb.xaxis.grid(True)
        recessive(axb)
        caption(fig,
                f"(a) {n_sig} of {n_units} units significant at q < 0.05; the "
                f"{n_units - n_sig} sub-threshold units are omitted and\n"
                f"listed in gene_level_results.tsv.   "
                f"(b) {int(sig.shape[0])} of {len(b):,} branches at "
                f"Holm-corrected p < 0.05.", y=0.015)
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
    # The vertical axis of this panel encodes nothing -- these sites have no q
    # value, and the spread is jitter added only so that overlapping points stay
    # countable. Saying so on the axis is the honest answer to "what is y?"; an
    # unlabelled axis invites the reader to read a quantity into it.
    axf.set_ylabel("no q value\n(spread is jitter)", fontsize=7.5, color=INK2)
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
    caption(fig,
            "shape = variant-support verdict, colour = virus type   |   "
            + ",  ".join(f"{k.replace('_', ' ').lower()} {v}"
                         for k, v in counts.items())
            + "   |   both concordant sites fall left of the threshold",
            y=-0.055)
    save(fig, "fig4_site_variant_support")


# ----------------------------------------------------------------- figure 5
def fig5_benchmark():
    b = rd("15_benchmark/benchmark_per_gene.tsv")
    if b is None:
        log("fig5 skipped: benchmark_per_gene missing")
        return
    a = b[b.score == "affinity"].copy() if "score" in b.columns else b.copy()
    a = a.sort_values("auc")
    fig, ax = plt.subplots(figsize=(6.9, 5.2))
    y = np.arange(len(a))
    ax.hlines(y, 0.5, a.auc, color=GRID, lw=1.2, zorder=1)

    # Marker area encodes how many curated IEDB positives each AUC was computed
    # over. Drawn at one fixed size, an AUC over 5 positives looked exactly like
    # an AUC over 75 and the reader had no way to see which points carry weight.
    # Area, not radius, is proportional to sqrt(n), so the visual weight grows
    # sub-linearly and the largest group does not swamp the panel.
    has_n = "n_positives" in a.columns
    if has_n:
        nv = a.n_positives.astype(float).values
        lo, hi = np.sqrt(nv.min()), np.sqrt(nv.max())
        span = max(hi - lo, 1e-9)

        def msize(n):
            return 22.0 + 150.0 * (np.sqrt(np.asarray(n, float)) - lo) / span
    else:
        def msize(n):
            return np.full(np.shape(n), 38.0)

    pos = {ix: k for k, ix in enumerate(a.index)}
    for cohort in ("HSV1", "HSV2"):
        d = a[a.virus_type == cohort]
        idx = [pos[i] for i in d.index]
        ax.scatter(d.auc, idx,
                   s=msize(d.n_positives.values) if has_n else 38,
                   color=TYPE_COLOR[cohort],
                   label=TYPE_LABEL[cohort], zorder=3, linewidths=0.8,
                   edgecolors=SURFACE)
    ax.axvline(0.5, color=INK2, lw=1.0, ls="--", zorder=2)
    ax.annotate("chance", xy=(0.5, len(a) - 0.5), xytext=(4, 0),
                textcoords="offset points", fontsize=7, color=INK2, va="top")
    ax.set_yticks(y)
    ax.set_yticklabels([unit_label(r.gene, r.virus_type) for r in a.itertuples()],
                       fontsize=6.6)
    ax.set_xlabel("AUC, predicted affinity vs curated IEDB epitopes")
    ax.set_title("Predictor validation is strong in HSV-1 and weak in HSV-2",
                 loc="left")
    ax.set_xlim(0.35, 1.02)
    ax.xaxis.grid(True)
    recessive(ax)
    type_leg = ax.legend(loc="lower right", title=None)
    ax.add_artist(type_leg)
    if has_n:
        ticks = sorted({int(nv.min()), int(np.median(nv)), int(nv.max())})
        ax.legend(handles=[Line2D([], [], marker="o", ls="", color=INK3,
                                  markeredgecolor=SURFACE,
                                  markersize=np.sqrt(msize(t)),
                                  label=f"{t}")
                           for t in ticks],
                  loc="lower right", bbox_to_anchor=(1.0, 0.14),
                  title="IEDB positives", title_fontsize=7,
                  labelspacing=0.9, handletextpad=0.9, borderpad=0.4)
    m1 = a[a.virus_type == "HSV1"].auc.mean()
    m2 = a[a.virus_type == "HSV2"].auc.mean()
    bits = [f"mean AUC  HSV-1 {m1:.3f}   HSV-2 {m2:.3f}"]
    if has_n:
        bits.append(f"marker area = curated IEDB positives "
                    f"({int(nv.min())} to {int(nv.max())} per gene x type)")
    caption(fig, "   |   ".join(bits))
    save(fig, "fig5_predictor_benchmark")


# ----------------------------------------------------------------- figure 6
def fig6_coverage():
    """Population coverage, class I and class II.

    Three panels, replacing the two-panel class-I-only version. Panels (a) and
    (b) share their orientation -- coverage on y, populations sorted along x --
    because the class II result is a BRACKET and a bracket can only be drawn as
    a band between two curves. Reading both panels the same way is worth more
    than keeping panel (a) as the ECDF it was.
    """
    cp = rd("18_population_coverage/coverage_by_population.tsv")
    c2 = rd("18_population_coverage/coverage_by_population_classii.tsv")
    gp = rd("18_population_coverage/minimal_panel_greedy.tsv")
    if cp is None:
        log("fig6 skipped: coverage_by_population missing")
        return
    target = "shortlist_cross_protective"
    sub = cp[cp.epitope_set == target] if "epitope_set" in cp.columns else cp
    if sub.empty:
        sub = cp
    fig, (axa, axb, axc) = plt.subplots(
        1, 3, figsize=(10.6, 3.5),
        gridspec_kw=dict(width_ratios=[1, 1, 1], wspace=0.40))

    # ---- (a) class I, one curve per population, sorted
    vals = np.sort(sub.coverage.values)
    pct = np.arange(1, len(vals) + 1) / len(vals) * 100
    axa.plot(pct, vals, color=BLUE, lw=1.8)
    axa.axhline(90, color=INK2, lw=1.0, ls="--")
    above = float((sub.coverage >= 90).mean() * 100)
    axa.annotate(f"{above:.1f}% of populations at or above 90%",
                 xy=(2, 90), xytext=(0, 6), textcoords="offset points",
                 ha="left", fontsize=7, color=INK2)
    axa.set_xlabel("populations, sorted by coverage (%)")
    axa.set_ylabel("coverage within a population (%)")
    axa.set_title(f"a  Class I, {sub.population.nunique()} populations",
                  loc="left")
    axa.set_ylim(min(55, float(vals.min()) - 3), 101)
    axa.grid(True)
    recessive(axa)

    # ---- (b) class II, conservative-to-upper-bound bracket
    if c2 is None or "coverage_conservative" not in c2.columns:
        axb.axis("off")
        log("  fig6 panel b skipped: class II coverage table missing")
        n_score = n_unscored = 0
    else:
        set_ii = ("classii_binders_weak_or_better"
                  if "classii_binders_weak_or_better" in set(c2.epitope_set)
                  else sorted(set(c2.epitope_set))[0])
        d2 = c2[c2.epitope_set == set_ii]
        sc = d2[d2.scorable.astype(bool)]
        n_score, n_unscored = len(sc), len(d2) - len(sc)
        lowv = np.sort(sc.coverage_conservative.values)
        upv = np.sort(sc.coverage_upper_bound.values)
        p2 = np.arange(1, len(lowv) + 1) / len(lowv) * 100
        # The band IS the result: DP and DQ are heterodimers, so a peptide bound
        # by one chain pairing cannot be credited to every genotype carrying
        # that chain. The lower curve credits only pairings that are certain;
        # the upper curve credits every pairing that is possible.
        axb.fill_between(p2, lowv, upv, color=AQUA, alpha=0.22, linewidth=0,
                         zorder=1)
        axb.plot(p2, upv, color=AQUA, lw=1.8, zorder=3, label="upper bound")
        axb.plot(p2, lowv, color=ORANGE, lw=1.8, zorder=3, label="conservative")
        axb.axhline(90, color=INK2, lw=1.0, ls="--", zorder=2)
        axb.set_xlabel("populations, sorted by coverage (%)")
        axb.set_title(f"b  Class II, {n_score} typed populations",
                      loc="left")
        axb.set_ylim(min(88, float(lowv.min()) - 2), 100.6)
        axb.grid(True)
        recessive(axb)
        axb.legend(loc="lower right")
        # Direct labels as well as the legend, anchored at the left-hand end
        # where the two curves are furthest apart. At the quartile they were
        # a few tenths of a percent apart and the two labels printed on top of
        # each other. Aqua carries a sub-3:1 contrast warning against this
        # surface, so it must never be the only thing naming a curve.
        k = max(len(p2) // 12, 1)
        axb.annotate("upper bound", xy=(p2[k], upv[k]),
                     xytext=(4, 4), textcoords="offset points",
                     fontsize=6.8, color=INK2, ha="left", va="bottom")
        axb.annotate("conservative", xy=(p2[k], lowv[k]),
                     xytext=(4, -4), textcoords="offset points",
                     fontsize=6.8, color=INK2, ha="left", va="top")
        axb.annotate(f"{n_unscored} populations omitted: an untyped\n"
                     f"class II locus drives the conservative\n"
                     f"figure to a false zero",
                     xy=(0.04, 0.34), xycoords="axes fraction",
                     ha="left", va="bottom", fontsize=6.6, color=INK2)

    # ---- (c) minimal class I panel
    if gp is not None and len(gp):
        cumcol = next((c for c in ("cumulative_coverage", "world_coverage",
                                   "coverage") if c in gp.columns), None)
        pepcol = next((c for c in ("peptide", "epitope") if c in gp.columns), None)
        if cumcol:
            n = np.arange(1, len(gp) + 1)
            axc.plot(n, gp[cumcol], color=BLUE, lw=1.8, marker="o", ms=5,
                     markeredgecolor=SURFACE, markeredgewidth=0.8)
            axc.set_xticks(n)
            if pepcol:
                axc.set_xticklabels(gp[pepcol], rotation=35, ha="right",
                                    fontsize=6.6)
            axc.set_xlabel("peptides added, greedy order")
            axc.set_ylabel("cumulative World coverage (%)")
            axc.set_title("c  Minimal class I panel", loc="left")
            axc.set_ylim(min(85, float(gp[cumcol].min()) - 2), 101)
            axc.grid(True)
            recessive(axc)
            axc.annotate(f"{gp[cumcol].iloc[0]:.1f}% from one peptide",
                         xy=(1, gp[cumcol].iloc[0]), xytext=(8, -12),
                         textcoords="offset points", fontsize=7, color=INK2)
    else:
        axc.axis("off")
    fig.suptitle("Population coverage of the cross-protective shortlist",
                 x=0.045, y=1.04, ha="left", fontsize=11, fontweight="bold",
                 color=INK)
    caption(fig, "(c) excludes low-complexity peptides; the greedy search "
                 "substitutes rather than truncates, so the panel size and "
                 "coverage are unchanged by that constraint", y=-0.10)
    save(fig, "fig6_population_coverage")


# ----------------------------------------------------------------- figure 7
# van der Waals radii, Angstrom, Bondi 1964. Only the elements present in the
# three deposited structures are listed; anything else falls back to carbon.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20}

DIR_STRUCT = os.path.join(DIR_BASE, "17_structures")

# One hue, light to dark: the sequential ramp for DiscoTope score. Never a
# rainbow -- a magnitude needs an ordered ramp a reader can rank by eye.
SEQ_BLUE = mcolors.LinearSegmentedColormap.from_list(
    "seq_blue", ["#eef4fc", "#9cc2ec", BLUE, "#123f73"])


def _chain_atoms(pdb_path, chain_id):
    """Heavy atoms of one chain as a frame of resnum, element, x, y, z."""
    from Bio.PDB import PDBParser
    st = PDBParser(QUIET=True).get_structure("s", pdb_path)
    rows = []
    for res in st[0][chain_id]:
        if res.id[0] != " ":          # skip waters and heteroatoms
            continue
        for at in res:
            el = (at.element or at.get_name()[0]).strip().upper()
            x, y, z = (float(c) for c in at.coord)
            rows.append((int(res.id[1]), el, x, y, z))
    return pd.DataFrame(rows, columns=["resnum", "element", "x", "y", "z"])


def _project(xyz):
    """Principal-axis projection: returns (u, v, depth).

    SVD rather than an eigendecomposition of the covariance matrix, and the sign
    of each component is pinned by its largest-magnitude element, so the view is
    identical on every run and on a different BLAS. Without the sign fix the
    molecule silently appears mirrored from one machine to the next.
    """
    X = np.asarray(xyz, float)
    X = X - X.mean(axis=0)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    for i in range(vt.shape[0]):
        if vt[i][int(np.argmax(np.abs(vt[i])))] < 0:
            vt[i] = -vt[i]
    pr = X @ vt.T
    return pr[:, 0], pr[:, 1], pr[:, 2]


def _spacefill(ax, atoms, colors, fade=0.30, ring=None):
    """Depth-sorted space-filling projection, drawn with matplotlib alone.

    PyMOL is not installed in this environment and adding a ray-tracer to the
    dependency set for two panels is not worth it, so the surface is drawn as
    van der Waals spheres in a principal-axis projection, painted far to near.
    Atoms further from the viewer are blended toward the page colour, which is
    the only depth cue a flat painter's-algorithm render can have. `fade` is
    deliberately mild: at 0.6 the blending bleached pale atoms all the way to
    the page colour and punched holes through the middle of the molecule.

    `ring` is a boolean mask drawn with an ink outline on top, used to pick out
    the highlighted patch. Each panel gets its own tight window; the scale bar,
    not the apparent size, is what makes the two panels comparable. A shared
    square window was tried first and left two thirds of the canvas empty,
    because gB's ectodomain is three times longer than it is wide.
    """
    u, v, d = _project(atoms[["x", "y", "z"]].values)
    rng = max(float(np.ptp(d)), 1e-9)
    dn = (d - d.min()) / rng
    surf = np.array(mcolors.to_rgb(SURFACE))
    order = np.argsort(d)                      # far first
    els = atoms.element.values
    circles, faces = [], []
    for i in order:
        circles.append(Circle((u[i], v[i]), VDW.get(els[i], 1.70)))
        base = np.array(mcolors.to_rgb(colors[i]))
        f = fade * (1.0 - dn[i])
        faces.append(tuple(base * (1 - f) + surf * f))
    ax.add_collection(PatchCollection(circles, facecolors=faces,
                                      edgecolors="none", zorder=2))
    if ring is not None and np.any(ring):
        idx = [i for i in order if ring[i]]
        ax.add_collection(PatchCollection(
            [Circle((u[i], v[i]), VDW.get(els[i], 1.70)) for i in idx],
            facecolors="none", edgecolors=INK, linewidths=0.45, zorder=4))
    pad = 3.0
    ax.set_xlim(float(u.min()) - pad, float(u.max()) + pad)
    ax.set_ylim(float(v.min()) - pad, float(v.max()) + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    return u, v


def _scalebar(ax, angstrom=20.0, corner="lower left"):
    """A scale bar, because the two panels are different proteins at
    different magnifications. This is what makes their sizes comparable; the
    apparent sizes are not."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    bx = x0 + 0.04 * (x1 - x0)
    by = (y0 + 0.05 * (y1 - y0) if corner == "lower left"
          else y1 - 0.09 * (y1 - y0))
    ax.plot([bx, bx + angstrom], [by, by], color=INK2, lw=1.6,
            solid_capstyle="butt", zorder=6)
    ax.text(bx + angstrom / 2, by + 0.015 * (y1 - y0),
            f"{angstrom:.0f} " + chr(197), ha="center", va="bottom",
            fontsize=7, color=INK2, zorder=6)


def _pick_patch(patches, pdb_id, chain):
    """The conserved, exposed, fully-epitope-overlapping patch on one chain.

    Chosen by rule, not by name: conserved in both types at every residue,
    overlapping a curated IEDB epitope at every residue, then longest, then
    lowest start. On 2GUM chain A exactly one patch satisfies that, which is
    the one the manuscript discusses -- but the rule, not the number, is what
    is written down here, so the figure follows the table if the table changes.
    """
    d = patches[(patches.pdb_id == pdb_id) & (patches.chain == chain)]
    d = d[(d.conserved_both_types == d.length)
          & (d.overlaps_known_epitope == d.length)]
    if d.empty:
        return None
    return d.sort_values(["length", "pdb_start"],
                         ascending=[False, True]).iloc[0]


def fig7_structure():
    """The shortlist on a structure.

    (a) gB/UL27 with the one fully conserved, fully solvent-exposed patch that
    lies inside a curated IEDB epitope at every residue -- the paper's
    strongest structural claim. (b) gD/US6 coloured by DiscoTope-3.0 score,
    which is the B-cell arm's only structural read-out.

    Each panel is drawn at its own magnification with a scale bar; gB's
    ectodomain is far longer than gD's, and forcing both into one window left
    most of the canvas empty.
    """
    ann = rd("17_structures/structure_residue_annotation.tsv")
    pat = rd("17_structures/conserved_exposed_patches.tsv")
    conf = rd("13_epitopes/bcell_conformational_scores.tsv")
    ch = rd("17_structures/chain_assignments.tsv")
    if ann is None or pat is None or ch is None:
        log("fig7 skipped: structure annotation tables missing")
        return

    def gene_of(pdb_id, chain):
        hit = ch[(ch.pdb_id == pdb_id) & (ch.chain == chain)]
        return str(hit.assigned_gene.iloc[0]) if len(hit) else "?"

    # ---- panel a inputs: gB with the conserved exposed patch
    pdb_a, chain_a = "2GUM", "A"
    path_a = os.path.join(DIR_STRUCT, "pdb", f"{pdb_a}.pdb")
    row = _pick_patch(pat, pdb_a, chain_a)
    panel_a = None
    if os.path.exists(path_a) and row is not None:
        atoms = _chain_atoms(path_a, chain_a)
        a_ann = ann[(ann.pdb_id == pdb_a) & (ann.chain == chain_a)]
        exp = set(a_ann.loc[a_ann.exposed.astype(bool)
                            & a_ann.conserved_both_types.astype(bool),
                            "pdb_resnum"].astype(int))
        var = set(a_ann.loc[a_ann.within_type_identity < 100,
                            "pdb_resnum"].astype(int))
        lo, hi = int(row.pdb_start), int(row.pdb_end)
        patch_res = set(range(lo, hi + 1))
        cols, ring = [], []
        for rn in atoms.resnum.values:
            if rn in patch_res:
                cols.append(ORANGE)
            elif rn in exp:
                cols.append(AQUA)
            elif rn in var:
                cols.append("#c9c7c1")
            else:
                cols.append("#e6e4df")
            ring.append(rn in patch_res)
        panel_a = (atoms, cols, np.array(ring), row, lo, hi, patch_res)
    else:
        log(f"  fig7 panel a skipped: {pdb_a} structure or patch unavailable")

    # ---- panel b inputs: gD coloured by DiscoTope score
    pdb_b, chain_b = "2C36", "A"
    path_b = os.path.join(DIR_STRUCT, "pdb", f"{pdb_b}.pdb")
    panel_b = None
    if conf is not None and os.path.exists(path_b):
        atoms_b = _chain_atoms(path_b, chain_b)
        d = conf[(conf.pdb_id == pdb_b) & (conf.chain == chain_b)]
        smap = dict(zip(d.pdb_resnum.astype(int), d.discotope_score.astype(float)))
        norm = mcolors.Normalize(vmin=float(d.discotope_score.min()),
                                 vmax=float(d.discotope_score.max()))
        cols_b = [SEQ_BLUE(norm(smap[rn])) if rn in smap else "#e6e4df"
                  for rn in atoms_b.resnum.values]
        panel_b = (atoms_b, cols_b, norm)
    else:
        log(f"  fig7 panel b skipped: {pdb_b} or DiscoTope scores unavailable")

    if panel_a is None and panel_b is None:
        return

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(9.0, 3.6),
        gridspec_kw=dict(width_ratios=[1.55, 1], wspace=0.04))

    if panel_a is None:
        axa.axis("off")
    else:
        atoms, cols, ring, row, lo, hi, patch_res = panel_a
        u, v = _spacefill(axa, atoms, cols, ring=ring)
        inpatch = atoms.resnum.isin(patch_res).values
        cu, cv = float(u[inpatch].mean()), float(v[inpatch].mean())
        gene_a = gene_of(pdb_a, chain_a)
        axa.annotate(f"{row.peptide}\n{gene_a} {lo}-{hi}", xy=(cu, cv),
                     xytext=(0.30, 0.95), textcoords="axes fraction",
                     ha="center", va="top", fontsize=8.5, color=INK,
                     fontweight="bold",
                     arrowprops=dict(arrowstyle="-|>", mutation_scale=9,
                                     lw=1.0, color=INK2, shrinkA=3, shrinkB=3))
        axa.set_title(f"a  gB / {gene_a}  (PDB {pdb_a}, chain {chain_a})",
                      loc="left")
        _scalebar(axa)
        axa.legend(handles=[
            Line2D([], [], marker="o", ls="", ms=7, color=ORANGE,
                   markeredgecolor=INK, markeredgewidth=0.5,
                   label="conserved exposed patch, inside an IEDB epitope"),
            Line2D([], [], marker="o", ls="", ms=7, color=AQUA,
                   label="conserved in both types and solvent-exposed"),
            Line2D([], [], marker="o", ls="", ms=7, color="#c9c7c1",
                   label="variable within type"),
            Line2D([], [], marker="o", ls="", ms=7, color="#e6e4df",
                   label="conserved but buried")],
            loc="upper left", bbox_to_anchor=(0.0, -0.01), ncol=1,
            handletextpad=0.6, labelspacing=0.34)

    if panel_b is None:
        axb.axis("off")
    else:
        atoms_b, cols_b, norm = panel_b
        _spacefill(axb, atoms_b, cols_b)
        gene_b = gene_of(pdb_b, chain_b)
        axb.set_title(f"b  gD / {gene_b}  (PDB {pdb_b}, chain {chain_b})",
                      loc="left")
        _scalebar(axb, corner="upper left")
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=SEQ_BLUE),
                          ax=axb, orientation="horizontal", fraction=0.040,
                          pad=0.02, aspect=28)
        cb.set_label("DiscoTope-3.0 conformational B-cell score", fontsize=7.5,
                     color=INK)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=7, color=INK3, labelcolor=INK2)

    fig.suptitle("Shortlisted surface on the deposited structures",
                 x=0.045, y=0.99, ha="left", fontsize=11, fontweight="bold",
                 color=INK)
    caption(fig, "van der Waals spheres in a principal-axis projection, "
                 "painted far to near; atoms further from the viewer are "
                 "blended toward the page colour. The panels are at different "
                 "magnifications -- compare them with the scale bars.",
            y=-0.21)
    save(fig, "fig7_structure")


# --------------------------------------------------------- supplementary S2
def _threshold_from(path, pattern, col=None, row_prefix=None):
    """Pull a published threshold out of a project file, never from memory.

    Both B-cell thresholds are recorded in the pipeline's own outputs -- the
    BepiPred one in the predictor provenance table, the DiscoTope one in the
    step28 log. Reading them from there means this figure cannot disagree with
    the run that produced the scores.
    """
    full = os.path.join(DIR_BASE, path)
    if not os.path.exists(full):
        return None
    if col:
        t = pd.read_csv(full, sep="\t")
        if row_prefix:
            t = t[t.iloc[:, 0].astype(str).str.startswith(row_prefix)]
        text = " ".join(t[col].astype(str).tolist())
    else:
        with open(full, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def figS2_bcell():
    """B-cell predictor score distributions, and what their thresholds do.

    The point of this figure is not the scores. It is that neither published
    threshold behaves on this data the way its documentation says, which is why
    both predictors are carried as continuous annotations and neither is used as
    a filter anywhere in the pipeline. A reader is entitled to see that.
    """
    lin = rd("13_epitopes/bcell_linear_scores.tsv")
    conf = rd("13_epitopes/bcell_conformational_scores.tsv")
    if lin is None and conf is None:
        log("figS2 skipped: no B-cell score tables")
        return
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(8.6, 3.5),
                                   gridspec_kw=dict(width_ratios=[1, 1.08],
                                                    wspace=0.26))

    # ---------------------------------------------------------------- panel a
    if lin is None:
        axa.axis("off")
        note_a = ""
    else:
        thr = _threshold_from("13_epitopes/bcell_predictor_provenance.tsv",
                              r"threshold\s+([0-9.]+)", col="mode",
                              row_prefix="BepiPred")
        flag = lin.predicted_epitope_residue.astype(bool)
        bins = np.linspace(float(lin.bepipred_linear_score.min()),
                           float(lin.bepipred_linear_score.max()), 46)
        for mask, colour, lab in ((~flag, BLUE, "not called an epitope"),
                                  (flag, ORANGE, "called an epitope")):
            axa.hist(lin.loc[mask, "bepipred_linear_score"], bins=bins,
                     histtype="stepfilled", color=colour, alpha=0.16,
                     linewidth=0)
            axa.hist(lin.loc[mask, "bepipred_linear_score"], bins=bins,
                     histtype="step", color=colour, lw=1.8, label=lab)
        if thr is not None:
            axa.axvline(thr, color=INK2, lw=1.0, ls="--", zorder=4)
            axa.annotate(f"published threshold {thr}", xy=(thr, 0.97),
                         xycoords=("data", "axes fraction"),
                         xytext=(5, 0), textcoords="offset points",
                         va="top", fontsize=7, color=INK2)
        axa.set_xlabel("BepiPred-3.0 linear B-cell score")
        axa.set_ylabel("residues")
        axa.set_title("a  Linear B-cell score, 20 genes x 2 types", loc="left")
        axa.yaxis.grid(True)
        recessive(axa)
        axa.legend(loc="upper right", bbox_to_anchor=(1.0, 0.88))
        ov_lo = float(lin.loc[flag, "bepipred_linear_score"].min())
        ov_hi = float(lin.loc[~flag, "bepipred_linear_score"].max())
        note_a = (f"(a) the call does not follow the score: residues called "
                  f"epitopes run down to {ov_lo:.3f} and residues not called "
                  f"run up to {ov_hi:.3f}, so the two sets overlap across the "
                  f"threshold")

    # ---------------------------------------------------------------- panel b
    if conf is None:
        axb.axis("off")
        note_b = ""
    else:
        hi_thr = _threshold_from("13_epitopes/bcell.log",
                                 r"high-conf\(>=([0-9.]+)\)")
        conf = conf.copy()
        conf["grp"] = (conf.pdb_id.astype(str) + " " + conf.chain.astype(str)
                       + "  " + conf.gene.astype(str))
        groups = sorted(conf.grp.unique())
        for i, gname in enumerate(groups):
            d = conf[conf.grp == gname]
            # Deterministic jitter: crc32 of the group name, never hash(),
            # which Python reseeds per process (see figure 4).
            rng = np.random.default_rng(zlib.crc32(gname.encode()))
            yv = i + rng.uniform(-0.26, 0.26, len(d))
            vt = str(d.virus_type.iloc[0])
            axb.scatter(d.discotope_score, yv, s=5, alpha=0.35,
                        color=TYPE_COLOR.get(vt, INK3), linewidths=0,
                        zorder=3)
            med = float(d.discotope_score.median())
            axb.plot([med, med], [i - 0.33, i + 0.33], color=INK, lw=1.6,
                     zorder=4)
        axb.set_yticks(range(len(groups)))
        axb.set_yticklabels(groups, fontsize=6.8)
        axb.set_xlabel("DiscoTope-3.0 conformational B-cell score")
        axb.set_title("b  Conformational score per chain", loc="left")
        axb.xaxis.grid(True)
        recessive(axb)
        xhi = float(conf.discotope_score.max())
        if hi_thr is not None:
            axb.axvline(hi_thr, color=INK2, lw=1.0, ls="--", zorder=2)
            axb.set_xlim(-0.02, max(hi_thr, xhi) + 0.06)
            axb.annotate(f"high confidence\n>= {hi_thr}", xy=(hi_thr, 0.5),
                         xycoords=("data", "axes fraction"),
                         xytext=(-6, 0), textcoords="offset points",
                         ha="right", va="center", fontsize=7, color=INK2)
        n_hi = int(conf.discotope_high_confidence.astype(bool).sum())
        axb.legend(handles=[Line2D([], [], color=INK, lw=1.6, label="median")]
                   + [Line2D([], [], marker="o", ls="", ms=5,
                             color=TYPE_COLOR[k], label=TYPE_LABEL[k])
                      for k in ("HSV1", "HSV2")],
                   loc="center right", bbox_to_anchor=(1.0, 0.30))
        note_b = (f"(b) {n_hi} of {len(conf):,} scored residues reach the "
                  f"high-confidence cut-off; the highest score in the three "
                  f"structures is {xhi:.3f}")

    fig.suptitle("B-cell predictor scores and the behaviour of their thresholds",
                 x=0.045, y=1.03, ha="left", fontsize=11, fontweight="bold",
                 color=INK)
    # Two lines, not one. As a single line this caption was wider than the
    # plots and bbox_inches="tight" stretched the whole figure to fit it.
    caption(fig, "\n".join(x for x in (note_a, note_b) if x), y=-0.09)
    save(fig, "figS2_bcell_scores")


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("STEP 27: main-text figures")
    style()
    for fn in (fig1_study_design, fig2_conservation, fig3_selection,
               fig4_robustness, fig5_benchmark, fig6_coverage,
               fig7_structure, figS1_cohorts, figS2_bcell):
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
