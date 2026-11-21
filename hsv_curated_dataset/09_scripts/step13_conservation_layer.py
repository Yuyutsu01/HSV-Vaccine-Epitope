"""
Step 13 (PENDING WORK B1): Per-residue conservation layer.

Proposal pipeline step 6: per-residue Shannon entropy and percent identity across all
isolates, within HSV-1, within HSV-2, and between the two types. The repository previously
held only gene-level means (07_qc_reports/09_conservation_metric_diagnostic.tsv); the
per-residue layer that the epitope ranking depends on did not exist.

Computed on the rebuilt, masked COMBINED alignments (step09b/09c), which is the one
legitimate use of a pooled HSV-1 + HSV-2 alignment: establishing residue correspondence
between the two types. No dN/dS or recombination inference is performed on pooled data
(CORRECTIONS_AND_PENDING_WORK.md, A2).

Both reference rows are present in those alignments, so each column carries an HSV-1
(NC_001806.2) and an HSV-2 (NC_001798.2) residue number.

Those residue numbers are counted on the UNMASKED combined alignment and indexed through
the retained-column list that step09c records in masking_report_v2.tsv, NOT counted on the
masked alignment. Counting on the masked alignment is wrong and was the original defect:
masking removes columns, so a reference residue in a dropped column is never counted and
every number after it is short by one, the drift accumulating down the protein. That put
7,343 of 22,955 rows on the wrong residue, worst in the heavily masked genes (RS1, RL2,
UL27, UL30), and it propagated into every downstream join on (gene, type, ref_residue) --
the shortlist conservation annotation, the cross-protective subset, the priority score, the
minimal panel and the structural mapping. `assert_reference_numbering` now checks every
emitted (ref_residue, ref_aa) pair against the GenBank translation and aborts the run on
any mismatch; that check is what the step was missing. Entropy is reported twice:
  - haplotype-weighted: one count per unique deduplicated sequence
  - isolate-weighted:   each sequence weighted by its cluster_size, recovering the true
                        population frequency of each residue

Outputs:
  - 11_conservation/per_residue_conservation.tsv   (master table, one row per gene/column)
  - 11_conservation/gene_level_conservation.tsv    (supersedes 09_conservation_metric_diagnostic.tsv)
  - 11_conservation/conserved_windows.tsv          (9-mer windows conserved in both types)
  - 11_conservation/conservation.log
"""

import os
import re
import sys
import glob
import math
from collections import Counter
from datetime import datetime

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_COMB = os.path.join(DIR_BASE, "06_alignments_v2", "masked", "combined")
DIR_UNMASKED = os.path.join(DIR_BASE, "06_alignments_v2", "combined")
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")
MASK_REPORT = os.path.join(DIR_BASE, "06_alignments_v2", "masking_report_v2.tsv")
DIR_OUT = os.path.join(DIR_BASE, "11_conservation")

REF_ACC = {"HSV1": "NC_001806.2", "HSV2": "NC_001798.2"}

# A window is reported as a cross-protective candidate when every residue in it is
# conserved in both types at >= this frequency and identical between the types.
WINDOW = 9
CONS_THRESHOLD = 0.95

LOG_PATH = os.path.join(DIR_OUT, "conservation.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def read_fasta(path):
    ids, seqs, cur = [], [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    seqs.append(cur)
                ids.append(line[1:])
                cur = ""
            elif cur is not None:
                cur += line
    if cur is not None:
        seqs.append(cur)
    return ids, seqs


def translate_codons(seq):
    """Codon-wise translation of an aligned row; '-' for a gap codon."""
    out = []
    for i in range(0, len(seq), 3):
        cod = seq[i:i + 3]
        if cod == "---":
            out.append("-")
        elif "-" in cod or "N" in cod:
            out.append("X")
        else:
            aa = str(Seq(cod).translate())
            out.append(aa if aa != "*" else "*")
    return out


def shannon(counter):
    """Shannon entropy in bits over non-gap residues."""
    total = sum(v for k, v in counter.items() if k not in ("-",))
    if total <= 0:
        return 0.0, 0, "-", 0.0
    ent = 0.0
    for k, v in counter.items():
        if k == "-":
            continue
        p = v / total
        if p > 0:
            ent -= p * math.log2(p)
    major, mcount = max(((k, v) for k, v in counter.items() if k != "-"),
                        key=lambda kv: kv[1])
    return round(ent, 4), len([k for k in counter if k != "-"]), major, round(mcount / total, 4)


def cluster_size(header):
    m = re.search(r"cluster_size=(\d+)", header)
    return int(m.group(1)) if m else 1


def retained_columns(gene):
    """1-based UNMASKED codon columns that survived masking, from step09c's own record.

    Reconstructed from `dropped_codon_ranges` rather than inferred by matching column
    content: content matching is ambiguous wherever two columns are identical, which in
    the tandem repeats of RL2 and RS1 is common.
    """
    rep = pd.read_csv(MASK_REPORT, sep="\t")
    rep = rep[(rep.cohort == "combined") & (rep.gene == gene)]
    if len(rep) != 1:
        raise RuntimeError(
            f"{gene}: expected exactly one combined masking record, found {len(rep)}")
    row = rep.iloc[0]
    total = int(row.unmasked_codons)
    txt = row.dropped_codon_ranges
    dropped = set()
    if isinstance(txt, str) and txt.strip() and txt.strip().lower() != "nan":
        for part in txt.split(";"):
            part = part.strip()
            if not part:
                continue
            m = re.match(r"^(\d+)-(\d+)$", part)
            if m:
                dropped.update(range(int(m.group(1)), int(m.group(2)) + 1))
            else:
                dropped.add(int(part))
    return [c for c in range(1, total + 1) if c not in dropped]


def reference_proteins():
    """{(type, gene): protein} from the two reference GenBank records."""
    out = {}
    for vt, acc in REF_ACC.items():
        rec = SeqIO.read(os.path.join(DIR_REF, f"{acc}.gbk"), "genbank")
        seen = set()
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            g = (feat.qualifiers.get("gene", [""])[0]
                 or feat.qualifiers.get("locus_tag", [""])[0])
            t = feat.qualifiers.get("translation", [""])[0]
            if g and t and g not in seen:
                out[(vt, g)] = t
                seen.add(g)
    return out


REF_PROTEINS = None


def assert_reference_numbering(gene, df):
    """Abort unless every emitted residue number carries the reference's own residue.

    This is the check whose absence allowed the numbering defect to ship. It is cheap and
    it is absolute: if ref_residue says 373 then the reference protein's residue 373 must
    be the letter in ref_aa.
    """
    global REF_PROTEINS
    if REF_PROTEINS is None:
        REF_PROTEINS = reference_proteins()
    for vt in ("hsv1", "hsv2"):
        prot = REF_PROTEINS.get((vt.upper(), gene))
        if prot is None:
            log(f"WARNING: {gene}/{vt.upper()} has no reference CDS; numbering unchecked")
            continue
        sub_df = df[df[f"{vt}_ref_residue"].notna()]
        bad = []
        for res, aa in zip(sub_df[f"{vt}_ref_residue"], sub_df[f"{vt}_ref_aa"]):
            i = int(res) - 1
            if not (0 <= i < len(prot)) or prot[i] != aa:
                bad.append((int(res), aa, prot[i] if 0 <= i < len(prot) else "?"))
        if bad:
            raise RuntimeError(
                f"{gene}/{vt.upper()}: {len(bad)} residue number(s) disagree with the "
                f"reference protein; first three {bad[:3]} "
                f"(ref_residue, ref_aa, protein residue)")


def process_gene(path, gene):
    ids, seqs = read_fasta(path)
    prots = [translate_codons(s) for s in seqs]
    ncol = len(prots[0])

    # Reference numbering comes from the unmasked alignment; see the module docstring.
    unm_path = os.path.join(DIR_UNMASKED, f"{gene}_combined_aligned.fasta")
    unm_ids, unm_seqs = read_fasta(unm_path)
    if unm_ids != ids:
        raise RuntimeError(f"{gene}: masked and unmasked alignments differ in row order")
    unm_prots = [translate_codons(s) for s in unm_seqs]
    retained = retained_columns(gene)
    if len(retained) != ncol:
        raise RuntimeError(
            f"{gene}: masking record retains {len(retained)} codons but the masked "
            f"alignment has {ncol}")

    ref_rows, groups = {}, {"HSV1": [], "HSV2": []}
    for k, header in enumerate(ids):
        vt = header.split("|")[0]
        if "REF_" in header:
            ref_rows[vt] = k
        elif vt in groups:
            groups[vt].append(k)

    if "HSV1" not in ref_rows or "HSV2" not in ref_rows:
        log(f"{gene}: missing a reference row, skipped")
        return None, None

    weights = {k: cluster_size(ids[k]) for k in range(len(ids))}

    # Reference residue numbering per type, counted on the UNMASKED alignment and then
    # read back through the retained columns. Counting on `prots` (masked) silently
    # renumbers every residue after each dropped column -- the original defect.
    ref_num = {}
    for vt, k in ref_rows.items():
        unm_num, n = [], 0
        for aa in unm_prots[k]:
            if aa != "-":
                n += 1
                unm_num.append(n)
            else:
                unm_num.append(None)
        # exactness: the masked reference row must be the unmasked one restricted to the
        # retained columns, or the retained list does not describe this alignment
        rebuilt = [unm_prots[k][c - 1] for c in retained]
        if rebuilt != list(prots[k]):
            raise RuntimeError(
                f"{gene}/{vt}: retained-column reconstruction does not reproduce the "
                f"masked reference row")
        ref_num[vt] = [unm_num[c - 1] for c in retained]

    rows = []
    for c in range(ncol):
        rec = dict(gene=gene, alignment_col=c + 1,
                   hsv1_ref_residue=ref_num["HSV1"][c],
                   hsv2_ref_residue=ref_num["HSV2"][c],
                   hsv1_ref_aa=prots[ref_rows["HSV1"]][c],
                   hsv2_ref_aa=prots[ref_rows["HSV2"]][c])

        for vt in ("HSV1", "HSV2"):
            col = [prots[k][c] for k in groups[vt]]
            hap = Counter(col)
            iso = Counter()
            for k, aa in zip(groups[vt], col):
                iso[aa] += weights[k]

            e_h, n_h, maj_h, f_h = shannon(hap)
            e_i, _n_i, maj_i, f_i = shannon(iso)
            occ = sum(v for k, v in hap.items() if k != "-")

            rec.update({
                f"{vt.lower()}_n_sequences": len(col),
                f"{vt.lower()}_occupancy": round(occ / len(col), 4) if col else 0.0,
                f"{vt.lower()}_consensus_aa": maj_h,
                f"{vt.lower()}_entropy_haplotype": e_h,
                f"{vt.lower()}_entropy_isolate": e_i,
                f"{vt.lower()}_n_distinct_aa": n_h,
                f"{vt.lower()}_pct_identity_haplotype": round(100 * f_h, 2),
                f"{vt.lower()}_pct_identity_isolate": round(100 * f_i, 2),
                f"{vt.lower()}_consensus_aa_isolate": maj_i,
            })

        # Two different questions, both worth recording, and they disagree in 55 of
        # 11,801 columns. `between_type_identical` compares the two COHORT CONSENSUS
        # residues -- the population-level statement. `between_type_reference_identical`
        # compares the two REFERENCE residues, which is the statement a peptide cut out
        # of a reference protein needs; using the consensus flag for that purpose put
        # peptides into the cross-protective shortlist that do not occur in the other
        # type at all.
        a, b = rec["hsv1_consensus_aa"], rec["hsv2_consensus_aa"]
        rec["between_type_identical"] = bool(a == b and a not in ("-", "X"))
        ra, rb = rec["hsv1_ref_aa"], rec["hsv2_ref_aa"]
        rec["between_type_reference_identical"] = bool(
            ra == rb and ra not in ("-", "X"))

        # Fraction of cross-type sequence pairs sharing a residue, isolate-weighted.
        c1 = Counter()
        for k in groups["HSV1"]:
            c1[prots[k][c]] += weights[k]
        c2 = Counter()
        for k in groups["HSV2"]:
            c2[prots[k][c]] += weights[k]
        t1 = sum(v for k, v in c1.items() if k not in ("-", "X"))
        t2 = sum(v for k, v in c2.items() if k not in ("-", "X"))
        if t1 and t2:
            shared = sum(c1[k] * c2[k] for k in set(c1) & set(c2) if k not in ("-", "X"))
            rec["between_type_pct_identity"] = round(100 * shared / (t1 * t2), 2)
        else:
            rec["between_type_pct_identity"] = float("nan")

        rec["conserved_both_types"] = bool(
            rec["hsv1_pct_identity_isolate"] >= 100 * CONS_THRESHOLD
            and rec["hsv2_pct_identity_isolate"] >= 100 * CONS_THRESHOLD
            and rec["between_type_identical"]
            and rec["between_type_reference_identical"])
        rows.append(rec)

    df = pd.DataFrame(rows)
    assert_reference_numbering(gene, df)

    # Gene-level means are taken only over columns the type actually occupies. Columns
    # absent from one type (e.g. the large HSV-2-specific insertion in gG/US4, where
    # HSV-1 gG is 238 aa and HSV-2 gG is 699 aa) would otherwise be scored as 0% identity
    # for that type and would misrepresent its real conservation.
    o1 = df[df.hsv1_occupancy >= 0.5]
    o2 = df[df.hsv2_occupancy >= 0.5]
    shared = df[(df.hsv1_occupancy >= 0.5) & (df.hsv2_occupancy >= 0.5)]

    gene_level = dict(
        gene=gene,
        hsv1_n=len(groups["HSV1"]), hsv2_n=len(groups["HSV2"]),
        alignment_codons=ncol,
        codons_occupied_hsv1=len(o1), codons_occupied_hsv2=len(o2),
        codons_shared_by_both_types=len(shared),
        mean_entropy_hsv1=round(o1.hsv1_entropy_isolate.mean(), 4),
        mean_entropy_hsv2=round(o2.hsv2_entropy_isolate.mean(), 4),
        mean_pct_identity_hsv1=round(o1.hsv1_pct_identity_isolate.mean(), 2),
        mean_pct_identity_hsv2=round(o2.hsv2_pct_identity_isolate.mean(), 2),
        mean_between_type_identity=round(
            shared.between_type_pct_identity.mean(skipna=True), 2),
        pct_cols_invariant_hsv1=round(100 * (o1.hsv1_entropy_isolate == 0).mean(), 2),
        pct_cols_invariant_hsv2=round(100 * (o2.hsv2_entropy_isolate == 0).mean(), 2),
        pct_shared_cols_identical_between_types=round(
            100 * shared.between_type_identical.mean(), 2),
        pct_shared_cols_conserved_both_types=round(
            100 * shared.conserved_both_types.mean(), 2),
        pct_alignment_cols_conserved_both_types=round(
            100 * df.conserved_both_types.mean(), 2),
    )
    return df, gene_level


def conserved_windows(df, gene):
    """Contiguous WINDOW-mers conserved in both types and identical between them."""
    flag = df.conserved_both_types.tolist()
    h1 = df.hsv1_ref_residue.tolist()
    aa = df.hsv1_consensus_aa.tolist()
    out = []
    run = 0
    for i, f in enumerate(flag + [False]):
        if f:
            run += 1
            continue
        if run >= WINDOW:
            s, e = i - run, i - 1
            out.append(dict(
                gene=gene, start_col=s + 1, end_col=e + 1,
                hsv1_ref_start=h1[s], hsv1_ref_end=h1[e], length=run,
                peptide="".join(aa[s:e + 1]),
                min_hsv1_identity=round(df.hsv1_pct_identity_isolate[s:e + 1].min(), 2),
                min_hsv2_identity=round(df.hsv2_pct_identity_isolate[s:e + 1].min(), 2)))
        run = 0
    return out


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 13: per-residue conservation layer (within HSV-1, within HSV-2, between types)")

    all_rows, gene_rows, windows = [], [], []
    for path in sorted(glob.glob(os.path.join(DIR_COMB, "*_combined_masked.fasta"))):
        gene = os.path.basename(path).replace("_combined_masked.fasta", "")
        df, gl = process_gene(path, gene)
        if df is None:
            continue
        all_rows.append(df)
        gene_rows.append(gl)
        w = conserved_windows(df, gene)
        windows.extend(w)
        log(f"{gene}: {gl['alignment_codons']} cols "
            f"(HSV1 occupies {gl['codons_occupied_hsv1']}, HSV2 {gl['codons_occupied_hsv2']}, "
            f"shared {gl['codons_shared_by_both_types']}) | "
            f"identity HSV1 {gl['mean_pct_identity_hsv1']}% / HSV2 {gl['mean_pct_identity_hsv2']}% / "
            f"between {gl['mean_between_type_identity']}% | "
            f"shared cols conserved in both types {gl['pct_shared_cols_conserved_both_types']}% | "
            f"{len(w)} cross-protective window(s) >= {WINDOW}aa")

    master = pd.concat(all_rows, ignore_index=True)
    master.to_csv(os.path.join(DIR_OUT, "per_residue_conservation.tsv"), sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(
        os.path.join(DIR_OUT, "gene_level_conservation.tsv"), sep="\t", index=False)
    pd.DataFrame(windows).to_csv(
        os.path.join(DIR_OUT, "conserved_windows.tsv"), sep="\t", index=False)

    log(f"Wrote per_residue_conservation.tsv ({len(master)} residue rows), "
        f"gene_level_conservation.tsv ({len(gene_rows)} genes), "
        f"conserved_windows.tsv ({len(windows)} windows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
