"""
Step 09b (CORRECTION): Rebuild true codon-aware multiple sequence alignments.

Supersedes the step-9 half of step08_dedup_and_step09_align.py, which performed
independent PAIRWISE alignments of each sequence against the reference protein and
projected codons into each sequence's own pairwise coordinate frame. Because indel
placement was never reconciled across sequences, every alignment column downstream of
a clinical indel was out of register (see CORRECTIONS_AND_PENDING_WORK.md, A1).

Method here (MACSE-style, as specified in the proposal):
  1. Translate each QC-passed, deduplicated CDS.
  2. Build a genuine progressive protein MSA with FAMSA (Deorowicz et al., Sci Rep 2016),
     so gap placement comes from a single shared profile.
  3. Back-project nucleotide codons onto the aligned protein, preserving the reading frame.

The type reference protein is included in every alignment as `REF_<accession>|<gene>` so
that every column can be mapped to reference residue numbering for the downstream
conservation and epitope layers. Selection analyses drop it (see step11 re-run).

Outputs:
  - 06_alignments_v2/HSV1/<gene>_aligned.fasta
  - 06_alignments_v2/HSV2/<gene>_aligned.fasta
  - 06_alignments_v2/combined/<gene>_combined_aligned.fasta   (CONSERVATION USE ONLY --
    never for dN/dS or GARD; see CORRECTIONS_AND_PENDING_WORK.md, A2)
  - 06_alignments_v2/alignment_qc_v2.tsv
  - 06_alignments_v2/alignment_rebuild.log
"""

import os
import sys
import random
import itertools
from datetime import datetime

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pyfamsa import Aligner, Sequence

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_QC_PASSED = os.path.join(DIR_BASE, "05_cds_qc_passed")
DIR_OUT = os.path.join(DIR_BASE, "06_alignments_v2")
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")

REF_GBK = {"HSV1": ("NC_001806.2", os.path.join(DIR_REF, "NC_001806.2.gbk")),
           "HSV2": ("NC_001798.2", os.path.join(DIR_REF, "NC_001798.2.gbk"))}

PRIORITY_20 = ["US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53",
               "UL10", "UL49A", "RL2", "RS1", "UL54", "UL48", "UL46", "UL47",
               "UL49", "UL23", "UL30"]

LOG_PATH = os.path.join(DIR_OUT, "alignment_rebuild.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def load_reference_cds(gbk_path):
    """Return {gene: (protein, cds_nt)} for the first CDS carrying each gene name."""
    record = SeqIO.read(gbk_path, "genbank")
    full = str(record.seq).upper()
    out = {}
    for feat in record.features:
        if feat.type != "CDS":
            continue
        gene = (feat.qualifiers.get("gene", [""])[0]
                or feat.qualifiers.get("locus_tag", [""])[0])
        if not gene or gene in out:
            continue
        prot = feat.qualifiers.get("translation", [""])[0]
        if not prot:
            continue
        parts = []
        for part in feat.location.parts:
            chunk = full[part.start:part.end]
            if feat.location.strand != 1:
                chunk = str(Seq(chunk).reverse_complement())
            parts.append(chunk)
        out[gene] = (prot, "".join(parts))
    return out


def read_cds(path):
    """Read CDS records, trim to a whole number of codons, translate."""
    entries = []
    if not os.path.exists(path):
        return entries
    for rec in SeqIO.parse(path, "fasta"):
        nt = str(rec.seq).upper().replace("-", "")
        nt = nt[:len(nt) - (len(nt) % 3)]
        if not nt:
            continue
        prot = str(Seq(nt).translate(to_stop=True))
        if not prot:
            continue
        entries.append((rec.id, nt, prot))
    return entries


def famsa_protein_msa(entries):
    """Progressive protein MSA. Returns {id: aligned_protein}."""
    seqs = [Sequence(sid.encode(), prot.encode()) for sid, _nt, prot in entries]
    aligner = Aligner(guide_tree="sl", tree_heuristic=None)
    msa = aligner.align(seqs)
    return {s.id.decode(): s.sequence.decode() for s in msa}


def backproject_codons(entries, aligned_prot):
    """Project codons onto the shared aligned protein profile."""
    rows = []
    for sid, nt, prot in entries:
        ap = aligned_prot[sid]
        codons = [nt[i:i + 3] for i in range(0, len(nt), 3)]
        out, ci = [], 0
        for aa in ap:
            if aa == "-":
                out.append("---")
            elif ci < len(codons):
                out.append(codons[ci])
                ci += 1
            else:
                out.append("---")
        rows.append((sid, "".join(out)))
    return rows


def write_fasta(rows, path):
    with open(path, "w", encoding="utf-8") as fh:
        for sid, seq in rows:
            fh.write(f">{sid}\n{seq}\n")


def alignment_stats(rows, sample=60, seed=0):
    """Mean pairwise nt identity over a sample, gap fraction, non-conserved columns."""
    if len(rows) < 2:
        return dict(mean_pairwise_nt_identity=float("nan"), gap_pct=0.0,
                    cols_major_lt90_pct=0.0)
    seqs = [s for _i, s in rows]
    n, L = len(seqs), len(seqs[0])
    gaps = sum(s.count("-") for s in seqs)

    random.seed(seed)
    idx = random.sample(range(n), min(sample, n))
    tot, cnt = 0.0, 0
    for a, b in itertools.combinations(idx, 2):
        A, B = seqs[a], seqs[b]
        m = sum(1 for x, y in zip(A, B) if x == y and x != "-")
        d = sum(1 for x, y in zip(A, B) if x != "-" and y != "-")
        if d:
            tot += m / d
            cnt += 1

    bad = 0
    for i in range(L):
        col = [s[i] for s in seqs]
        nb = [c for c in col if c != "-"]
        if nb:
            top = max(nb.count(c) for c in set(nb))
            if top / len(nb) < 0.9:
                bad += 1

    return dict(mean_pairwise_nt_identity=round(100 * tot / cnt, 2) if cnt else float("nan"),
                gap_pct=round(100 * gaps / (n * L), 2),
                cols_major_lt90_pct=round(100 * bad / L, 2))


def validate(rows, label):
    """Hard invariants: rectangular, in frame, no internal stops."""
    problems = []
    lengths = {len(s) for _i, s in rows}
    if len(lengths) != 1:
        problems.append(f"ragged alignment: {len(lengths)} distinct row lengths")
    L = next(iter(lengths))
    if L % 3:
        problems.append(f"length {L} not divisible by 3")
    for sid, seq in rows:
        prot = str(Seq(seq.replace("-", "")).translate())
        if "*" in prot[:-1]:
            problems.append(f"internal stop in {sid}")
            break
    for p in problems:
        log(f"  !! {label}: {p}")
    return not problems


def main():
    for sub in ("HSV1", "HSV2", "combined"):
        os.makedirs(os.path.join(DIR_OUT, sub), exist_ok=True)

    log("=" * 78)
    log("STEP 09b: rebuilding codon-aware MSAs with FAMSA (supersedes step-9 star alignment)")

    refs = {}
    for vt, (acc, path) in REF_GBK.items():
        refs[vt] = load_reference_cds(path)
        log(f"Loaded {len(refs[vt])} reference CDS from {acc} ({vt})")

    qc_rows = []
    ok_all = True

    for gene in PRIORITY_20:
        per_type = {}

        for vt in ("HSV1", "HSV2"):
            entries = read_cds(os.path.join(DIR_QC_PASSED, vt, f"{gene}_dedup.fasta"))
            if not entries:
                log(f"{gene} {vt}: no dedup input, skipped")
                continue

            acc, _ = REF_GBK[vt]
            ref = refs[vt].get(gene)
            if ref:
                ref_prot, ref_nt = ref
                ref_nt = ref_nt[:len(ref_nt) - (len(ref_nt) % 3)]
                entries = [(f"REF_{acc}|{gene}", ref_nt, ref_prot)] + entries
            else:
                log(f"{gene} {vt}: WARNING no reference CDS found, alignment has no anchor")

            aligned_prot = famsa_protein_msa(entries)
            rows = backproject_codons(entries, aligned_prot)
            ok = validate(rows, f"{gene}/{vt}")
            ok_all = ok_all and ok

            out = os.path.join(DIR_OUT, vt, f"{gene}_aligned.fasta")
            write_fasta(rows, out)
            per_type[vt] = entries

            st = alignment_stats(rows)
            qc_rows.append(dict(gene=gene, cohort=vt, n_sequences=len(rows),
                                includes_reference=bool(ref),
                                alignment_length_nt=len(rows[0][1]),
                                alignment_length_codons=len(rows[0][1]) // 3,
                                validated=ok, **st))
            log(f"{gene} {vt}: n={len(rows)} len={len(rows[0][1])}nt "
                f"ident={st['mean_pairwise_nt_identity']}% gap={st['gap_pct']}% "
                f"var_cols={st['cols_major_lt90_pct']}%")

        # Combined alignment -- conservation / cross-protective epitope use ONLY.
        if len(per_type) == 2:
            combined = ([(f"HSV1|{i}", nt, p) for i, nt, p in per_type["HSV1"]]
                        + [(f"HSV2|{i}", nt, p) for i, nt, p in per_type["HSV2"]])
            aligned_prot = famsa_protein_msa(combined)
            rows = backproject_codons(combined, aligned_prot)
            ok = validate(rows, f"{gene}/combined")
            ok_all = ok_all and ok
            write_fasta(rows, os.path.join(DIR_OUT, "combined",
                                           f"{gene}_combined_aligned.fasta"))
            st = alignment_stats(rows)
            qc_rows.append(dict(gene=gene, cohort="combined", n_sequences=len(rows),
                                includes_reference=True,
                                alignment_length_nt=len(rows[0][1]),
                                alignment_length_codons=len(rows[0][1]) // 3,
                                validated=ok, **st))
            log(f"{gene} combined: n={len(rows)} len={len(rows[0][1])}nt "
                f"ident={st['mean_pairwise_nt_identity']}% gap={st['gap_pct']}%")

    df = pd.DataFrame(qc_rows)
    out_tsv = os.path.join(DIR_OUT, "alignment_qc_v2.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote {out_tsv} ({len(df)} alignments)")
    log(f"VALIDATION: {'ALL PASSED' if ok_all else 'FAILURES PRESENT -- see !! lines above'}")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
