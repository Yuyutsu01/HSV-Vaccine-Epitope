"""
Step 22 (PENDING WORK B8): Map conservation, known epitopes and candidates onto 3D structures.

Proposal pipeline step 10. Downloads the experimental structures that exist for the principal
HSV entry glycoproteins, assigns each chain to a reference protein by sequence (rather than by
assumption), computes per-residue relative solvent accessibility, and joins the conservation
layer, the curated IEDB epitopes and the ranked shortlist onto structural coordinates.

Structures used:
  2C36  HSV-1 gD ectodomain in complex with nectin-1
  2GUM  HSV-1 gB ectodomain
  3M1C  HSV-2 gH/gL heterodimer

Surface exposure matters for the B-cell arm: a conserved residue buried in the core is not an
antibody target, so relative solvent accessibility is reported alongside conservation rather
than conservation alone.

Outputs:
  - 17_structures/pdb/<ID>.pdb                      downloaded coordinates
  - 17_structures/chain_assignments.tsv             which chain is which reference protein
  - 17_structures/structure_residue_annotation.tsv  per residue: RSA, conservation, epitopes
  - 17_structures/conserved_exposed_patches.tsv     conserved AND exposed surface patches
  - 17_structures/figures/<ID>_<gene>.pml           PyMOL script for the paper figures
  - 17_structures/structures.log

The .pml scripts are opened in PyMOL or ChimeraX to render the figures; no local PyMOL
installation is needed to produce the tables.
"""

import os
import sys
import urllib.request
from datetime import datetime

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")
DIR_CONS = os.path.join(DIR_BASE, "11_conservation")
DIR_IEDB = os.path.join(DIR_BASE, "12_iedb")
DIR_SHORT = os.path.join(DIR_BASE, "16_shortlist")
DIR_SEL = os.path.join(DIR_BASE, "14_selection")
DIR_OUT = os.path.join(DIR_BASE, "17_structures")
DIR_PDB = os.path.join(DIR_OUT, "pdb")
DIR_FIG = os.path.join(DIR_OUT, "figures")

PDB_IDS = {"2C36": "HSV-1 gD ectodomain with nectin-1",
           "2GUM": "HSV-1 gB ectodomain",
           "3M1C": "HSV-2 gH/gL heterodimer"}

REF_ACC = {"HSV1": "NC_001806.2", "HSV2": "NC_001798.2"}
PRIORITY_20 = ["US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53",
               "UL10", "UL49A", "RL2", "RS1", "UL54", "UL48", "UL46", "UL47",
               "UL49", "UL23", "UL30"]

MIN_CHAIN_IDENTITY = 0.80      # a chain must match a reference protein this well
EXPOSED_RSA = 0.25             # conventional surface threshold
CONSERVED_IDENTITY = 99.0
PATCH_MIN = 5                  # minimum run length for a reported surface patch

# Tien et al. (2013) theoretical maximum accessible surface area, A^2.
MAX_ASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
           "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
           "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}

LOG_PATH = os.path.join(DIR_OUT, "structures.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def fetch_pdb(pdb_id):
    os.makedirs(DIR_PDB, exist_ok=True)
    path = os.path.join(DIR_PDB, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 10000:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    req = urllib.request.Request(url, headers={"User-Agent": "research-client/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(path, "wb") as fh:
        fh.write(resp.read())
    log(f"  downloaded {pdb_id} ({os.path.getsize(path)} bytes)")
    return path


def load_reference_proteins():
    out = {}
    for vt, acc in REF_ACC.items():
        rec = SeqIO.read(os.path.join(DIR_REF, f"{acc}.gbk"), "genbank")
        seen = set()
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            gene = (feat.qualifiers.get("gene", [""])[0]
                    or feat.qualifiers.get("locus_tag", [""])[0])
            if gene in PRIORITY_20 and gene not in seen:
                prot = feat.qualifiers.get("translation", [""])[0]
                if prot:
                    out[(vt, gene)] = prot
                    seen.add(gene)
    return out


def three_to_one(resname):
    return protein_letters_3to1.get(resname.capitalize(), "X")


def chain_residues(chain):
    """Ordered (residue_object, one_letter, pdb_number) for standard amino acids."""
    out = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        out.append((res, three_to_one(res.get_resname()), res.id[1]))
    return out


def make_aligner():
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    return al


def assign_chain(seq, refs, aligner):
    """Best-matching reference protein for an observed chain sequence."""
    best = (None, 0.0)
    for key, prot in refs.items():
        if not prot or len(seq) < 25:
            continue
        aln = next(iter(aligner.align(prot, seq)))
        a, b = aln[0], aln[1]
        matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
        ident = matches / len(seq)
        if ident > best[1]:
            best = (key, ident)
    return best


def map_chain_to_reference(seq, prot, aligner):
    """{index in chain sequence -> reference residue number}"""
    aln = next(iter(aligner.align(prot, seq)))
    a, b = aln[0], aln[1]
    mapping, ri, qi = {}, 0, 0
    for x, y in zip(a, b):
        if x != "-":
            ri += 1
        if y != "-":
            if x != "-":
                mapping[qi] = ri
            qi += 1
    return mapping


def load_tables():
    cons = {}
    path = os.path.join(DIR_CONS, "per_residue_conservation.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        for vt, rescol, idcol in (("HSV1", "hsv1_ref_residue", "hsv1_pct_identity_isolate"),
                                 ("HSV2", "hsv2_ref_residue", "hsv2_pct_identity_isolate")):
            sub = df[df[rescol].notna()]
            for gene, res, ident, both, ent in zip(
                    sub.gene, sub[rescol], sub[idcol], sub.conserved_both_types,
                    sub[f"{vt.lower()}_entropy_isolate"]):
                cons[(vt, gene, int(res))] = (float(ident), bool(both), float(ent))

    known = {}
    path = os.path.join(DIR_IEDB, "iedb_residue_flags.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        known = {(r.virus_type, r.gene, int(r.ref_residue)):
                 (bool(r.known_epitope), int(r.n_bcell_epitopes), int(r.n_tcell_epitopes))
                 for r in df.itertuples()}

    short = {}
    path = os.path.join(DIR_SHORT, "shortlist_primary.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        for r in df.itertuples():
            for i in range(int(r.ref_start), int(r.ref_end) + 1):
                key = (r.virus_type, r.gene, i)
                short[key] = short.get(key, 0) + 1

    selected = set()
    path = os.path.join(DIR_SEL, "selected_sites_final.tsv")
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        selected = {(r.cohort, r.gene, int(r.ref_residue)) for r in df.itertuples()}

    return cons, known, short, selected


def runs_of(flags, min_len):
    out, run = [], 0
    for i, f in enumerate(list(flags) + [False]):
        if f:
            run += 1
            continue
        if run >= min_len:
            out.append((i - run, i - 1))
        run = 0
    return out


def pymol_script(pdb_id, gene, vt, rows, path):
    """Colour by conservation, highlight known epitopes and shortlisted candidates."""
    cons_res = [r["pdb_resnum"] for r in rows
                if r["within_type_identity"] is not None
                and r["within_type_identity"] >= CONSERVED_IDENTITY]
    var_res = [r["pdb_resnum"] for r in rows
               if r["within_type_identity"] is not None
               and r["within_type_identity"] < 99.0]
    known_res = [r["pdb_resnum"] for r in rows if r["known_epitope"]]
    short_res = [r["pdb_resnum"] for r in rows if r["shortlist_peptides"]]
    both_res = [r["pdb_resnum"] for r in rows if r["conserved_both_types"]]
    sel_res = [r["pdb_resnum"] for r in rows if r["positively_selected"]]

    def sel(name, resl, chain):
        if not resl:
            return f"# {name}: none\n"
        return (f"select {name}, chain {chain} and resi "
                + "+".join(str(x) for x in sorted(set(resl))) + "\n")

    chain = rows[0]["chain"] if rows else "A"
    lines = [
        f"# {pdb_id} -- {vt} {gene}: conservation and epitope mapping",
        f"# generated by step22_structure_mapping.py on {datetime.now():%Y-%m-%d}",
        f"load ../pdb/{pdb_id}.pdb, {pdb_id}",
        "hide everything", f"show surface, {pdb_id}", "bg_color white",
        "set surface_quality, 1", "set transparency, 0.1",
        f"color grey80, {pdb_id}",
        "",
        sel("variable_sites", var_res, chain),
        sel("conserved_within_type", cons_res, chain),
        sel("conserved_both_types", both_res, chain),
        sel("known_iedb_epitopes", known_res, chain),
        sel("shortlisted_candidates", short_res, chain),
        sel("positively_selected", sel_res, chain),
        "",
        "color salmon, variable_sites",
        "color palecyan, conserved_within_type",
        "color deepteal, conserved_both_types",
        "color yellow, known_iedb_epitopes",
        "color orange, shortlisted_candidates",
        "color red, positively_selected",
        "show sticks, positively_selected",
        "",
        "set ray_opaque_background, 0",
        "orient", f"png {pdb_id}_{gene}_conservation.png, width=2400, height=1800, dpi=300, ray=1",
    ]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    os.makedirs(DIR_FIG, exist_ok=True)
    log("=" * 78)
    log("STEP 22: structural mapping of conservation, epitopes and selected sites")

    refs = load_reference_proteins()
    cons, known, short, selected = load_tables()
    log(f"Loaded {len(refs)} reference proteins, {len(cons)} conservation rows, "
        f"{len(known)} IEDB residue flags, {len(short)} shortlisted residues, "
        f"{len(selected)} selected sites"
        + ("" if selected else " (selection still running -- column will be all False)"))

    aligner = make_aligner()
    parser = PDBParser(QUIET=True)
    sr = ShrakeRupley()

    assignments, annotation, patches = [], [], []

    for pdb_id, desc in PDB_IDS.items():
        log(f"{pdb_id}: {desc}")
        try:
            path = fetch_pdb(pdb_id)
        except Exception as exc:
            log(f"  download failed: {type(exc).__name__} {exc}")
            continue

        structure = parser.get_structure(pdb_id, path)
        model = next(iter(structure))
        try:
            sr.compute(model, level="R")
        except Exception as exc:
            log(f"  SASA failed: {type(exc).__name__}; RSA will be blank")

        for chain in model:
            residues = chain_residues(chain)
            if len(residues) < 25:
                continue
            seq = "".join(r[1] for r in residues)
            (key, ident) = assign_chain(seq, refs, aligner)
            if key is None or ident < MIN_CHAIN_IDENTITY:
                log(f"  chain {chain.id}: {len(residues)} aa, no confident assignment "
                    f"(best {ident:.2f}) -- skipped")
                continue
            vt, gene = key
            assignments.append(dict(pdb_id=pdb_id, description=desc, chain=chain.id,
                                    observed_residues=len(residues), assigned_type=vt,
                                    assigned_gene=gene, sequence_identity=round(ident, 4)))
            log(f"  chain {chain.id}: {len(residues)} aa -> {vt} {gene} "
                f"(identity {ident:.3f})")

            mapping = map_chain_to_reference(seq, refs[key], aligner)
            rows = []
            for i, (res, aa, num) in enumerate(residues):
                ref_res = mapping.get(i)
                sasa = getattr(res, "sasa", None)
                rsa = (round(sasa / MAX_ASA[aa], 4)
                       if sasa is not None and aa in MAX_ASA else None)
                c = cons.get((vt, gene, ref_res)) if ref_res else None
                k = known.get((vt, gene, ref_res)) if ref_res else None
                rows.append(dict(
                    pdb_id=pdb_id, chain=chain.id, pdb_resnum=num, aa=aa,
                    virus_type=vt, gene=gene, ref_residue=ref_res,
                    sasa=round(sasa, 2) if sasa is not None else None, rsa=rsa,
                    exposed=(rsa is not None and rsa >= EXPOSED_RSA),
                    within_type_identity=c[0] if c else None,
                    conserved_both_types=bool(c[1]) if c else False,
                    entropy=c[2] if c else None,
                    known_epitope=bool(k[0]) if k else False,
                    n_bcell_epitopes=k[1] if k else 0,
                    n_tcell_epitopes=k[2] if k else 0,
                    shortlist_peptides=short.get((vt, gene, ref_res), 0) if ref_res else 0,
                    positively_selected=((vt, gene, ref_res) in selected) if ref_res else False,
                ))
            annotation.extend(rows)

            # Conserved AND exposed surface patches: the B-cell-relevant target regions.
            # Segments must be contiguous in BOTH the structure and reference numbering;
            # unresolved loops otherwise let one "patch" span a gap of tens of residues.
            segments, cur = [], []
            for r in rows:
                good = (bool(r["exposed"]) and r["within_type_identity"] is not None
                        and r["within_type_identity"] >= CONSERVED_IDENTITY
                        and r["ref_residue"] is not None)
                contiguous = (not cur) or (r["ref_residue"] == cur[-1]["ref_residue"] + 1
                                           and r["pdb_resnum"] == cur[-1]["pdb_resnum"] + 1)
                if good and contiguous:
                    cur.append(r)
                else:
                    if len(cur) >= PATCH_MIN:
                        segments.append(cur)
                    cur = [r] if good else []
            if len(cur) >= PATCH_MIN:
                segments.append(cur)

            for seg in segments:
                refs_in = [r["ref_residue"] for r in seg if r["ref_residue"]]
                patches.append(dict(
                    pdb_id=pdb_id, chain=chain.id, virus_type=vt, gene=gene,
                    pdb_start=seg[0]["pdb_resnum"], pdb_end=seg[-1]["pdb_resnum"],
                    ref_start=min(refs_in) if refs_in else None,
                    ref_end=max(refs_in) if refs_in else None,
                    length=len(seg), peptide="".join(r["aa"] for r in seg),
                    mean_rsa=round(float(np.mean([r["rsa"] for r in seg
                                                  if r["rsa"] is not None])), 3),
                    min_identity=round(min(r["within_type_identity"] for r in seg), 2),
                    conserved_both_types=sum(1 for r in seg if r["conserved_both_types"]),
                    overlaps_known_epitope=sum(1 for r in seg if r["known_epitope"]),
                ))

            pymol_script(pdb_id, gene, vt, rows,
                         os.path.join(DIR_FIG, f"{pdb_id}_{gene}.pml"))

    if not annotation:
        log("No chains could be annotated.")
        return 1

    pd.DataFrame(assignments).to_csv(
        os.path.join(DIR_OUT, "chain_assignments.tsv"), sep="\t", index=False)
    ann = pd.DataFrame(annotation)
    ann.to_csv(os.path.join(DIR_OUT, "structure_residue_annotation.tsv"),
               sep="\t", index=False)
    pat = pd.DataFrame(patches)
    pat.to_csv(os.path.join(DIR_OUT, "conserved_exposed_patches.tsv"), sep="\t", index=False)

    log(f"Wrote structure_residue_annotation.tsv ({len(ann)} residues), "
        f"chain_assignments.tsv ({len(assignments)} chains), "
        f"conserved_exposed_patches.tsv ({len(pat)} patches)")

    for (pdb_id, gene), grp in ann.groupby(["pdb_id", "gene"]):
        mapped = grp[grp.ref_residue.notna()]
        exposed = mapped[mapped.exposed]
        cons_exp = exposed[exposed.within_type_identity >= CONSERVED_IDENTITY]
        log(f"  {pdb_id} {gene}: {len(mapped)} residues mapped to reference numbering, "
            f"{len(exposed)} exposed (RSA>={EXPOSED_RSA}), "
            f"{len(cons_exp)} conserved AND exposed, "
            f"{int(mapped.known_epitope.sum())} inside a known IEDB epitope")

    if len(pat):
        log(f"Top conserved, exposed surface patches (>= {PATCH_MIN} residues):")
        for r in pat.sort_values("length", ascending=False).head(12).to_dict("records"):
            log(f"    {r['pdb_id']} {r['gene']:<5} ref {r['ref_start']}-{r['ref_end']} "
                f"len={r['length']:<3} {r['peptide'][:24]:<24} "
                f"mean RSA={r['mean_rsa']:.2f} min ident={r['min_identity']:.1f}% "
                f"known-epitope residues={r['overlaps_known_epitope']}")

    log(f"PyMOL scripts in {os.path.relpath(DIR_FIG, DIR_BASE)} -- open in PyMOL or ChimeraX "
        "to render the figures.")
    if not selected:
        log("Re-run this step after step18 so positively selected sites appear on the "
            "structures and in the .pml scripts.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
