"""
Step 14 (PENDING WORK B2): IEDB epitope download and curation.

Proposal pipeline step 9 / curation step 10. Fetches HSV-1 (NCBITaxon:10298) and
HSV-2 (NCBITaxon:10310) epitopes from the IEDB Query API, then:
  - keeps human-host records and positive assays, recording assay type, MHC class and allele
  - maps each linear peptide to the reference proteins of the 20 priority genes by EXACT
    substring match, in reference residue numbering
  - drops peptides that do not map cleanly to a priority gene
  - keeps negative and non-human records in a separate file, never as training negatives:
    absence of evidence is not a negative label (curation step 10)

Outputs:
  - 12_iedb/iedb_raw_hsv1.jsonl, iedb_raw_hsv2.jsonl       (unfiltered API payloads)
  - 12_iedb/iedb_curated_epitopes.tsv                      (human + positive + mapped)
  - 12_iedb/iedb_excluded_epitopes.tsv                     (with the reason for exclusion)
  - 12_iedb/iedb_residue_flags.tsv                         (per reference residue: covered by a known epitope?)
  - 12_iedb/iedb_curation_summary.tsv
  - 12_iedb/iedb.log
"""

import os
import sys
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

DIR_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR_REF = os.path.join(DIR_BASE, "01_reference_genomes")
DIR_OUT = os.path.join(DIR_BASE, "12_iedb")

API = "https://query-api.iedb.org/epitope_search"
# Generic client string: no personal or project identifying information is sent.
HEADERS = {"User-Agent": "research-client/1.0", "Accept": "application/json"}
PAGE = 500

TAXA = {"HSV1": "10298", "HSV2": "10310"}
REF_ACC = {"HSV1": "NC_001806.2", "HSV2": "NC_001798.2"}

PRIORITY_20 = ["US6", "UL27", "UL44", "US4", "UL22", "UL1", "US8", "US7", "UL53",
               "UL10", "UL49A", "RL2", "RS1", "UL54", "UL48", "UL46", "UL47",
               "UL49", "UL23", "UL30"]

MIN_PEPTIDE = 8
LOG_PATH = os.path.join(DIR_OUT, "iedb.log")


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    os.makedirs(DIR_OUT, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def api_get(params, prefer=None):
    url = API + "?" + urllib.parse.urlencode(params, safe="{}:,.*")
    headers = dict(HEADERS)
    if prefer:
        headers["Prefer"] = prefer
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp, json.loads(resp.read().decode())


def fetch_all(taxon, out_path):
    """Page through every epitope for a taxon, caching to JSONL."""
    if os.path.exists(out_path):
        with open(out_path, encoding="utf-8") as fh:
            recs = [json.loads(l) for l in fh if l.strip()]
        log(f"taxon {taxon}: {len(recs)} records from cache {os.path.basename(out_path)}")
        return recs

    flt = "cs.{NCBITaxon:%s}" % taxon
    resp, _ = api_get({"source_organism_iri_search": flt, "limit": 1,
                       "order": "structure_id"}, prefer="count=exact")
    total = int(resp.headers.get("content-range", "0-0/0").split("/")[-1])
    log(f"taxon {taxon}: {total} epitope records to fetch")

    recs = []
    offset = 0
    while offset < total:
        for attempt in range(4):
            try:
                _r, page = api_get({"source_organism_iri_search": flt,
                                    "limit": PAGE, "offset": offset,
                                    "order": "structure_id"})
                break
            except Exception as exc:                       # transient API errors
                if attempt == 3:
                    raise
                log(f"  retry {attempt + 1} at offset {offset}: {type(exc).__name__}")
                time.sleep(3 * (attempt + 1))
        if not page:
            break
        recs.extend(page)
        offset += len(page)
        log(f"  fetched {len(recs)}/{total}")
        time.sleep(0.3)

    with open(out_path, "w", encoding="utf-8") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    return recs


def load_reference_proteins():
    """{(type, gene): protein} for the priority genes."""
    out = {}
    for vt, acc in REF_ACC.items():
        record = SeqIO.read(os.path.join(DIR_REF, f"{acc}.gbk"), "genbank")
        seen = set()
        for feat in record.features:
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


def first(lst):
    if isinstance(lst, list):
        return "; ".join(str(x) for x in lst if x is not None) if lst else ""
    return "" if lst is None else str(lst)


def curate(recs, vt, refs):
    kept, excluded = [], []
    for r in recs:
        pep = (r.get("linear_sequence") or "").strip().upper()
        hosts = first(r.get("host_organism_names"))
        quals = r.get("qualitative_measures") or []
        assays = first(r.get("assay_names"))
        alleles = first(r.get("mhc_allele_names"))
        classes = first(r.get("mhc_classes"))

        base = dict(virus_type=vt, structure_id=r.get("structure_id"),
                    iedb_iri=r.get("structure_iri"), peptide=pep,
                    peptide_length=len(pep), host=hosts,
                    qualitative_measures=first(quals), assay_names=assays,
                    mhc_classes=classes, mhc_alleles=alleles,
                    source_antigen=first(r.get("parent_source_antigen_names")),
                    source_organism=first(r.get("source_organism_names")),
                    pubmed_ids=first(r.get("pubmed_ids")),
                    n_tcell_assays=len(r.get("tcell_ids") or []),
                    n_bcell_assays=len(r.get("bcell_ids") or []),
                    n_mhc_ligand_assays=len(r.get("elution_ids") or []))

        if not pep or len(pep) < MIN_PEPTIDE:
            excluded.append({**base, "exclusion_reason":
                             "no linear sequence" if not pep
                             else f"peptide shorter than {MIN_PEPTIDE} aa"})
            continue
        if "Homo sapiens" not in hosts:
            excluded.append({**base, "exclusion_reason": "non-human host"})
            continue
        if not any("Positive" in str(q) for q in quals):
            excluded.append({**base, "exclusion_reason": "no positive assay"})
            continue

        hits = []
        for gene in PRIORITY_20:
            prot = refs.get((vt, gene))
            if not prot:
                continue
            pos = prot.find(pep)
            while pos != -1:
                hits.append((gene, pos + 1, pos + len(pep)))
                pos = prot.find(pep, pos + 1)

        if not hits:
            excluded.append({**base, "exclusion_reason":
                             "no exact match in a priority-gene reference protein"})
            continue

        for gene, start, end in hits:
            kept.append({**base, "gene": gene, "ref_start": start, "ref_end": end,
                         "n_reference_matches": len(hits),
                         "epitope_class": ("T-cell" if base["n_tcell_assays"] else
                                           "B-cell" if base["n_bcell_assays"] else
                                           "MHC-ligand")})
    return kept, excluded


def residue_flags(kept, refs):
    """Per reference residue: is it inside a known, positively assayed human epitope?"""
    cov = defaultdict(lambda: {"T-cell": 0, "B-cell": 0, "MHC-ligand": 0})
    for row in kept:
        key = (row["virus_type"], row["gene"])
        for res in range(row["ref_start"], row["ref_end"] + 1):
            cov[(key, res)][row["epitope_class"]] += 1

    out = []
    for (vt, gene), prot in sorted(refs.items()):
        for res in range(1, len(prot) + 1):
            c = cov.get(((vt, gene), res))
            out.append(dict(virus_type=vt, gene=gene, ref_residue=res, ref_aa=prot[res - 1],
                            known_epitope=bool(c),
                            n_tcell_epitopes=c["T-cell"] if c else 0,
                            n_bcell_epitopes=c["B-cell"] if c else 0,
                            n_mhc_ligand_epitopes=c["MHC-ligand"] if c else 0))
    return out


def main():
    os.makedirs(DIR_OUT, exist_ok=True)
    log("=" * 78)
    log("STEP 14: IEDB epitope download and curation")

    refs = load_reference_proteins()
    log(f"Loaded {len(refs)} reference proteins for {len(PRIORITY_20)} priority genes")
    missing = [(vt, g) for vt in REF_ACC for g in PRIORITY_20 if (vt, g) not in refs]
    if missing:
        log(f"  note: no reference CDS for {missing}")

    all_kept, all_excl, summary = [], [], []
    for vt, taxon in TAXA.items():
        recs = fetch_all(taxon, os.path.join(DIR_OUT, f"iedb_raw_{vt.lower()}.jsonl"))
        kept, excl = curate(recs, vt, refs)
        all_kept.extend(kept)
        all_excl.extend(excl)

        reasons = pd.Series([e["exclusion_reason"] for e in excl]).value_counts().to_dict()
        summary.append(dict(virus_type=vt, taxon=taxon, records_fetched=len(recs),
                            curated_mappings=len(kept),
                            unique_peptides_kept=len({k["peptide"] for k in kept}),
                            excluded=len(excl), **reasons))
        log(f"{vt}: {len(recs)} fetched -> {len({k['peptide'] for k in kept})} unique "
            f"curated peptides ({len(kept)} gene mappings), {len(excl)} excluded")
        for reason, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
            log(f"    excluded {n}: {reason}")

    pd.DataFrame(all_kept).to_csv(
        os.path.join(DIR_OUT, "iedb_curated_epitopes.tsv"), sep="\t", index=False)
    pd.DataFrame(all_excl).to_csv(
        os.path.join(DIR_OUT, "iedb_excluded_epitopes.tsv"), sep="\t", index=False)
    pd.DataFrame(residue_flags(all_kept, refs)).to_csv(
        os.path.join(DIR_OUT, "iedb_residue_flags.tsv"), sep="\t", index=False)
    pd.DataFrame(summary).to_csv(
        os.path.join(DIR_OUT, "iedb_curation_summary.tsv"), sep="\t", index=False)

    if all_kept:
        df = pd.DataFrame(all_kept)
        log("Curated epitopes per gene (unique peptides):")
        for (vt, gene), grp in df.groupby(["virus_type", "gene"]):
            log(f"    {vt} {gene}: {grp.peptide.nunique()} "
                f"(T-cell {grp[grp.epitope_class == 'T-cell'].peptide.nunique()}, "
                f"B-cell {grp[grp.epitope_class == 'B-cell'].peptide.nunique()}, "
                f"MHC-ligand {grp[grp.epitope_class == 'MHC-ligand'].peptide.nunique()})")
    log("Note: excluded records are NOT negatives. Untested regions have no label.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
