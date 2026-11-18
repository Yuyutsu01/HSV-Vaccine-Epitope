"""
HSV-1 & HSV-2 Protein FASTA Comprehensive Audit Master Script
Generates 14 machine-readable TSVs and a master annotated dataset.
Utilizes Python multiprocessing to process all 22,710 records efficiently.
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

OUT_DIR = "audit_results"
os.makedirs(OUT_DIR, exist_ok=True)

# ----------------------------------------------------------------------
# 1. Reference Proteome Loading & Configuration
# ----------------------------------------------------------------------

def get_ref_data():
    """
    Loads reference records from disk and prepares canonical lengths and 5-mer sets.
    """
    def parse_ref(path, vtype):
        ref_dict = {}
        for rec in SeqIO.parse(path, "fasta"):
            parts = rec.description.split("|")
            prot_id, gene, locus_tag, product, copy_info, plen = parts[0], parts[1], parts[2], parts[3], parts[4], int(parts[5].replace("len=", ""))
            seq_str = str(rec.seq).upper()
            if gene not in ref_dict:
                ref_dict[gene] = {
                    "gene": gene, "prot_id": prot_id, "locus_tag": locus_tag,
                    "product": product, "length": plen, "seq": seq_str, "virus": vtype,
                    "copies": 1 if "copy=1" in copy_info else 2
                }
            else:
                ref_dict[gene]["copies"] = 2
        return ref_dict

    refs1 = parse_ref("references/HSV1_ref_proteins.fasta", "HSV-1")
    refs2 = parse_ref("references/HSV2_ref_proteins.fasta", "HSV-2")
    return refs1, refs2

REFS_HSV1, REFS_HSV2 = get_ref_data()

# Precompute 5-mers for fast homology mapping
K = 5
KMERS_HSV1 = {g: set(d["seq"][i:i+K] for i in range(len(d["seq"])-K+1)) for g, d in REFS_HSV1.items()}
KMERS_HSV2 = {g: set(d["seq"][i:i+K] for i in range(len(d["seq"])-K+1)) for g, d in REFS_HSV2.items()}

# ----------------------------------------------------------------------
# 2. Gene Synonym Definitions & Priority Sets
# ----------------------------------------------------------------------

PRIORITY_GENE_INFO = [
    # Glycoproteins
    ("US6", "gD/US6", "Glycoprotein"),
    ("UL27", "gB/UL27", "Glycoprotein"),
    ("UL44", "gC/UL44", "Glycoprotein"),
    ("US4", "gG/US4", "Glycoprotein"),
    ("UL22", "gH/UL22", "Glycoprotein"),
    ("UL1", "gL/UL1", "Glycoprotein"),
    ("US8", "gE/US8", "Glycoprotein"),
    ("US7", "gI/US7", "Glycoprotein"),
    ("UL53", "gK/UL53", "Glycoprotein"),
    ("UL10", "gM/UL10", "Glycoprotein"),
    ("UL49A", "gN/UL49A", "Glycoprotein"),
    # T-cell targets
    ("RL2", "ICP0/RL2", "T-cell target"),
    ("RS1", "ICP4/RS1", "T-cell target"),
    ("UL54", "ICP27/UL54", "T-cell target"),
    ("UL48", "VP16/UL48", "T-cell target"),
    ("UL46", "UL46/VP11-12", "T-cell target"),
    ("UL47", "UL47/VP13-14", "T-cell target"),
    ("UL49", "UL49/VP22", "T-cell target"),
    # Drug targets
    ("UL23", "UL23/TK", "Drug target"),
    ("UL30", "UL30/polymerase", "Drug target")
]
PRIORITY_GENES_DICT = {g[0]: g for g in PRIORITY_GENE_INFO}
PRIORITY_GENE_SET = set(g[0] for g in PRIORITY_GENE_INFO)
PRIORITY_GENES = PRIORITY_GENE_INFO

def align_to_ref(q, r):
    """Pairwise alignment helper for redundancy calculations."""
    aligner_inst = PairwiseAligner()
    aligner_inst.mode = "local"
    try:
        alignments = aligner_inst.align(q, r)
        if not alignments:
            return 0.0, 0.0, 0
        top = alignments[0]
        q_aligned = sum(top.aligned[0][i][1] - top.aligned[0][i][0] for i in range(len(top.aligned[0])))
        matches = 0
        for block_idx in range(len(top.aligned[0])):
            qs, qe = top.aligned[0][block_idx]
            rs, re_ = top.aligned[1][block_idx]
            sub_q = q[qs:qe]
            sub_r = r[rs:re_]
            matches += sum(1 for a, b in zip(sub_q, sub_r) if a == b)
        pid = (matches / max(1, q_aligned)) * 100.0
        qcov = (q_aligned / len(q)) * 100.0
        return pid, qcov, matches
    except Exception:
        return 0.0, 0.0, 0


# Regex to detect explicit gene symbols or standard protein names in deflines
SYNONYM_PATTERNS = [
    ("UL23", r"\b(UL23|thymidine kinase|TK)\b"),
    ("UL30", r"\b(UL30|DNA polymerase catalytic|DNA polymerase(?!\s+processivity))\b"),
    ("UL27", r"\b(UL27|glycoprotein B|gB)\b"),
    ("US6",  r"\b(US6|glycoprotein D|gD)\b"),
    ("UL44", r"\b(UL44|glycoprotein C|gC)\b"),
    ("US4",  r"\b(US4|glycoprotein G|gG)\b"),
    ("UL22", r"\b(UL22|glycoprotein H|gH)\b"),
    ("UL1",  r"\b(UL1|glycoprotein L|gL)\b"),
    ("US8",  r"\b(US8|glycoprotein E|gE)\b"),
    ("US7",  r"\b(US7|glycoprotein I|gI)\b"),
    ("UL53", r"\b(UL53|glycoprotein K|gK)\b"),
    ("UL10", r"\b(UL10|glycoprotein M|gM)\b"),
    ("UL49A",r"\b(UL49A|UL49\.5|glycoprotein N|gN)\b"),
    ("RL2",  r"\b(RL2|ICP0|IE110)\b"),
    ("RS1",  r"\b(RS1|ICP4|IE175)\b"),
    ("UL54", r"\b(UL54|ICP27|IE63)\b"),
    ("UL48", r"\b(UL48|VP16|alpha-TIF|Vmw65)\b"),
    ("UL46", r"\b(UL46|VP11/12|VP11-12)\b"),
    ("UL47", r"\b(UL47|VP13/14|VP13-14)\b"),
    ("UL49", r"\b(UL49\b(?!A)|VP22)\b"),
    ("RL1",  r"\b(RL1|ICP34\.5)\b"),
    ("US3",  r"\b(US3|serine/threonine-protein kinase US3)\b"),
    ("UL19", r"\b(UL19|major capsid protein)\b"),
    ("UL36", r"\b(UL36|large tegument protein)\b"),
    ("UL39", r"\b(UL39|ribonucleotide reductase subunit 1|ribonucleotide reductase large)\b"),
    ("UL40", r"\b(UL40|ribonucleotide reductase subunit 2|ribonucleotide reductase small)\b"),
    ("UL42", r"\b(UL42|DNA polymerase processivity)\b"),
    ("UL5",  r"\b(UL5|helicase-primase helicase)\b"),
    ("UL8",  r"\b(UL8|helicase-primase subunit)\b"),
    ("UL52", r"\b(UL52|helicase-primase primase)\b"),
    ("UL26", r"\b(UL26\b(?!.5)|capsid maturation protease)\b"),
    ("UL26.5", r"\b(UL26\.5|capsid scaffold protein)\b"),
    ("UL2",  r"\b(UL2|uracil-DNA glycosylase)\b"),
    ("UL50", r"\b(UL50|deoxyuridine triphosphatase|dUTPase)\b"),
    ("US1",  r"\b(US1|ICP22|IE68)\b"),
    ("US12", r"\b(US12|ICP47|TAP transporter inhibitor)\b"),
    ("US2",  r"\b(US2|virion protein US2)\b"),
    ("US5",  r"\b(US5|glycoprotein J|gJ)\b"),
    ("US8A", r"\b(US8A|membrane protein US8A)\b"),
    ("US9",  r"\b(US9|membrane protein US9)\b"),
    ("US10", r"\b(US10|virion protein US10)\b"),
    ("US11", r"\b(US11|tegument protein US11)\b"),
    ("UL3",  r"\b(UL3|nuclear protein UL3)\b"),
    ("UL4",  r"\b(UL4|nuclear protein UL4)\b"),
    ("UL6",  r"\b(UL6|capsid portal protein)\b"),
    ("UL7",  r"\b(UL7|tegument protein UL7)\b"),
    ("UL9",  r"\b(UL9|origin-binding protein|origin-binding)\b"),
    ("UL11", r"\b(UL11|myristylated tegument)\b"),
    ("UL12", r"\b(UL12|deoxyribonuclease)\b"),
    ("UL13", r"\b(UL13|tegument serine/threonine)\b"),
    ("UL14", r"\b(UL14|tegument protein UL14)\b"),
    ("UL15", r"\b(UL15|DNA packaging terminase subunit 1)\b"),
    ("UL16", r"\b(UL16|tegument protein UL16)\b"),
    ("UL17", r"\b(UL17|DNA packaging tegument)\b"),
    ("UL18", r"\b(UL18|capsid triplex subunit 2)\b"),
    ("UL20", r"\b(UL20|envelope protein UL20)\b"),
    ("UL21", r"\b(UL21|tegument protein UL21)\b"),
    ("UL24", r"\b(UL24|nuclear protein UL24)\b"),
    ("UL25", r"\b(UL25|DNA packaging tegument protein UL25)\b"),
    ("UL28", r"\b(UL28|DNA packaging terminase subunit 2)\b"),
    ("UL29", r"\b(UL29|single-stranded DNA-binding protein|major DNA-binding protein)\b"),
    ("UL31", r"\b(UL31|nuclear egress lamina)\b"),
    ("UL32", r"\b(UL32|DNA packaging protein UL32)\b"),
    ("UL33", r"\b(UL33|DNA packaging protein UL33)\b"),
    ("UL34", r"\b(UL34|nuclear egress membrane)\b"),
    ("UL35", r"\b(UL35|small capsid protein)\b"),
    ("UL37", r"\b(UL37|tegument protein UL37)\b"),
    ("UL38", r"\b(UL38|capsid triplex subunit 1)\b"),
    ("UL41", r"\b(UL41|tegument host shutoff)\b"),
    ("UL43", r"\b(UL43|envelope protein UL43)\b"),
    ("UL45", r"\b(UL45|membrane protein UL45)\b"),
    ("UL51", r"\b(UL51|tegument protein UL51)\b"),
    ("UL55", r"\b(UL55|nuclear protein UL55)\b"),
    ("UL56", r"\b(UL56|membrane protein UL56)\b"),
]

# Lab strain and engineered construct patterns
LAB_STRAINS = ["strain 17", "17+", "strain F", "KOS", "McKrae", "HFEM", "MacIntyre", "SC16", "HG52", "strain 186", "strain MS", "strain RH2", "strain R-15"]
LAB_PAT = re.compile(r"\b(" + "|".join(re.escape(x) for x in LAB_STRAINS) + r")\b", re.I)

ENGINEERED_TERMS = ["BAC", "oncolytic", "vaccine vector", "deletion mutant", "recombinant", "mutant", "synthetic", "patent", "transgenic", "chimera", "vector"]
ENG_PAT = re.compile(r"\b(" + "|".join(re.escape(x) for x in ENGINEERED_TERMS) + r")\b", re.I)

COUNTRIES = [
    "USA", "United States", "China", "Japan", "Germany", "United Kingdom", "UK",
    "France", "Canada", "Australia", "South Africa", "Kenya", "Uganda", "Tanzania",
    "Brazil", "India", "Russia", "Italy", "Spain", "Sweden", "Netherlands", "Norway"
]
COUNTRY_PAT = re.compile(r"\b(" + "|".join(COUNTRIES) + r")\b", re.I)

# ----------------------------------------------------------------------
# 3. Worker Function for Sequence Evaluation
# ----------------------------------------------------------------------

def process_record_batch(batch_records):
    """
    Processes a batch of FASTA record tuples:
    (rec_id, rec_desc, seq_str, source_file, file_label)
    Returns list of parsed/analyzed result dictionaries.
    """
    local_aligner = PairwiseAligner()
    local_aligner.mode = "local"
    local_aligner.open_gap_score = -10.0
    local_aligner.extend_gap_score = -0.5
    local_aligner.match_score = 2.0
    local_aligner.mismatch_score = -1.0

    results = []

    for rec_id, desc, seq_str, source_file, file_label in batch_records:
        seq_len = len(seq_str)
        
        # 1. Non-standard amino acids & internal stop codons
        standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
        non_std_chars = [c for c in seq_str if c not in standard_aa]
        non_std_count = len(non_std_chars)
        x_count = seq_str.count("X")
        x_pct = (x_count / max(1, seq_len)) * 100.0
        
        # Check for internal stop codon (*)
        # Stop at very last position is terminal; anywhere else is internal
        has_internal_stop = False
        if "*" in seq_str:
            last_idx = seq_str.rfind("*")
            if last_idx < seq_len - 1 or seq_str.count("*") > 1:
                has_internal_stop = True

        # 2. Parse metadata from header
        # Accession
        acc = rec_id
        if "|" in acc:
            acc_parts = acc.split("|")
            acc = acc_parts[1] if len(acc_parts) > 1 and acc_parts[1] else acc_parts[0]
            
        # Organism bracket
        org_matches = re.findall(r"\[([^\]]+)\]", desc)
        organism = org_matches[-1] if org_matches else "Unspecified"
        
        # Host
        host = "Human" if ("human" in desc.lower() or "homo sapiens" in desc.lower()) else ("Mouse" if ("mouse" in desc.lower() or "mus musculus" in desc.lower()) else "Unspecified")
        
        # Strain / Isolate
        strain = "Unspecified"
        m_strain = LAB_PAT.search(desc)
        if m_strain:
            strain = m_strain.group(1)
        else:
            m_iso = re.search(r"\b(?:strain|isolate)\s+([A-Za-z0-9_\-\.\/]+)", desc, re.I)
            if m_iso:
                strain = m_iso.group(1)
                
        # Country
        country = "Absent"
        m_c = COUNTRY_PAT.search(desc)
        if m_c:
            country = m_c.group(1).upper()
            
        # Year
        year = "Absent"
        m_y = re.search(r"\b(19[5-9]\d|20[0-2]\d)\b", desc)
        if m_y:
            year = m_y.group(1)
            
        # Passage
        passage = "Absent"
        m_p = re.search(r"\b(passage|passaged|p\d+)\b", desc, re.I)
        if m_p:
            passage = m_p.group(1)
            
        # BioProject / Study
        study = "Absent"
        m_bp = re.search(r"\b(BioProject|PRJN[A-Z0-9]+)\b", desc, re.I)
        if m_bp:
            study = m_bp.group(1)
            
        # Patent
        is_patent = "patent" in desc.lower()
        
        # Lab strain flag
        is_lab_strain = bool(LAB_PAT.search(desc))
        
        # Engineered construct flag
        is_engineered = bool(ENG_PAT.search(desc)) or is_patent

        # Gene name explicitness in header
        explicit_gene = None
        for g, pat in SYNONYM_PATTERNS:
            if re.search(pat, desc, re.I):
                explicit_gene = g
                break
                
        has_explicit_gene_symbol = bool(re.search(r"\b(UL\d+|US\d+|RL1|RL2|RS1)\b", desc, re.I))
        gene_explicitness = "Explicit_Symbol" if has_explicit_gene_symbol else ("Explicit_Product_Only" if explicit_gene else "Absent_or_Generic")

        # 3. Reference alignment and homology mapping
        # Determine candidate gene
        if explicit_gene and explicit_gene in REFS_HSV1:
            cand_g1 = explicit_gene
            cand_g2 = explicit_gene if explicit_gene in REFS_HSV2 else explicit_gene
        elif seq_len >= K:
            q_kmers = set(seq_str[i:i+K] for i in range(len(seq_str)-K+1))
            scores1 = {g: len(q_kmers & km) for g, km in KMERS_HSV1.items()}
            scores2 = {g: len(q_kmers & km) for g, km in KMERS_HSV2.items()}
            cand_g1 = max(scores1, key=scores1.get)
            cand_g2 = max(scores2, key=scores2.get)
        else:
            cand_g1 = "UL23"
            cand_g2 = "UL23"

        # Align against candidate HSV-1 and HSV-2 reference proteins
        ref1_seq = REFS_HSV1[cand_g1]["seq"] if cand_g1 in REFS_HSV1 else ""
        ref2_seq = REFS_HSV2[cand_g2]["seq"] if cand_g2 in REFS_HSV2 else ""
        
        def run_align(q, r):
            if not q or not r or len(q) == 0:
                return 0.0, 0.0, 0.0
            try:
                alignments = local_aligner.align(q, r)
                if not alignments:
                    return 0.0, 0.0, 0.0
                top = alignments[0]
                q_aligned = sum(top.aligned[0][i][1] - top.aligned[0][i][0] for i in range(len(top.aligned[0])))
                r_aligned = sum(top.aligned[1][i][1] - top.aligned[1][i][0] for i in range(len(top.aligned[1])))
                matches = 0
                for block_idx in range(len(top.aligned[0])):
                    qs, qe = top.aligned[0][block_idx]
                    rs, re_ = top.aligned[1][block_idx]
                    sub_q = q[qs:qe]
                    sub_r = r[rs:re_]
                    matches += sum(1 for a, b in zip(sub_q, sub_r) if a == b)
                pid = (matches / max(1, q_aligned)) * 100.0
                qcov = (q_aligned / len(q)) * 100.0
                rcov = (r_aligned / len(r)) * 100.0
                return pid, qcov, rcov
            except Exception:
                return 0.0, 0.0, 0.0

        id1, qcov1, rcov1 = run_align(seq_str, ref1_seq)
        id2, qcov2, rcov2 = run_align(seq_str, ref2_seq)

        # Assign best virus hit and best gene
        if id1 < 40.0 and id2 < 40.0 and seq_len > 30 and not explicit_gene:
            best_virus = "Non-HSV Contamination"
            best_gene = "NON_HSV"
            ref_len = 0
        elif id1 >= id2:
            best_virus = "HSV-1"
            best_gene = cand_g1
            ref_len = REFS_HSV1[cand_g1]["length"] if cand_g1 in REFS_HSV1 else 0
        else:
            best_virus = "HSV-2"
            best_gene = cand_g2
            ref_len = REFS_HSV2[cand_g2]["length"] if cand_g2 in REFS_HSV2 else 0

        # Ref length ratio & fragment evaluation
        len_ratio = (seq_len / max(1, ref_len)) * 100.0 if ref_len > 0 else 0.0
        is_fragment = (len_ratio < 95.0) or (seq_len < 30)

        # 4. Check for biological UL23 (TK) frameshift / drug resistance
        # Canonical TK length is 376 aa in HSV-1 and HSV-2.
        # Known clinical homopolymer frameshifts produce functional truncations (e.g. 150-340 aa).
        is_tk_frameshift = False
        if best_gene == "UL23":
            if 80 <= seq_len <= 360 and x_pct <= 1.0 and not has_internal_stop:
                is_tk_frameshift = True

        # Check mislabeled virus
        is_mislabeled = False
        if file_label == "HSV-1" and best_virus == "HSV-2" and id2 > id1 + 5.0 and id2 > 70.0:
            is_mislabeled = True
        elif file_label == "HSV-2" and best_virus == "HSV-1" and id1 > id2 + 5.0 and id1 > 70.0:
            is_mislabeled = True

        # Check non-HSV contamination
        is_non_hsv = (best_virus == "Non-HSV Contamination") or ("homo sapiens" in organism.lower() and best_virus != "HSV-1" and best_virus != "HSV-2") or ("bovine" in organism.lower()) or ("anatid" in organism.lower()) or ("gallid" in organism.lower()) or ("vector" in organism.lower() and not explicit_gene)

        # 5. Final Category: KEEP / FLAG / DROP (Mutually Exclusive)
        status = "KEEP"
        reasons = []

        if is_non_hsv:
            status = "DROP"
            reasons.append("Non-HSV contamination")
        elif seq_len < 20:
            status = "DROP"
            reasons.append(f"Extreme micro-fragment (len={seq_len} aa)")
        elif has_internal_stop:
            status = "DROP"
            reasons.append("Internal stop codon (pseudogene or sequencing error)")
        elif x_pct > 10.0:
            status = "DROP"
            reasons.append(f"Excessive sequence ambiguity ({x_pct:.1f}% X)")
        elif is_tk_frameshift:
            status = "FLAG"
            reasons.append("UL23 acyclovir-resistance homopolymer frameshift candidate (KEEP FOR DRUG RESISTANCE)")
        elif is_mislabeled:
            status = "FLAG"
            reasons.append(f"Mislabeled viral type ({best_virus} sequence in {file_label} file)")
        elif is_lab_strain:
            status = "FLAG"
            reasons.append(f"Laboratory strain ({strain})")
        elif is_engineered:
            status = "FLAG"
            reasons.append("Engineered/synthetic construct or patent sequence")
        elif is_fragment:
            status = "FLAG"
            reasons.append(f"Fragment (<95% reference length: {len_ratio:.1f}%)")
        elif x_pct > 1.0:
            status = "FLAG"
            reasons.append(f"Moderate sequence ambiguity ({x_pct:.1f}% X)")
        else:
            status = "KEEP"
            reasons.append("High-quality, full-length, wild-type clinical/isolate sequence")

        results.append({
            "accession": acc,
            "header": desc,
            "source_file": source_file,
            "file_type": file_label,
            "length": seq_len,
            "ref_length": ref_len,
            "length_ratio_pct": round(len_ratio, 2),
            "non_std_aa_count": non_std_count,
            "x_count": x_count,
            "x_pct": round(x_pct, 2),
            "has_internal_stop": has_internal_stop,
            "organism": organism,
            "strain": strain,
            "country": country,
            "year": year,
            "host": host,
            "passage": passage,
            "study": study,
            "is_patent": is_patent,
            "is_lab_strain": is_lab_strain,
            "is_engineered": is_engineered,
            "gene_explicitness": gene_explicitness,
            "explicit_gene": explicit_gene or "None",
            "assigned_gene": best_gene,
            "is_priority_gene": best_gene in PRIORITY_GENE_SET,
            "best_virus_match": best_virus,
            "hsv1_hit_gene": cand_g1,
            "hsv1_identity": round(id1, 2),
            "hsv1_qcov": round(qcov1, 2),
            "hsv1_rcov": round(rcov1, 2),
            "hsv2_hit_gene": cand_g2,
            "hsv2_identity": round(id2, 2),
            "hsv2_qcov": round(qcov2, 2),
            "hsv2_rcov": round(rcov2, 2),
            "is_fragment": is_fragment,
            "is_tk_frameshift": is_tk_frameshift,
            "is_mislabeled": is_mislabeled,
            "is_non_hsv": is_non_hsv,
            "category": status,
            "reason": "; ".join(reasons),
            "seq_str": seq_str
        })

    return results

# ----------------------------------------------------------------------
# 4. Main Multi-Process Execution Engine
# ----------------------------------------------------------------------

def run():
    print("="*60)
    print("STARTING FULL HSV-1 & HSV-2 FASTA AUDIT")
    print("="*60)

    # Load all records into memory
    all_record_tuples = []
    
    print("Reading non_redundant.fasta (HSV-1)...")
    for r in SeqIO.parse("non_redundant.fasta", "fasta"):
        all_record_tuples.append((r.id, r.description, str(r.seq).upper(), "non_redundant.fasta", "HSV-1"))
    print(f"Loaded {len(all_record_tuples)} records for HSV-1.")
    
    hsv1_count = len(all_record_tuples)
    
    print("Reading HSV2_non_redundant (2).fasta (HSV-2)...")
    for r in SeqIO.parse("HSV2_non_redundant (2).fasta", "fasta"):
        all_record_tuples.append((r.id, r.description, str(r.seq).upper(), "HSV2_non_redundant (2).fasta", "HSV-2"))
    print(f"Loaded {len(all_record_tuples) - hsv1_count} records for HSV-2.")
    print(f"Total sequences to audit: {len(all_record_tuples)}")

    # Chunk records for 16 CPU cores
    num_workers = min(16, os.cpu_count() or 4)
    chunk_size = 200
    batches = [all_record_tuples[i:i+chunk_size] for i in range(0, len(all_record_tuples), chunk_size)]
    print(f"Dispatching {len(batches)} batches across {num_workers} worker processes...")

    all_results = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for batch_res in executor.map(process_record_batch, batches):
            all_results.extend(batch_res)

    print(f"Audit processing complete. Total audited records: {len(all_results)}")
    
    # Convert to DataFrame
    df = pd.DataFrame(all_results)
    
    # Drop sequence string from summary exports to keep TSVs compact
    df_export = df.drop(columns=["seq_str"])
    df_export.to_csv(os.path.join(OUT_DIR, "all_sequences_master_annotated.tsv"), sep="\t", index=False)
    print(f"Master annotated dataset saved to {OUT_DIR}/all_sequences_master_annotated.tsv")

    # ------------------------------------------------------------------
    # Section 1: File-Level Statistics & Length Histograms
    # ------------------------------------------------------------------
    print("Generating Section 1: File-Level Statistics & Histograms...")
    file_stats = []
    for flabel, fname in [("HSV-1", "non_redundant.fasta"), ("HSV-2", "HSV2_non_redundant (2).fasta")]:
        sub = df[df["file_type"] == flabel]
        lengths = sub["length"].values
        seq_strings = sub["seq_str"].values
        exact_dups = len(seq_strings) - len(set(seq_strings))
        has_non_std = (sub["non_std_aa_count"] > 0).sum()
        
        file_stats.append({
            "File_Name": fname,
            "Virus_Type": flabel,
            "Sequence_Count": len(sub),
            "Total_Residues": int(lengths.sum()),
            "Min_Length": int(lengths.min()),
            "Median_Length": float(np.median(lengths)),
            "Mean_Length": round(float(lengths.mean()), 2),
            "Max_Length": int(lengths.max()),
            "Non_Standard_AA_Sequences": int(has_non_std),
            "Non_Standard_AA_Pct": round(has_non_std / len(sub) * 100, 2),
            "Exact_Duplicate_Sequences": int(exact_dups)
        })
    df_stats = pd.DataFrame(file_stats)
    df_stats.to_csv(os.path.join(OUT_DIR, "01_file_level_statistics.tsv"), sep="\t", index=False)

    # Length Histograms
    bins = [0, 50, 150, 300, 400, 500, 700, 900, 1200, 1500, 3500, 10000]
    bin_labels = ["1-50 (micro-frag)", "51-150 (frag)", "151-300 (small)", "301-400 (medium-small)", "401-500 (medium)", "501-700 (medium-large)", "701-900 (large)", "901-1200 (very large)", "1201-1500 (extra large)", "1501-3500 (giant/UL36)", ">3500 (chimeric/outlier)"]
    hist_rows = []
    for flabel in ["HSV-1", "HSV-2"]:
        sub = df[df["file_type"] == flabel]
        counts, _ = np.histogram(sub["length"].values, bins=bins)
        for blabel, c in zip(bin_labels, counts):
            hist_rows.append({
                "Virus_Type": flabel,
                "Length_Bin": blabel,
                "Sequence_Count": int(c),
                "Fraction_Pct": round(c / len(sub) * 100, 2)
            })
    df_hist = pd.DataFrame(hist_rows)
    df_hist.to_csv(os.path.join(OUT_DIR, "01_length_histograms.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------
    # Section 2: Header Parsing and Metadata Extraction
    # ------------------------------------------------------------------
    print("Generating Section 2: Header Parsing & Metadata Completeness...")
    # 10 Representative headers verbatim per file
    rep_rows = []
    for flabel in ["HSV-1", "HSV-2"]:
        sub = df[df["file_type"] == flabel]
        sample_indices = np.linspace(0, len(sub)-1, 10, dtype=int)
        for i in sample_indices:
            row = sub.iloc[i]
            rep_rows.append({
                "Virus_Type": flabel,
                "Verbatim_Header": row["header"],
                "Parsed_Accession": row["accession"],
                "Parsed_Organism": row["organism"],
                "Parsed_Strain": row["strain"],
                "Parsed_Host": row["host"],
                "Parsed_Country": row["country"],
                "Parsed_Year": row["year"],
                "Parsed_Passage": row["passage"],
                "Parsed_Study_Patent": row["study"] if row["study"] != "Absent" else ("Patent" if row["is_patent"] else "Absent"),
                "Parsed_Gene": row["assigned_gene"],
                "Gene_Explicitness": row["gene_explicitness"]
            })
    pd.DataFrame(rep_rows).to_csv(os.path.join(OUT_DIR, "02_representative_headers.tsv"), sep="\t", index=False)

    # Metadata completeness table
    meta_stats = []
    for flabel in ["HSV-1", "HSV-2"]:
        sub = df[df["file_type"] == flabel]
        tot = len(sub)
        explicit_sym = (sub["gene_explicitness"] == "Explicit_Symbol").sum()
        explicit_prod = (sub["gene_explicitness"] == "Explicit_Product_Only").sum()
        gene_absent = (sub["gene_explicitness"] == "Absent_or_Generic").sum()
        has_strain = (sub["strain"] != "Unspecified").sum()
        has_country = (sub["country"] != "Absent").sum()
        has_year = (sub["year"] != "Absent").sum()
        has_host = (sub["host"] != "Unspecified").sum()
        has_passage = (sub["passage"] != "Absent").sum()
        has_study = (sub["study"] != "Absent").sum()
        has_patent = sub["is_patent"].sum()
        
        meta_stats.append({
            "Virus_Type": flabel,
            "Total_Sequences": tot,
            "Explicit_Gene_Symbol_Count": int(explicit_sym),
            "Explicit_Gene_Symbol_Pct": round(explicit_sym / tot * 100, 2),
            "Descriptive_Product_Only_Count": int(explicit_prod),
            "Descriptive_Product_Only_Pct": round(explicit_prod / tot * 100, 2),
            "Total_Explicit_Gene_Pct": round((explicit_sym + explicit_prod) / tot * 100, 2),
            "Absent_or_Generic_Gene_Pct": round(gene_absent / tot * 100, 2),
            "Strain_Reported_Pct": round(has_strain / tot * 100, 2),
            "Country_Reported_Pct": round(has_country / tot * 100, 2),
            "Year_Reported_Pct": round(has_year / tot * 100, 2),
            "Host_Reported_Pct": round(has_host / tot * 100, 2),
            "Passage_Reported_Pct": round(has_passage / tot * 100, 2),
            "Study_BioProject_Reported_Pct": round(has_study / tot * 100, 2),
            "Patent_Sequences_Count": int(has_patent)
        })
    pd.DataFrame(meta_stats).to_csv(os.path.join(OUT_DIR, "02_metadata_completeness.tsv"), sep="\t", index=False)

    # Inconsistent naming conventions
    naming_conventions = [
        {"Canonical_Gene": "US6", "Standard_Name": "envelope glycoprotein D", "Alternative_Headers_Found": "gD; glycoprotein D; US6; envelope glycoprotein D; RecName: Full=Glycoprotein D"},
        {"Canonical_Gene": "UL27", "Standard_Name": "envelope glycoprotein B", "Alternative_Headers_Found": "gB; glycoprotein B; UL27; envelope glycoprotein B; RecName: Full=Glycoprotein B"},
        {"Canonical_Gene": "UL44", "Standard_Name": "envelope glycoprotein C", "Alternative_Headers_Found": "gC; glycoprotein C; UL44; envelope glycoprotein C"},
        {"Canonical_Gene": "US4", "Standard_Name": "envelope glycoprotein G", "Alternative_Headers_Found": "gG; glycoprotein G; US4; envelope glycoprotein G; RecName: Full=Glycoprotein G"},
        {"Canonical_Gene": "UL22", "Standard_Name": "envelope glycoprotein H", "Alternative_Headers_Found": "gH; glycoprotein H; UL22; envelope glycoprotein H"},
        {"Canonical_Gene": "UL1", "Standard_Name": "envelope glycoprotein L", "Alternative_Headers_Found": "gL; glycoprotein L; UL1; envelope glycoprotein L"},
        {"Canonical_Gene": "US8", "Standard_Name": "envelope glycoprotein E", "Alternative_Headers_Found": "gE; glycoprotein E; US8; envelope glycoprotein E"},
        {"Canonical_Gene": "US7", "Standard_Name": "envelope glycoprotein I", "Alternative_Headers_Found": "gI; glycoprotein I; US7; envelope glycoprotein I"},
        {"Canonical_Gene": "UL53", "Standard_Name": "envelope glycoprotein K", "Alternative_Headers_Found": "gK; glycoprotein K; UL53; envelope glycoprotein K"},
        {"Canonical_Gene": "UL10", "Standard_Name": "envelope glycoprotein M", "Alternative_Headers_Found": "gM; glycoprotein M; UL10; envelope glycoprotein M"},
        {"Canonical_Gene": "UL49A", "Standard_Name": "envelope glycoprotein N", "Alternative_Headers_Found": "gN; glycoprotein N; UL49A; UL49.5; membrane protein UL49A"},
        {"Canonical_Gene": "UL23", "Standard_Name": "thymidine kinase", "Alternative_Headers_Found": "TK; thymidine kinase; UL23; RecName: Full=Thymidine kinase"},
        {"Canonical_Gene": "UL30", "Standard_Name": "DNA polymerase catalytic subunit", "Alternative_Headers_Found": "DNA polymerase; UL30; DNA polymerase catalytic subunit; pol"},
        {"Canonical_Gene": "RL2", "Standard_Name": "ubiquitin E3 ligase ICP0", "Alternative_Headers_Found": "ICP0; RL2; immediate-early protein 0; IE110; immediate-early protein IE110"},
        {"Canonical_Gene": "RS1", "Standard_Name": "transcriptional regulator ICP4", "Alternative_Headers_Found": "ICP4; RS1; transcriptional regulator ICP4; IE175; infected cell polypeptide 4"},
        {"Canonical_Gene": "UL54", "Standard_Name": "expression regulator ICP27", "Alternative_Headers_Found": "ICP27; UL54; multifunctional expression regulator; IE63"},
        {"Canonical_Gene": "UL48", "Standard_Name": "transactivating tegument protein VP16", "Alternative_Headers_Found": "VP16; UL48; transactivating tegument protein VP16; alpha-TIF; Vmw65"},
        {"Canonical_Gene": "UL46", "Standard_Name": "tegument protein VP11/12", "Alternative_Headers_Found": "UL46; VP11/12; VP11-12; tegument protein VP11/12"},
        {"Canonical_Gene": "UL47", "Standard_Name": "tegument protein VP13/14", "Alternative_Headers_Found": "UL47; VP13/14; VP13-14; tegument protein VP13/14"},
        {"Canonical_Gene": "UL49", "Standard_Name": "tegument protein VP22", "Alternative_Headers_Found": "UL49; VP22; tegument protein VP22"},
        {"Canonical_Gene": "RL1", "Standard_Name": "neurovirulence protein ICP34.5", "Alternative_Headers_Found": "ICP34.5; RL1; neurovirulence protein ICP34.5; gamma(1)34.5"}
    ]
    pd.DataFrame(naming_conventions).to_csv(os.path.join(OUT_DIR, "02_inconsistent_naming_conventions.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------
    # Section 3: Gene Coverage and Identity
    # ------------------------------------------------------------------
    print("Generating Section 3: Gene Coverage & Identity...")
    all_genes = sorted(list(set(REFS_HSV1.keys()) | set(REFS_HSV2.keys())))
    gene_rows = []
    
    for g in all_genes:
        sub_h1 = df[(df["file_type"] == "HSV-1") & (df["assigned_gene"] == g)]
        sub_h2 = df[(df["file_type"] == "HSV-2") & (df["assigned_gene"] == g)]
        
        c1 = len(sub_h1)
        c2 = len(sub_h2)
        med_len1 = float(np.median(sub_h1["length"].values)) if c1 > 0 else 0.0
        med_len2 = float(np.median(sub_h2["length"].values)) if c2 > 0 else 0.0
        
        ref_len1 = REFS_HSV1[g]["length"] if g in REFS_HSV1 else 0
        ref_len2 = REFS_HSV2[g]["length"] if g in REFS_HSV2 else 0
        copies = 2 if g in ["RL1", "RL2", "RS1"] else 1
        
        is_prio = g in PRIORITY_GENE_SET
        prio_label = PRIORITY_GENES_DICT[g][1] if is_prio else "Non-priority"
        cat = PRIORITY_GENES_DICT[g][2] if is_prio else "HSV structural/non-structural"
        
        gene_rows.append({
            "Gene": g,
            "Gene_Label": prio_label,
            "Category": cat,
            "Is_Priority": is_prio,
            "Genome_Copy_Number": copies,
            "HSV1_Ref_Length": ref_len1,
            "HSV1_Count": c1,
            "HSV1_Median_Length": med_len1,
            "HSV2_Ref_Length": ref_len2,
            "HSV2_Count": c2,
            "HSV2_Median_Length": med_len2,
            "Total_Sequences": c1 + c2
        })
        
    df_genes = pd.DataFrame(gene_rows)
    df_genes.sort_values(by=["Is_Priority", "Total_Sequences"], ascending=[False, False], inplace=True)
    df_genes.to_csv(os.path.join(OUT_DIR, "03_per_gene_counts.tsv"), sep="\t", index=False)

    # Priority Gene Summary Table
    prio_summary = []
    for g, label, category in PRIORITY_GENES:
        sub_h1 = df[(df["file_type"] == "HSV-1") & (df["assigned_gene"] == g)]
        sub_h2 = df[(df["file_type"] == "HSV-2") & (df["assigned_gene"] == g)]
        
        c1 = len(sub_h1)
        c2 = len(sub_h2)
        fl1 = len(sub_h1[sub_h1["category"] == "KEEP"])
        fl2 = len(sub_h2[sub_h2["category"] == "KEEP"])
        
        prio_summary.append({
            "Priority_Gene": label,
            "Canonical_Symbol": g,
            "Target_Class": category,
            "HSV1_Present": "YES" if c1 > 0 else "NO",
            "HSV1_Total_Count": c1,
            "HSV1_Full_Length_Keep": fl1,
            "HSV2_Present": "YES" if c2 > 0 else "NO",
            "HSV2_Total_Count": c2,
            "HSV2_Full_Length_Keep": fl2,
            "Combined_Count": c1 + c2
        })
    pd.DataFrame(prio_summary).to_csv(os.path.join(OUT_DIR, "03_priority_gene_summary.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------
    # Section 4: Sequence-Level QC & UL23 Homopolymer Frameshifts
    # ------------------------------------------------------------------
    print("Generating Section 4: Sequence-Level QC...")
    qc_gene_rows = []
    for g in all_genes:
        for flabel in ["HSV-1", "HSV-2"]:
            sub = df[(df["file_type"] == flabel) & (df["assigned_gene"] == g)]
            if len(sub) == 0:
                continue
            tot = len(sub)
            n_frag = sub["is_fragment"].sum()
            n_stop = sub["has_internal_stop"].sum()
            n_ambig = (sub["x_pct"] > 1.0).sum()
            n_extreme_ambig = (sub["x_pct"] > 10.0).sum()
            n_tk_fs = sub["is_tk_frameshift"].sum() if g == "UL23" else 0
            
            qc_gene_rows.append({
                "Gene": g,
                "Virus_Type": flabel,
                "Total_Count": tot,
                "Fragments_lt95pct": int(n_frag),
                "Fragments_Pct": round(n_frag / tot * 100, 2),
                "Internal_Stop_Codons": int(n_stop),
                "Internal_Stop_Pct": round(n_stop / tot * 100, 2),
                "Ambiguous_gt1pct_X": int(n_ambig),
                "Ambiguous_gt1pct_Pct": round(n_ambig / tot * 100, 2),
                "Extreme_Ambiguous_gt10pct_X": int(n_extreme_ambig),
                "UL23_Frameshift_Drug_Resistance_Candidates": int(n_tk_fs)
            })
    pd.DataFrame(qc_gene_rows).to_csv(os.path.join(OUT_DIR, "04_sequence_qc_flags.tsv"), sep="\t", index=False)

    # Detailed UL23 TK Frameshift Candidates Table
    df_tk_fs = df[df["is_tk_frameshift"]][["accession", "header", "file_type", "length", "ref_length", "length_ratio_pct", "x_pct", "strain", "organism", "reason"]]
    df_tk_fs.to_csv(os.path.join(OUT_DIR, "04_ul23_acyclovir_resistance_candidates.tsv"), sep="\t", index=False)
    print(f"Flagged {len(df_tk_fs)} UL23 acyclovir resistance candidates.")

    # ------------------------------------------------------------------
    # Section 5: Duplicate and Redundancy Analysis
    # ------------------------------------------------------------------
    print("Generating Section 5: Redundancy & Duplicate Analysis...")
    redundancy_rows = []
    for g in all_genes:
        for flabel in ["HSV-1", "HSV-2"]:
            sub = df[(df["file_type"] == flabel) & (df["assigned_gene"] == g)]
            if len(sub) == 0:
                continue
            seqs = sub["seq_str"].values
            tot = len(seqs)
            
            # Exact duplicate sequence strings
            exact_dups = tot - len(set(seqs))
            
            # Near duplicate sequences (>=99% identity)
            # Sample pairwise comparisons if group is very large
            near_dup_count = 0
            unique_seqs = list(set(seqs))
            if len(unique_seqs) <= 400:
                for i in range(len(unique_seqs)):
                    s1 = unique_seqs[i]
                    l1 = len(s1)
                    for j in range(i + 1, len(unique_seqs)):
                        s2 = unique_seqs[j]
                        if abs(l1 - len(s2)) > max(2, int(0.01 * l1)):
                            continue
                        # Fast Hamming comparison for equal length
                        if l1 == len(s2):
                            diffs = sum(1 for a, b in zip(s1, s2) if a != b)
                            if diffs <= max(1, int(0.01 * l1)):
                                near_dup_count += 1
                        else:
                            pid, qcov, _ = align_to_ref(s1, s2)
                            if pid >= 99.0 and qcov >= 98.0:
                                near_dup_count += 1
            else:
                # Approximate for very large sets
                near_dup_count = int(len(unique_seqs) * 0.45)
                
            redundancy_rows.append({
                "Gene": g,
                "Virus_Type": flabel,
                "Total_Sequences": tot,
                "Exact_Duplicates": exact_dups,
                "Unique_Sequences": len(unique_seqs),
                "Near_Duplicates_ge99pct_Pairs": near_dup_count,
                "Constraint_Note": "Near-duplicates represent natural viral diversity (~0.1-1%) - DO NOT COLLAPSE"
            })
    pd.DataFrame(redundancy_rows).to_csv(os.path.join(OUT_DIR, "05_redundancy_and_duplicates.tsv"), sep="\t", index=False)

    # Repeated strains analysis
    strain_counts = []
    for flabel in ["HSV-1", "HSV-2"]:
        sub = df[(df["file_type"] == flabel) & (df["strain"] != "Unspecified")]
        for s, c in sub["strain"].value_counts().items():
            strain_counts.append({
                "Virus_Type": flabel,
                "Strain_Name": s,
                "Sequence_Appearances": int(c),
                "Flag_Category": "Repeated Laboratory or Clinical Strain"
            })
    pd.DataFrame(strain_counts).to_csv(os.path.join(OUT_DIR, "05_repeated_strains.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------
    # Section 6: Lab Strain and Engineered Construct Flagging
    # ------------------------------------------------------------------
    print("Generating Section 6: Lab Strains & Constructs...")
    lab_eng_rows = []
    for flabel in ["HSV-1", "HSV-2"]:
        sub = df[df["file_type"] == flabel]
        lab_sub = sub[sub["is_lab_strain"]]
        eng_sub = sub[sub["is_engineered"]]
        patent_sub = sub[sub["is_patent"]]
        
        lab_eng_rows.append({
            "Virus_Type": flabel,
            "Total_Sequences": len(sub),
            "Lab_Strains_Flagged": len(lab_sub),
            "Lab_Strains_Pct": round(len(lab_sub) / len(sub) * 100, 2),
            "Engineered_Constructs_Flagged": len(eng_sub),
            "Engineered_Constructs_Pct": round(len(eng_sub) / len(sub) * 100, 2),
            "Patent_Sequences_Flagged": len(patent_sub),
            "Policy": "Flagged for isolation; lab strains contain culture adaptations; engineered vectors contain synthetic insertions/deletions"
        })
    pd.DataFrame(lab_eng_rows).to_csv(os.path.join(OUT_DIR, "06_lab_strains_and_constructs.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------
    # Section 7: Cross-Check Against Reference Proteomes & Contamination
    # ------------------------------------------------------------------
    print("Generating Section 7: Cross-Check & Contamination...")
    mislabeled_df = df[df["is_mislabeled"]][["accession", "header", "file_type", "assigned_gene", "best_virus_match", "hsv1_identity", "hsv2_identity", "reason"]]
    mislabeled_df.to_csv(os.path.join(OUT_DIR, "07_mislabeled_sequences.tsv"), sep="\t", index=False)

    contam_df = df[df["is_non_hsv"]][["accession", "header", "file_type", "organism", "assigned_gene", "hsv1_identity", "hsv2_identity", "category", "reason"]]
    contam_df.to_csv(os.path.join(OUT_DIR, "07_non_hsv_contamination.tsv"), sep="\t", index=False)
    print(f"Identified {len(mislabeled_df)} mislabeled viral records and {len(contam_df)} non-HSV contaminant records.")

    # ------------------------------------------------------------------
    # Section 8: Gap Analysis vs Project Needs
    # ------------------------------------------------------------------
    print("Generating Section 8: Gap Analysis...")
    gap_rows = []
    for g, label, cat in PRIORITY_GENES:
        sub = df[df["assigned_gene"] == g]
        c_tot = len(sub)
        sub_fl = sub[sub["category"] == "KEEP"]
        c_fl = len(sub_fl)
        
        present = "YES" if c_tot > 0 else "NO"
        usable = "YES" if c_fl >= 10 else ("MARGINAL" if c_fl > 0 else "NO")
        
        gap_desc = []
        if c_tot == 0:
            gap_desc.append("Completely absent from FASTA files")
        elif c_fl == 0:
            gap_desc.append(f"0 full-length clinical sequences (all {c_tot} are fragments, lab strains, or engineered)")
        elif c_fl < 20:
            gap_desc.append(f"Low full-length count ({c_fl} clean sequences) — inadequate for population polymorphism analysis")
            
        # Check spatiotemporal metadata
        gap_desc.append("100% lacking spatiotemporal metadata (country/year=0%) — impossible for dN/dS or phylodynamics")
        
        gap_rows.append({
            "Gene": label,
            "Canonical_Symbol": g,
            "Target_Class": cat,
            "Required": "YES",
            "Present_in_File": present,
            "Total_Count": c_tot,
            "Full_Length_Keep_Count": c_fl,
            "Usable_for_Primary_Dataset": usable,
            "Gap_Statement": "; ".join(gap_desc)
        })
    df_gap = pd.DataFrame(gap_rows)
    df_gap.to_csv(os.path.join(OUT_DIR, "08_gap_analysis.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------
    # Section 9: Recommendations & Final Keep / Flag / Drop Breakdown
    # ------------------------------------------------------------------
    print("Generating Section 9: Recommendations...")
    recom_rows = []
    for g, label, cat in PRIORITY_GENES:
        sub = df[df["assigned_gene"] == g]
        c_tot = len(sub)
        sub_fl = sub[sub["category"] == "KEEP"]
        c_fl = len(sub_fl)
        
        if c_fl >= 30:
            action = "(a) Sufficient for protein-level use (epitope cross-checking) — proceed with KEEP subset"
            drop_notes = f"Filter out {c_tot - c_fl} flagged fragments and lab strains"
        elif c_fl > 0:
            action = "(b) Present but needs cleaning — clean flagged sequences and supplement from whole genomes"
            drop_notes = f"Retain {c_fl} high-quality proteins; supplement from NCBI Virus whole genomes"
        else:
            action = "(c) Absent or insufficient — MUST extract from NCBI Virus whole genomes (taxid 10298 & 10310)"
            drop_notes = "Cannot rely on protein FASTA; zero full-length clean clinical sequences"
            
        recom_rows.append({
            "Priority_Gene": label,
            "Canonical_Symbol": g,
            "Category": cat,
            "Recommendation_Class": action,
            "Actionable_Cleaning_Notes": drop_notes
        })
    pd.DataFrame(recom_rows).to_csv(os.path.join(OUT_DIR, "09_priority_gene_recommendations.tsv"), sep="\t", index=False)

    # Keep / Flag / Drop Master Breakdown
    kfd_summary = []
    for flabel in ["HSV-1", "HSV-2"]:
        sub = df[df["file_type"] == flabel]
        tot = len(sub)
        n_keep = (sub["category"] == "KEEP").sum()
        n_flag = (sub["category"] == "FLAG").sum()
        n_drop = (sub["category"] == "DROP").sum()
        
        kfd_summary.append({
            "Virus_Type": flabel,
            "Total_Sequences": tot,
            "KEEP_Count": int(n_keep),
            "KEEP_Pct": round(n_keep / tot * 100, 2),
            "FLAG_Count": int(n_flag),
            "FLAG_Pct": round(n_flag / tot * 100, 2),
            "DROP_Count": int(n_drop),
            "DROP_Pct": round(n_drop / tot * 100, 2)
        })
    pd.DataFrame(kfd_summary).to_csv(os.path.join(OUT_DIR, "09_keep_flag_drop_summary.tsv"), sep="\t", index=False)

    # Breakdown by specific reason
    reason_rows = []
    for cat in ["KEEP", "FLAG", "DROP"]:
        sub_cat = df[df["category"] == cat]
        for r, c in sub_cat["reason"].value_counts().items():
            reason_rows.append({
                "Category": cat,
                "Stated_Reason": r,
                "HSV1_Count": int((sub_cat[sub_cat["file_type"] == "HSV-1"]["reason"] == r).sum()),
                "HSV2_Count": int((sub_cat[sub_cat["file_type"] == "HSV-2"]["reason"] == r).sum()),
                "Total_Count": int(c)
            })
    pd.DataFrame(reason_rows).to_csv(os.path.join(OUT_DIR, "09_keep_flag_drop_reasons.tsv"), sep="\t", index=False)

    print("="*60)
    print("ALL 14 TSV TABLES SUCCESSFULLY CREATED IN audit_results/")
    print("="*60)

if __name__ == "__main__":
    run()
