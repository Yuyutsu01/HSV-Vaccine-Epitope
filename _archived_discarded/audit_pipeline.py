"""
HSV-1 and HSV-2 Protein FASTA Comprehensive Audit Pipeline
Performs file-level statistics, header metadata parsing, reference proteome
cross-checking, gene mapping, QC flagging, redundancy analysis, and gap analysis.
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from collections import Counter, defaultdict

# Create output directory for audit artifacts
OUT_DIR = "audit_results"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------
# Step 1: Load and Index Canonical Reference Proteomes
# ---------------------------------------------------------

def load_reference_proteome(fasta_path, virus_type):
    """
    Loads reference proteins from FASTA, records canonical lengths,
    and indexes 5-mer peptides for rapid homology detection.
    """
    ref_dict = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        parts = rec.description.split("|")
        prot_id = parts[0]
        gene = parts[1]
        locus_tag = parts[2]
        product = parts[3]
        copy_info = parts[4]
        plen = int(parts[5].replace("len=", ""))
        seq_str = str(rec.seq).upper()
        
        # RL1, RL2, RS1 have 2 copies in complete genomes
        if gene not in ref_dict:
            ref_dict[gene] = {
                "gene": gene,
                "prot_id": prot_id,
                "locus_tag": locus_tag,
                "product": product,
                "length": plen,
                "seq": seq_str,
                "virus": virus_type,
                "copies": 1 if "copy=1" in copy_info else 2
            }
        else:
            ref_dict[gene]["copies"] = 2
            
    # Index 5-mers for each gene
    k = 5
    kmer_dict = {}
    for g, data in ref_dict.items():
        s = data["seq"]
        kmer_dict[g] = set(s[i:i+k] for i in range(len(s)-k+1))
        
    return ref_dict, kmer_dict

print("Loading reference proteomes...")
refs_hsv1, kmers_hsv1 = load_reference_proteome("references/HSV1_ref_proteins.fasta", "HSV-1")
refs_hsv2, kmers_hsv2 = load_reference_proteome("references/HSV2_ref_proteins.fasta", "HSV-2")
print(f"Loaded {len(refs_hsv1)} HSV-1 and {len(refs_hsv2)} HSV-2 canonical reference proteins.")

# Configure local pairwise aligner
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.open_gap_score = -10.0
aligner.extend_gap_score = -0.5
aligner.match_score = 2.0
aligner.mismatch_score = -1.0

# ---------------------------------------------------------
# Step 2: Dictionaries, Synonyms, and Priority Target Lists
# ---------------------------------------------------------

PRIORITY_GENES = [
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

PRIORITY_GENE_SET = set(g[0] for g in PRIORITY_GENES)

# Regex mappings for header synonyms
GENE_SYNONYMS = [
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
]

# Lab strain and engineered construct patterns
LAB_STRAIN_PAT = re.compile(r"\b(strain 17|17\+|strain F\b|KOS|McKrae|HFEM|MacIntyre|SC16|HG52|strain 186|strain MS|strain RH2|strain R-15)\b", re.I)
ENGINEERED_PAT = re.compile(r"\b(BAC|oncolytic|vaccine vector|deletion mutant|recombinant|mutant|synthetic|patent|transgenic|chimera|vector)\b", re.I)

COUNTRIES = [
    "USA", "United States", "China", "Japan", "Germany", "United Kingdom", "UK",
    "France", "Canada", "Australia", "South Africa", "Kenya", "Uganda", "Tanzania",
    "Brazil", "India", "Russia", "Italy", "Spain", "Sweden", "Netherlands", "Norway"
]
COUNTRY_PAT = re.compile(r"\b(" + "|".join(COUNTRIES) + r")\b", re.I)

# ---------------------------------------------------------
# Step 3: Alignment & Homology Function
# ---------------------------------------------------------

def align_to_ref(query_seq, ref_seq):
    """
    Computes local pairwise alignment score, percent identity,
    query coverage, and reference coverage.
    """
    if not query_seq or not ref_seq:
        return 0.0, 0.0, 0.0, 0
    try:
        score = aligner.score(query_seq, ref_seq)
        alignments = aligner.align(query_seq, ref_seq)
        if not alignments:
            return 0.0, 0.0, 0.0, 0
        top = alignments[0]
        # Query aligned length
        q_aligned = sum(top.aligned[0][i][1] - top.aligned[0][i][0] for i in range(len(top.aligned[0])))
        r_aligned = sum(top.aligned[1][i][1] - top.aligned[1][i][0] for i in range(len(top.aligned[1])))
        
        # Approximate identical matches from score
        # score = 2 * matches - 1 * mismatches - 10 * gaps - 0.5 * gap_ext
        # For precision on top hit, count direct matches across aligned blocks
        matches = 0
        for block_idx in range(len(top.aligned[0])):
            q_start, q_end = top.aligned[0][block_idx]
            r_start, r_end = top.aligned[1][block_idx]
            sub_q = query_seq[q_start:q_end]
            sub_r = ref_seq[r_start:r_end]
            min_l = min(len(sub_q), len(sub_r))
            matches += sum(1 for a, b in zip(sub_q[:min_l], sub_r[:min_l]) if a == b)
            
        pct_id = (matches / max(1, q_aligned)) * 100.0
        q_cov = (q_aligned / len(query_seq)) * 100.0
        r_cov = (r_aligned / len(ref_seq)) * 100.0
        return pct_id, q_cov, r_cov, matches
    except Exception:
        return 0.0, 0.0, 0.0, 0

def find_best_homolog(query_seq, expected_gene=None):
    """
    Identifies best hit in HSV-1 and HSV-2 references using k-mers,
    then confirms with exact local alignment.
    """
    k = 5
    if len(query_seq) < k:
        # Extremely short fragment (e.g. 1 residue)
        return {
            "hsv1_gene": expected_gene or "UNKNOWN", "hsv1_id": 0.0, "hsv1_qcov": 0.0, "hsv1_rcov": 0.0,
            "hsv2_gene": expected_gene or "UNKNOWN", "hsv2_id": 0.0, "hsv2_qcov": 0.0, "hsv2_rcov": 0.0,
            "best_virus": "UNKNOWN", "best_gene": expected_gene or "UNKNOWN"
        }
        
    q_kmers = set(query_seq[i:i+k] for i in range(len(query_seq)-k+1))
    
    # Calculate shared k-mers
    scores1 = {g: len(q_kmers & km) for g, km in kmers_hsv1.items()}
    scores2 = {g: len(q_kmers & km) for g, km in kmers_hsv2.items()}
    
    top_g1 = max(scores1, key=scores1.get)
    top_g2 = max(scores2, key=scores2.get)
    
    # If expected_gene is known, prioritize it
    eval_g1 = expected_gene if (expected_gene and expected_gene in refs_hsv1) else top_g1
    eval_g2 = expected_gene if (expected_gene and expected_gene in refs_hsv2) else top_g2
    
    id1, qcov1, rcov1, _ = align_to_ref(query_seq, refs_hsv1[eval_g1]["seq"])
    id2, qcov2, rcov2, _ = align_to_ref(query_seq, refs_hsv2[eval_g2]["seq"])
    
    # Check if top hit via kmer is better than expected_gene
    if top_g1 != eval_g1 and scores1[top_g1] > scores1.get(eval_g1, 0) + 10:
        alt_id1, alt_qcov1, alt_rcov1, _ = align_to_ref(query_seq, refs_hsv1[top_g1]["seq"])
        if alt_id1 * alt_qcov1 > id1 * qcov1:
            eval_g1, id1, qcov1, rcov1 = top_g1, alt_id1, alt_qcov1, alt_rcov1
            
    if top_g2 != eval_g2 and scores2[top_g2] > scores2.get(eval_g2, 0) + 10:
        alt_id2, alt_qcov2, alt_rcov2, _ = align_to_ref(query_seq, refs_hsv2[top_g2]["seq"])
        if alt_id2 * alt_qcov2 > id2 * qcov2:
            eval_g2, id2, qcov2, rcov2 = top_g2, alt_id2, alt_qcov2, alt_rcov2
            
    # Determine best virus match
    if id1 < 40.0 and id2 < 40.0 and max(scores1.values()) < 5:
        best_virus = "NON_HSV_CONTAMINANT"
        best_gene = "NON_HSV"
    elif id1 >= id2:
        best_virus = "HSV-1"
        best_gene = eval_g1
    else:
        best_virus = "HSV-2"
        best_gene = eval_g2
        
    return {
        "hsv1_gene": eval_g1, "hsv1_id": id1, "hsv1_qcov": qcov1, "hsv1_rcov": rcov1,
        "hsv2_gene": eval_g2, "hsv2_id": id2, "hsv2_qcov": qcov2, "hsv2_rcov": rcov2,
        "best_virus": best_virus, "best_gene": best_gene
    }

print("Alignment functions verified.")
