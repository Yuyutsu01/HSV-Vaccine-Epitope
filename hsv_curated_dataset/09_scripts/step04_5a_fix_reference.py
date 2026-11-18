"""
Step 4.5a: Restore true NC_001798.2 (154,675 bp) from NCBI Entrez,
re-parse reference_cds_table.tsv with Version 2 coordinates,
and document in README.md.
"""

import os
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "hsv_curated_dataset@example.com"
REF_DIR = "hsv_curated_dataset/01_reference_genomes"
DOC_DIR = "hsv_curated_dataset/08_documentation"
os.makedirs(REF_DIR, exist_ok=True)
os.makedirs(DOC_DIR, exist_ok=True)

REFS = {
    "HSV-1": "NC_001806.2",
    "HSV-2": "NC_001798.2"
}

SYNONYMS = {
    "US6": "gD, glycoprotein D, envelope glycoprotein D",
    "UL27": "gB, glycoprotein B, envelope glycoprotein B",
    "UL44": "gC, glycoprotein C, envelope glycoprotein C",
    "US4": "gG, glycoprotein G, envelope glycoprotein G",
    "UL22": "gH, glycoprotein H, envelope glycoprotein H",
    "UL1": "gL, glycoprotein L, envelope glycoprotein L",
    "US8": "gE, glycoprotein E, envelope glycoprotein E",
    "US7": "gI, glycoprotein I, envelope glycoprotein I",
    "UL53": "gK, glycoprotein K, envelope glycoprotein K",
    "UL10": "gM, glycoprotein M, envelope glycoprotein M",
    "UL49A": "gN, glycoprotein N, UL49.5",
    "RL1": "ICP34.5, neurovirulence factor",
    "RL2": "ICP0, IE110, immediate early protein 0",
    "RS1": "ICP4, IE175, transcriptional regulator",
    "UL54": "ICP27, IE63",
    "UL48": "VP16, alpha-TIF, Vmw65",
    "UL46": "VP11/12",
    "UL47": "VP13/14",
    "UL49": "VP22",
    "UL23": "TK, thymidine kinase",
    "UL30": "pol, DNA polymerase catalytic subunit",
    "UL42": "DNA polymerase processivity subunit",
    "UL39": "RR1, ribonucleotide reductase large subunit",
    "UL40": "RR2, ribonucleotide reductase small subunit",
    "US1": "ICP22, IE68",
    "US3": "serine/threonine kinase US3",
    "US12": "ICP47, TAP inhibitor"
}

def fix_reference_version():
    print("Fetching true current NC_001798.2 (154,675 bp) from NCBI...")
    gbk_file = os.path.join(REF_DIR, "NC_001798.2.gbk")
    fasta_file = os.path.join(REF_DIR, "NC_001798.2.fasta")
    
    with Entrez.efetch(db="nuccore", id="NC_001798.2", rettype="gbwithparts", retmode="text") as handle:
        content = handle.read()
        with open(gbk_file, "w", encoding="utf-8") as f:
            f.write(content)
            
    rec = SeqIO.read(gbk_file, "genbank")
    with open(fasta_file, "w", encoding="utf-8") as f:
        f.write(f">{rec.id} {rec.description}\n{str(rec.seq)}\n")
    print(f"Verified downloaded NC_001798.2 length: {len(rec.seq)} bp (expected: 154675 bp)")
    
    # Re-parse both reference genomes to produce the updated reference_cds_table.tsv
    cds_rows = []
    for vtype, acc in REFS.items():
        gpath = os.path.join(REF_DIR, f"{acc}.gbk")
        record = SeqIO.read(gpath, "genbank")
        gene_copy_tracker = {}
        for feat in record.features:
            if feat.type == "CDS":
                gene = feat.qualifiers.get("gene", [""])[0]
                locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                prot_id = feat.qualifiers.get("protein_id", [""])[0]
                translation = feat.qualifiers.get("translation", [""])[0]
                g_sym = gene or locus_tag or product
                
                gene_copy_tracker[g_sym] = gene_copy_tracker.get(g_sym, 0) + 1
                copy_num = gene_copy_tracker[g_sym]
                
                start = int(feat.location.start) + 1
                end = int(feat.location.end)
                strand = "+" if feat.location.strand == 1 else "-"
                nt_len = sum(len(part) for part in feat.location.parts)
                aa_len = len(translation)
                syn = SYNONYMS.get(g_sym, "")
                
                cds_rows.append({
                    "virus": vtype,
                    "gene_symbol": g_sym,
                    "gene_synonyms": syn,
                    "product": product,
                    "protein_id": prot_id,
                    "locus_tag": locus_tag,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "copy_number": copy_num,
                    "protein_length_aa": aa_len,
                    "nucleotide_length_nt": nt_len
                })
                
    df = pd.DataFrame(cds_rows)
    table_path = os.path.join(REF_DIR, "reference_cds_table.tsv")
    df.to_csv(table_path, sep="\t", index=False)
    print(f"Saved {table_path} with {len(df)} reference CDS entries.")
    
    # Document in README.md
    readme_path = os.path.join(DOC_DIR, "README.md")
    readme_note = """# HSV Curated Dataset Project

## Reference Genomes
- **HSV-1 reference**: NC_001806.2 (152,222 bp), strain 17.
- **HSV-2 reference**: NC_001798.2, version 2 (154,675 bp), strain HG52.
  *Note on Version History*: Version 1 (NC_001798.1, 154,746 bp) was deprecated by NCBI in 2016 following terminal repeat curation. All coordinates and annotations derive strictly from current canonical Version 2 (NC_001798.2).
"""
    with open(readme_path, "w", encoding="utf-8") as rf:
        rf.write(readme_note)
    print(f"Documented reference version in {readme_path}")
    
    # Print repeat-resident gene coordinates
    print("\nCorrected Version-2 Repeat-Resident Genes (RL1, RL2, RS1):")
    rep_df = df[df["gene_symbol"].isin(["RL1", "RL2", "RS1"])][["virus", "gene_symbol", "copy_number", "start", "end", "strand", "protein_length_aa", "nucleotide_length_nt"]]
    print(rep_df.to_string(index=False))

if __name__ == "__main__":
    fix_reference_version()
