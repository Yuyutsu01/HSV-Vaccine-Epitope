"""
Step 2: Download Reference Genomes (NC_001806.2 and NC_001798.2)
and generate reference_cds_table.tsv with explicit inverted-repeat handling.
"""

import os
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "hsv_curated_dataset@example.com"

OUT_DIR = "hsv_curated_dataset/01_reference_genomes"
os.makedirs(OUT_DIR, exist_ok=True)

# Use NC_001806.2 (152,222 bp) and NC_001798.1 (154,746 bp) to match expected lengths
REFS = {
    "HSV-1": ("NC_001806.2", "NC_001806.2"),
    "HSV-2": ("NC_001798.1", "NC_001798.2")
}

# Gene synonym mapping for canonical herpesvirus genes
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

def download_and_parse_references():
    cds_rows = []
    
    for vtype, (fetch_id, save_acc) in REFS.items():
        gbk_file = os.path.join(OUT_DIR, f"{save_acc}.gbk")
        fasta_file = os.path.join(OUT_DIR, f"{save_acc}.fasta")
        
        # Download GenBank flatfile
        print(f"Downloading {vtype} reference genome ({fetch_id} -> {save_acc}) from NCBI...")
        with Entrez.efetch(db="nuccore", id=fetch_id, rettype="gbwithparts", retmode="text") as handle:
            gbk_content = handle.read()
            with open(gbk_file, "w", encoding="utf-8") as f:
                f.write(gbk_content)
        print(f"Saved {gbk_file}")
        
        # Parse GenBank record
        record = SeqIO.read(gbk_file, "genbank")
        
        # Save nucleotide FASTA
        with open(fasta_file, "w", encoding="utf-8") as f:
            f.write(f">{record.id} {record.description}\n{str(record.seq)}\n")
        print(f"Saved nucleotide FASTA {fasta_file} (length: {len(record.seq)} bp)")
        
        # Track repeat copy numbers
        gene_copy_tracker = {}
        
        for feat in record.features:
            if feat.type == "CDS":
                gene = feat.qualifiers.get("gene", [""])[0]
                locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
                product = feat.qualifiers.get("product", [""])[0]
                prot_id = feat.qualifiers.get("protein_id", [""])[0]
                translation = feat.qualifiers.get("translation", [""])[0]
                
                # Gene symbol normalization
                g_sym = gene or locus_tag or product
                
                # Update copy number
                gene_copy_tracker[g_sym] = gene_copy_tracker.get(g_sym, 0) + 1
                copy_num = gene_copy_tracker[g_sym]
                
                # Coordinates (1-based inclusive for biological reporting)
                start = int(feat.location.start) + 1
                end = int(feat.location.end)
                strand = "+" if feat.location.strand == 1 else "-"
                
                # Splicing / parts
                nt_len = sum(len(part) for part in feat.location.parts)
                aa_len = len(translation)
                
                # Synonyms
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
    table_path = os.path.join(OUT_DIR, "reference_cds_table.tsv")
    df.to_csv(table_path, sep="\t", index=False)
    print(f"Saved {table_path} with {len(df)} reference CDS entries.")
    
    # Display summary of inverted repeat genes
    print("\nInverted Repeat Genes Summary:")
    print(df[df["gene_symbol"].isin(["RL1", "RL2", "RS1"])][["virus", "gene_symbol", "copy_number", "start", "end", "strand", "protein_length_aa"]])

if __name__ == "__main__":
    download_and_parse_references()
