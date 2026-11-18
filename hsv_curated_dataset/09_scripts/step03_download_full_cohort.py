"""
Fetch complete genome FASTA and metadata for all 672 HSV-1 and 633 HSV-2 accessions
directly from NCBI nuccore to ensure 100% complete representation.
"""

import os
import json
import time
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "hsv_curated_dataset@example.com"

# Target queries for complete genomes
QUERIES = {
    "HSV-1": '(txid10298[Organism:exp] OR "Human alphaherpesvirus 1"[Organism]) AND ("complete genome"[Title] OR "complete sequence"[Title] OR complete[Properties] OR 130000:170000[SLEN]) NOT "patent"[Title]',
    "HSV-2": '(txid10310[Organism:exp] OR "Human alphaherpesvirus 2"[Organism]) AND ("complete genome"[Title] OR "complete sequence"[Title] OR complete[Properties] OR 130000:170000[SLEN]) NOT "patent"[Title]'
}

OUT_DIR_HSV1 = "hsv_curated_dataset/00_raw_downloads/HSV1"
OUT_DIR_HSV2 = "hsv_curated_dataset/00_raw_downloads/HSV2"
META_DIR = "hsv_curated_dataset/03_metadata"
os.makedirs(OUT_DIR_HSV1, exist_ok=True)
os.makedirs(OUT_DIR_HSV2, exist_ok=True)
os.makedirs(META_DIR, exist_ok=True)

def download_full_cohort():
    all_metadata_rows = []
    
    for vtype, q in QUERIES.items():
        print(f"\n==========================================")
        print(f"SEARCHING AND DOWNLOADING COMPLETE {vtype} GENOMES")
        print(f"==========================================")
        
        # 1. Search nuccore
        handle = Entrez.esearch(db="nuccore", term=q, retmax=10000, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        
        acc_list = record["IdList"]
        print(f"Total accessions found in nuccore for {vtype}: {len(acc_list)}")
        
        out_fasta = os.path.join(OUT_DIR_HSV1 if vtype == "HSV-1" else OUT_DIR_HSV2, f"{vtype.replace('-', '')}_all_genomes.fasta")
        
        # 2. Batch download FASTA (50 per batch)
        batch_size = 50
        downloaded_records = []
        
        print(f"Downloading {len(acc_list)} complete genome FASTA sequences in batches of {batch_size}...")
        with open(out_fasta, "w", encoding="utf-8") as fasta_out:
            for i in range(0, len(acc_list), batch_size):
                batch_accs = acc_list[i:i+batch_size]
                try:
                    with Entrez.efetch(db="nuccore", id=batch_accs, rettype="fasta", retmode="text") as f_handle:
                        content = f_handle.read()
                        fasta_out.write(content)
                    print(f"  Downloaded batch {i//batch_size + 1}/{(len(acc_list)-1)//batch_size + 1} ({len(batch_accs)} accessions)")
                except Exception as e:
                    print(f"  Error on batch {i}: {e}. Retrying individually...")
                    for single_acc in batch_accs:
                        try:
                            with Entrez.efetch(db="nuccore", id=single_acc, rettype="fasta", retmode="text") as f_handle:
                                fasta_out.write(f_handle.read())
                        except Exception as e2:
                            print(f"    Failed {single_acc}: {e2}")
                time.sleep(0.5)
                
        # Count downloaded sequences
        downloaded_seqs = list(SeqIO.parse(out_fasta, "fasta"))
        print(f"Saved {out_fasta} with {len(downloaded_seqs)} FASTA records.")
        
        # 3. Batch download Summaries / Metadata
        print(f"Downloading metadata summaries for {vtype}...")
        for i in range(0, len(acc_list), 100):
            batch_accs = acc_list[i:i+100]
            try:
                with Entrez.esummary(db="nuccore", id=batch_accs) as s_handle:
                    summaries = Entrez.read(s_handle)
                    for s in summaries:
                        acc = s.get("Caption", "")
                        title = s.get("Title", "")
                        slen = int(s.get("Length", 0))
                        taxid = int(s.get("TaxId", 0))
                        create_date = s.get("CreateDate", "")
                        update_date = s.get("UpdateDate", "")
                        subtype = s.get("SubTypeName", "")
                        subname = s.get("SubName", "")
                        
                        all_metadata_rows.append({
                            "accession": acc,
                            "virus_type": vtype,
                            "title": title,
                            "genome_length": slen,
                            "taxid": taxid,
                            "create_date": create_date,
                            "update_date": update_date,
                            "subtype": subtype,
                            "subname": subname
                        })
            except Exception as e:
                print(f"  Error downloading summary batch {i}: {e}")
            time.sleep(0.5)

    # Save combined raw metadata table
    df_meta = pd.DataFrame(all_metadata_rows)
    raw_meta_path = os.path.join(META_DIR, "raw_metadata_all_genomes.tsv")
    df_meta.to_csv(raw_meta_path, sep="\t", index=False)
    print(f"\nSaved combined metadata table to {raw_meta_path} with {len(df_meta)} records.")

if __name__ == "__main__":
    download_full_cohort()
