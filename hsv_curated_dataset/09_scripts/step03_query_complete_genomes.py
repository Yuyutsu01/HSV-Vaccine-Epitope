"""
Step 3 Query: Search NCBI nuccore for all complete HSV-1 and HSV-2 genomes
and report accession counts and metadata summaries.
"""

import time
import pandas as pd
from Bio import Entrez

Entrez.email = "hsv_curated_dataset@example.com"

# Target queries for complete genomes
QUERIES = {
    "HSV-1": '(txid10298[Organism:exp] OR "Human alphaherpesvirus 1"[Organism]) AND ("complete genome"[Title] OR "complete sequence"[Title] OR complete[Properties] OR 130000:170000[SLEN]) NOT "patent"[Title]',
    "HSV-2": '(txid10310[Organism:exp] OR "Human alphaherpesvirus 2"[Organism]) AND ("complete genome"[Title] OR "complete sequence"[Title] OR complete[Properties] OR 130000:170000[SLEN]) NOT "patent"[Title]'
}

def query_ncbi_virus():
    counts = {}
    id_lists = {}
    
    for vtype, q in QUERIES.items():
        print(f"Querying NCBI for {vtype}...")
        # Search nuccore
        handle = Entrez.esearch(db="nuccore", term=q, retmax=10000, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        
        acc_list = record["IdList"]
        count = int(record["Count"])
        counts[vtype] = count
        id_lists[vtype] = acc_list
        print(f"  {vtype}: Found {count} records matching query (returned {len(acc_list)} accessions).")
        time.sleep(1)
        
    return counts, id_lists

if __name__ == "__main__":
    counts, id_lists = query_ncbi_virus()
    print("\n=== SUMMARY OF COMPLETE GENOME COUNTS ===")
    for vtype, c in counts.items():
        print(f"{vtype}: {c} complete/near-complete genomes found")
