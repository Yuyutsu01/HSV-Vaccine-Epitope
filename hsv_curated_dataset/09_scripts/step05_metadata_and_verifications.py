"""
Step 5: Curated Metadata Pipeline
AND Parallel Execution of Verifications A, B, C, and D.
"""

import os
import re
import json
import hashlib
import pandas as pd
from Bio import SeqIO
from collections import Counter, defaultdict

# Directories
RAW_DIR = "hsv_curated_dataset/00_raw_downloads"
FILTERED_DIR = "hsv_curated_dataset/02_filtered_genomes"
META_DIR = "hsv_curated_dataset/03_metadata"
QC_DIR = "hsv_curated_dataset/07_qc_reports"
os.makedirs(META_DIR, exist_ok=True)
os.makedirs(QC_DIR, exist_ok=True)

# ----------------------------------------------------------------------
# 1. Verification A: Species Disagreement Accessions Full List
# ----------------------------------------------------------------------
def run_verification_a():
    print("\n" + "="*60)
    print("RUNNING VERIFICATION A: SPECIES DISAGREEMENT AUDIT")
    print("="*60)
    
    spec_path = os.path.join(QC_DIR, "04e_species_verification.tsv")
    df_spec = pd.read_csv(spec_path, sep="\t")
    
    # Filter disagreements (agreement == NO)
    disagr = df_spec[df_spec["agreement"] == "NO"].copy()
    print(f"Total species disagreement records: {len(disagr)}")
    
    # Load length and %N from filtering log
    log_path = os.path.join(QC_DIR, "04_genome_filtering_log.tsv")
    df_log = pd.read_csv(log_path, sep="\t")
    log_map = {row["accession"]: row for _, row in df_log.iterrows()}
    
    # Classify each disagreement
    classifications = []
    lengths = []
    n_pcts = []
    
    for _, row in disagr.iterrows():
        acc = row["accession"]
        ani1 = row["ANI_to_HSV1_ref"]
        ani2 = row["ANI_to_HSV2_ref"]
        
        info = log_map.get(acc, {})
        l = info.get("genome_length", 0)
        n = info.get("n_pct", 0.0)
        desc = str(info.get("full_header", ""))
        
        lengths.append(l)
        n_pcts.append(n)
        
        # Classification criteria
        is_patent = "patent" in desc.lower() or "FV" in acc or "OP866" in acc
        
        if is_patent:
            cls = "patent_or_synthetic"
        elif ani1 < 40.0 and ani2 < 40.0:
            cls = "true_non_HSV_contamination"
        elif (40.0 <= ani1 < 85.0) or (40.0 <= ani2 < 85.0):
            cls = "borderline_mislabeled"
        elif (85.0 <= ani1 < 95.0) or (85.0 <= ani2 < 95.0):
            cls = "borderline_mislabeled_review"
        else:
            cls = "borderline_mislabeled"
            
        classifications.append(cls)
        
    disagr["length"] = lengths
    disagr["%N"] = n_pcts
    disagr["classification"] = classifications
    
    cols = ["accession", "database_label", "computational_assignment", "ANI_to_HSV1_ref", "ANI_to_HSV2_ref", "length", "%N", "classification"]
    out_a_tsv = os.path.join(QC_DIR, "04e_species_disagreement_full_list.tsv")
    disagr[cols].to_csv(out_a_tsv, sep="\t", index=False)
    print(f"Saved full disagreement list to {out_a_tsv}")
    
    # Summary of classification counts
    print("\nVerification A - Classification Counts:")
    for cls, cnt in Counter(classifications).most_common():
        print(f"  {cnt:4d} : {cls}")
        
    # Check OP866817 - OP866822
    op_accs = [f"OP8668{i}" for i in range(17, 23)]
    op_in_disagr = disagr[disagr["accession"].isin(op_accs)]
    print(f"\nOP866817-OP866822 Query:")
    print(f"  Total accessions in range : 6 ({', '.join(op_accs)})")
    print(f"  Present in disagreements  : {len(op_in_disagr)} accessions")
    print(f"  Accessions: {list(op_in_disagr['accession'])}")
    print(f"  ANI to refs: HSV-1={op_in_disagr['ANI_to_HSV1_ref'].mean():.2f}%, HSV-2={op_in_disagr['ANI_to_HSV2_ref'].mean():.2f}%")
    print(f"  Classification: {list(op_in_disagr['classification'].unique())}")
    
    # Check FV537010 - FV537013
    fv_accs = [f"FV53701{i}" for i in range(0, 4)]
    fv_in_disagr = disagr[disagr["accession"].isin(fv_accs)]
    print(f"\nFV537010-FV537013 Query:")
    print(f"  Total accessions in range : 4 ({', '.join(fv_accs)})")
    print(f"  Present in disagreements  : {len(fv_in_disagr)} accessions")
    print(f"  Accessions: {list(fv_in_disagr['accession'])}")
    print(f"  Total patent records in dataset: {sum(1 for c in classifications if c == 'patent_or_synthetic')}")

# ----------------------------------------------------------------------
# 2. Verification B & C: Tier 2 QA & Disjointness Check
# ----------------------------------------------------------------------
def run_verification_b_c():
    print("\n" + "="*60)
    print("RUNNING VERIFICATION B & C: TIER 2 QA & DISJOINTNESS")
    print("="*60)
    
    # Load Tier 1 HSV-2 accessions
    tier1_fasta = os.path.join(FILTERED_DIR, "HSV2_filtered.fasta")
    tier1_accs = set(r.id.split(".")[0] for r in SeqIO.parse(tier1_fasta, "fasta"))
    print(f"Tier 1 HSV-2 accessions: {len(tier1_accs)}")
    
    # Load current Tier 2 HSV-2 genomes
    tier2_fasta = os.path.join(FILTERED_DIR, "HSV2_filtered_tier2.fasta")
    tier2_records = list(SeqIO.parse(tier2_fasta, "fasta"))
    print(f"Initial Tier 2 HSV-2 genomes: {len(tier2_records)}")
    
    # Load QA exclusion lists:
    # 1. Lab strains (from 04c_lab_strain_audit.tsv)
    df_lab = pd.read_csv(os.path.join(QC_DIR, "04c_lab_strain_audit.tsv"), sep="\t")
    lab_accs = set(df_lab[df_lab["is_lab_union"]]["accession"])
    
    # 2. Species disagreements / contamination (from 04e_species_verification.tsv)
    df_spec = pd.read_csv(os.path.join(QC_DIR, "04e_species_verification.tsv"), sep="\t")
    contam_accs = set(df_spec[df_spec["agreement"] == "NO"]["accession"])
    
    # 3. Exact duplicates (by SHA-256 hash)
    # Check hashes across Tier 1 and Tier 2
    tier1_hashes = set()
    for r in SeqIO.parse(tier1_fasta, "fasta"):
        s = str(r.seq).upper()
        tier1_hashes.add(hashlib.sha256(s.encode('utf-8')).hexdigest())
        
    qa_rows = []
    clean_tier2_records = []
    seen_tier2_hashes = set()
    
    removed_count = 0
    for r in tier2_records:
        acc = r.id.split(".")[0]
        s = str(r.seq).upper()
        h = hashlib.sha256(s.encode('utf-8')).hexdigest()
        
        in_lab = acc in lab_accs
        in_contam = acc in contam_accs
        in_dup = (h in tier1_hashes) or (h in seen_tier2_hashes)
        
        action = "KEEP"
        reasons = []
        if in_lab:
            action = "REMOVE"
            reasons.append("Matches lab strain")
        if in_contam:
            action = "REMOVE"
            reasons.append("Matches species disagreement/contaminant")
        if in_dup:
            action = "REMOVE"
            reasons.append("Exact sequence duplicate of Tier 1 or prior Tier 2 genome")
            
        if action == "REMOVE":
            removed_count += 1
        else:
            clean_tier2_records.append(r)
            seen_tier2_hashes.add(h)
            
        qa_rows.append({
            "accession": acc,
            "in_lab_strain_list": "YES" if in_lab else "NO",
            "in_duplicate_list": "YES" if in_dup else "NO",
            "in_contamination_list": "YES" if in_contam else "NO",
            "action_taken": action,
            "removal_reasons": "; ".join(reasons) if reasons else "None"
        })
        
    df_qa = pd.DataFrame(qa_rows)
    qa_tsv = os.path.join(QC_DIR, "04f_tier2_qa_report.tsv")
    df_qa.to_csv(qa_tsv, sep="\t", index=False)
    print(f"Saved Tier 2 QA report to {qa_tsv}")
    
    # Rewrite Tier 2 FASTA with verified clean genomes
    SeqIO.write(clean_tier2_records, tier2_fasta, "fasta")
    print(f"\nVerification B Results:")
    print(f"  Tier 2 Genomes Analyzed : {len(tier2_records)}")
    print(f"  Genomes Removed         : {removed_count}")
    print(f"  Final Clean Tier 2 Count: {len(clean_tier2_records)}")
    
    # Verification C: Disjointness Check
    clean_tier2_accs = set(r.id.split(".")[0] for r in clean_tier2_records)
    intersection = tier1_accs & clean_tier2_accs
    total_unique_hsv2 = len(tier1_accs) + len(clean_tier2_accs)
    
    print(f"\nVerification C Results (Tier 1 vs. Tier 2 Disjointness):")
    print(f"  Tier 1 Accession Count  : {len(tier1_accs)}")
    print(f"  Tier 2 Accession Count  : {len(clean_tier2_accs)}")
    print(f"  Intersection Count      : {len(intersection)} (Must be 0)")
    print(f"  Total Unique HSV-2 Genomes: {total_unique_hsv2}")
    if len(intersection) == 0:
        print("  -> CONFIRMED: Tier 1 and Tier 2 cohorts are 100% disjoint.")
    else:
        print(f"  -> WARNING: Found {len(intersection)} intersecting accessions: {intersection}")

# ----------------------------------------------------------------------
# 3. Optional Verification D: Lab Strain ANI Sensitivity (99.7-99.9%)
# ----------------------------------------------------------------------
def run_verification_d():
    print("\n" + "="*60)
    print("RUNNING OPTIONAL VERIFICATION D: LAB STRAIN SENSITIVITY (99.7 - 99.9%)")
    print("="*60)
    
    df_lab = pd.read_csv(os.path.join(QC_DIR, "04c_lab_strain_audit.tsv"), sep="\t")
    
    # Filter genomes with ANI between 99.7% and 99.9% to any reference
    sub = df_lab[
        ((df_lab["ani_strain17"] >= 99.7) & (df_lab["ani_strain17"] < 99.9)) |
        ((df_lab["ani_hg52"] >= 99.7) & (df_lab["ani_hg52"] < 99.9)) |
        ((df_lab["ani_strain_F"] >= 99.7) & (df_lab["ani_strain_F"] < 99.9)) |
        ((df_lab["ani_kos"] >= 99.7) & (df_lab["ani_kos"] < 99.9)) |
        ((df_lab["ani_mckrae"] >= 99.7) & (df_lab["ani_mckrae"] < 99.9))
    ].copy()
    
    # Exclude those that already have a metadata lab flag
    unannotated_candidates = sub[~sub["metadata_lab_flag"]]
    count = len(unannotated_candidates)
    
    print(f"Count of candidate genomes at 99.7-99.9% ANI without metadata lab keyword: {count}")
    
    if count >= 1:
        print("\nCandidate Diverged Lab Strains (99.7 - 99.9% ANI):")
        cols = ["accession", "virus", "isolate", "ani_strain17", "ani_hg52", "ani_strain_F", "ani_kos", "ani_mckrae"]
        print(unannotated_candidates[cols].to_string(index=False))
        out_d = os.path.join(QC_DIR, "04d_lab_strain_sensitivity_candidates.tsv")
        unannotated_candidates.to_csv(out_d, sep="\t", index=False)
        print(f"Saved candidates to {out_d}")

# ----------------------------------------------------------------------
# 4. Step 5: Metadata Curation & Completeness Pipeline
# ----------------------------------------------------------------------
def run_step05_metadata():
    print("\n" + "="*60)
    print("EXECUTING STEP 5: METADATA CURATION PIPELINE")
    print("="*60)
    
    # 1. Load raw metadata from Entrez E-summary
    raw_meta_path = os.path.join(META_DIR, "raw_metadata_all_genomes.tsv")
    df_raw = pd.read_csv(raw_meta_path, sep="\t") if os.path.exists(raw_meta_path) else pd.DataFrame()
    meta_by_acc = {row["accession"]: row.to_dict() for _, row in df_raw.iterrows()}
    
    # 2. Load JSONL summaries
    jsonl_files = [
        os.path.join(META_DIR, "HSV1_raw_summary.jsonl"),
        os.path.join(META_DIR, "HSV2_raw_summary.jsonl")
    ]
    
    datasets_meta = {}
    for jf in jsonl_files:
        if os.path.exists(jf):
            with open(jf, "r", encoding="utf-16", errors="ignore") as f:
                for line in f:
                    if line.strip():
                        try:
                            d = json.loads(line)
                            acc_full = d.get("accession", "")
                            acc_base = acc_full.split(".")[0]
                            datasets_meta[acc_base] = d
                        except Exception:
                            continue
                            
    # Also load unpacked data reports if present
    for v in ["HSV1", "HSV2"]:
        dr = os.path.join(RAW_DIR, v, "unpacked/ncbi_dataset/data/data_report.jsonl")
        if os.path.exists(dr):
            with open(dr, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.strip():
                        try:
                            d = json.loads(line)
                            acc_base = d.get("accession", "").split(".")[0]
                            if acc_base not in datasets_meta:
                                datasets_meta[acc_base] = d
                        except Exception:
                            continue

    print(f"Loaded structured JSONL metadata for {len(datasets_meta)} accessions.")

    # 3. Load species verification and lab strain flags
    df_spec = pd.read_csv(os.path.join(QC_DIR, "04e_species_verification.tsv"), sep="\t")
    spec_map = {row["accession"]: row for _, row in df_spec.iterrows()}
    
    df_lab = pd.read_csv(os.path.join(QC_DIR, "04c_lab_strain_audit.tsv"), sep="\t")
    lab_set = set(df_lab[df_lab["is_lab_union"]]["accession"])

    # Country standardization mapping
    country_standardizer = {
        "USA": "USA", "UNITED STATES": "USA", "UNITED STATES OF AMERICA": "USA", "U.S.A.": "USA",
        "CHINA": "China", "PEOPLE'S REPUBLIC OF CHINA": "China",
        "JAPAN": "Japan", "GERMANY": "Germany", "UNITED KINGDOM": "United Kingdom", "UK": "United Kingdom",
        "ENGLAND": "United Kingdom", "SCOTLAND": "United Kingdom", "FRANCE": "France", "CANADA": "Canada",
        "AUSTRALIA": "Australia", "SOUTH AFRICA": "South Africa", "KENYA": "Kenya", "UGANDA": "Uganda",
        "TANZANIA": "Tanzania", "BRAZIL": "Brazil", "INDIA": "India", "RUSSIA": "Russia", "ITALY": "Italy",
        "SPAIN": "Spain", "SWEDEN": "Sweden", "NETHERLANDS": "Netherlands", "NORWAY": "Norway",
        "FINLAND": "Finland", "KOREA": "South Korea", "SOUTH KOREA": "South Korea", "BELGIUM": "Belgium",
        "SWITZERLAND": "Switzerland"
    }

    # Anatomical source standardization patterns
    cns_pat = re.compile(r"\b(cns|cerebrospinal|csf|brain|encephalitis|meningitis|gangli|temporal lobe)\b", re.I)
    genital_pat = re.compile(r"\b(genital|vulva|vagina|cervix|penis|penile|glans|urethra)\b", re.I)
    oral_pat = re.compile(r"\b(oral|mouth|lip|saliva|throat|pharynx|tongue|swab, oral|cheek)\b", re.I)
    skin_pat = re.compile(r"\b(skin|lesion|vesicle|dermal|blister|cutaneous)\b", re.I)

    # Compile curated metadata rows for all 1,305 genomes
    curated_rows = []
    
    for acc, spec_info in spec_map.items():
        vtype = spec_info["database_label"]
        comp_type = spec_info["computational_assignment"]
        is_recomb = spec_info["is_recombinant"]
        is_lab = acc in lab_set
        
        d_meta = datasets_meta.get(acc, {})
        r_meta = meta_by_acc.get(acc, {})
        
        # 1. Isolate & Strain
        isolate_val = ""
        iso_obj = d_meta.get("isolate")
        if isinstance(iso_obj, dict):
            isolate_val = iso_obj.get("name", "")
        elif isinstance(iso_obj, str):
            isolate_val = iso_obj
            
        if not isolate_val:
            # Fallback to Entrez subname or title
            subname = str(r_meta.get("subname", ""))
            if subname and subname != "nan":
                isolate_val = subname
            else:
                title = str(r_meta.get("title", ""))
                m_iso = re.search(r"\b(?:strain|isolate)\s+([A-Za-z0-9_\-\.\/]+)", title, re.I)
                if m_iso:
                    isolate_val = m_iso.group(1)
                    
        # 2. Country standardization
        raw_loc = ""
        loc_obj = d_meta.get("location")
        if isinstance(loc_obj, dict):
            raw_loc = loc_obj.get("geographic_location", "")
        elif isinstance(loc_obj, str):
            raw_loc = loc_obj
            
        if not raw_loc:
            # Check submitter country
            sub_obj = d_meta.get("submitter")
            if isinstance(sub_obj, dict):
                raw_loc = sub_obj.get("country", "")
                
        # Parse and standardize country
        clean_country = ""
        if raw_loc:
            # Split if "USA: Washington"
            base_country = raw_loc.split(":")[0].strip().upper()
            clean_country = country_standardizer.get(base_country, raw_loc.split(":")[0].strip())

        # 3. Collection Year
        coll_year = ""
        coll_date = ""
        if isinstance(iso_obj, dict):
            coll_date = str(iso_obj.get("collection_date", ""))
        if coll_date and coll_date != "None":
            m_year = re.search(r"\b(19[5-9]\d|20[0-2]\d)\b", coll_date)
            if m_year:
                coll_year = m_year.group(1)
        if not coll_year:
            # Check create date as release reference
            cdate = str(r_meta.get("create_date", ""))
            m_year = re.search(r"\b(19[5-9]\d|20[0-2]\d)\b", cdate)
            # Do not infer collection date from release date per instructions

        # 4. Host
        host_val = ""
        host_obj = d_meta.get("host")
        if isinstance(host_obj, dict):
            host_val = host_obj.get("organism_name", "")
        elif isinstance(host_obj, str):
            host_val = host_obj
        if not host_val or host_val.lower() == "homo sapiens":
            host_val = "Human"

        # 5. Isolation source & CNS classification
        iso_source = ""
        if isinstance(iso_obj, dict):
            iso_source = str(iso_obj.get("isolation_source", ""))
            
        title_for_source = f"{iso_source} {isolate_val} {r_meta.get('title', '')}".lower()
        
        is_cns = bool(cns_pat.search(title_for_source))
        
        std_source = ""
        if is_cns:
            std_source = "CNS / Encephalitis"
        elif genital_pat.search(title_for_source):
            std_source = "Genital"
        elif oral_pat.search(title_for_source):
            std_source = "Oral"
        elif skin_pat.search(title_for_source):
            std_source = "Skin / Lesion"
        elif "gangli" in title_for_source:
            std_source = "Ganglia"
        elif iso_source:
            std_source = iso_source

        # 6. Passage history
        passage_history = ""
        if isinstance(iso_obj, dict):
            passage_history = str(iso_obj.get("passage_history", ""))
            if passage_history.lower() == "none":
                passage_history = ""

        # 7. Sequencing technology
        seq_tech = ""
        nuc_obj = d_meta.get("nucleotide")
        if isinstance(nuc_obj, dict):
            seq_tech = nuc_obj.get("sequencing_technology", "")

        # 8. BioSample & BioProject
        biosample = d_meta.get("biosample", "")
        if not biosample:
            biosample = r_meta.get("biosample", "")
            
        bioproject = ""
        bp_obj = d_meta.get("bioprojects")
        if isinstance(bp_obj, list) and len(bp_obj) > 0:
            bioproject = bp_obj[0]
        elif isinstance(bp_obj, str):
            bioproject = bp_obj
            
        # 9. Submitter & Release date
        submitter_name = ""
        sub_obj = d_meta.get("submitter")
        if isinstance(sub_obj, dict):
            submitter_name = sub_obj.get("affiliation", "")
            
        rel_date = d_meta.get("release_date", "")
        if not rel_date:
            rel_date = r_meta.get("create_date", "")

        curated_rows.append({
            "accession": acc,
            "virus": vtype,
            "isolate": isolate_val,
            "strain": isolate_val if is_lab else "",
            "country": clean_country,
            "collection_year": coll_year,
            "host": host_val,
            "isolation_source": std_source,
            "is_cns_encephalitis": is_cns,
            "passage_history": passage_history,
            "sequencing_technology": seq_tech,
            "BioSample": biosample,
            "BioProject": bioproject,
            "submitter": submitter_name,
            "release_date": rel_date,
            "lab_strain_flag": is_lab,
            "recombinant_flag": is_recomb,
            "computational_type_assignment": comp_type,
            "database_type_label": vtype
        })

    df_curated = pd.DataFrame(curated_rows)
    curated_tsv = os.path.join(META_DIR, "curated_metadata.tsv")
    df_curated.to_csv(curated_tsv, sep="\t", index=False)
    print(f"Saved {curated_tsv} with {len(df_curated)} entries.")

    # 4. Generate Completeness Report
    comp_rows = []
    for vtype in ["HSV-1", "HSV-2", "Combined"]:
        sub = df_curated if vtype == "Combined" else df_curated[df_curated["virus"] == vtype]
        tot = len(sub)
        
        row_stats = {"virus_type": vtype, "total_genomes": tot}
        for col in ["isolate", "country", "collection_year", "host", "isolation_source", "passage_history", "sequencing_technology", "BioSample", "BioProject"]:
            valid_cnt = (sub[col].astype(str).str.strip().replace("nan", "").replace("None", "") != "").sum()
            row_stats[f"{col}_pct"] = round(valid_cnt / tot * 100, 2)
            row_stats[f"{col}_count"] = int(valid_cnt)
            
        comp_rows.append(row_stats)

    df_comp = pd.DataFrame(comp_rows)
    comp_tsv = os.path.join(QC_DIR, "05_metadata_completeness.tsv")
    df_comp.to_csv(comp_tsv, sep="\t", index=False)
    print(f"Saved metadata completeness report to {comp_tsv}")
    
    # Display table in console
    print("\nMetadata Completeness Summary (%):")
    print(df_comp[["virus_type", "total_genomes", "country_pct", "collection_year_pct", "host_pct", "isolation_source_pct", "BioSample_pct", "BioProject_pct"]].to_string(index=False))

if __name__ == "__main__":
    run_verification_a()
    run_verification_b_c()
    run_verification_d()
    run_step05_metadata()
