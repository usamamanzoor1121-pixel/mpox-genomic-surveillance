"""
Script 1: Fetch MPXV (Monkeypox virus) genomes from NCBI
Searches for complete genomes across clades and years
Downloads FASTA sequences + metadata (collection date, country, host)
"""

from Bio import Entrez, SeqIO
import pandas as pd
import time
import os
import re

Entrez.email = "usamamanzoor1121@gmail.com"  # required by NCBI, any valid format works

OUTPUT_DIR = "data/raw"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Search query ───────────────────────────────────────────────────────────────
# Complete MPXV genomes, broad date range to capture clade Ia/Ib/IIa/IIb diversity
# Genome size filter: MPXV genomes are ~190,000-210,000 bp

SEARCH_TERM = (
    '("Monkeypox virus"[Organism]) '
    'AND complete genome[Title] '
    'AND 180000:215000[SLEN]'
)

MAX_RECORDS = 300  # manageable size at ~197kb each = ~60MB total

def search_genomes():
    print(f"Searching NCBI: {SEARCH_TERM}")
    handle = Entrez.esearch(
        db="nucleotide",
        term=SEARCH_TERM,
        retmax=MAX_RECORDS,
        sort="date",  # most recent first - captures latest outbreak strains
    )
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    print(f"Found {len(ids)} genome records (showing most recent {MAX_RECORDS})")
    return ids


def fetch_metadata_and_sequences(ids, batch_size=20):
    """Fetch GenBank records in batches - extract metadata + sequence."""
    
    metadata_records = []
    fasta_path = os.path.join(OUTPUT_DIR, "mpox_genomes.fasta")
    
    with open(fasta_path, "w") as fasta_out:
        for i in range(0, len(ids), batch_size):
            batch_ids = ids[i:i + batch_size]
            print(f"  Fetching records {i+1}-{i+len(batch_ids)} of {len(ids)}...")
            
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=batch_ids,
                    rettype="gb",
                    retmode="text",
                )
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
            except Exception as e:
                print(f"    Error fetching batch: {e}")
                time.sleep(2)
                continue
            
            for rec in records:
                # Write FASTA
                fasta_out.write(f">{rec.id}\n{str(rec.seq)}\n")
                
                # Extract metadata from source feature
                country = "Unknown"
                collection_date = "Unknown"
                host = "Unknown"
                isolate = "Unknown"
                
                for feature in rec.features:
                    if feature.type == "source":
                        qualifiers = feature.qualifiers
                        if "country" in qualifiers:
                            country = qualifiers["country"][0]
                        if "collection_date" in qualifiers:
                            collection_date = qualifiers["collection_date"][0]
                        if "host" in qualifiers:
                            host = qualifiers["host"][0]
                        if "isolate" in qualifiers:
                            isolate = qualifiers["isolate"][0]
                
                metadata_records.append({
                    "accession": rec.id,
                    "description": rec.description,
                    "length": len(rec.seq),
                    "country": country,
                    "collection_date": collection_date,
                    "host": host,
                    "isolate": isolate,
                })
            
            time.sleep(0.5)  # be polite to NCBI
    
    return pd.DataFrame(metadata_records)


def main():
    print("=" * 55)
    print("MPXV Genome Fetch — NCBI")
    print("=" * 55)
    
    ids = search_genomes()
    if not ids:
        print("No records found.")
        return
    
    print("\nFetching metadata and sequences...")
    df = fetch_metadata_and_sequences(ids)
    
    # Save metadata
    meta_path = os.path.join(OUTPUT_DIR, "mpox_metadata.csv")
    df.to_csv(meta_path, index=False)
    
    print("\n" + "=" * 55)
    print("SUMMARY")
    print("=" * 55)
    print(f"Total genomes fetched: {len(df)}")
    print(f"\nCountry distribution:")
    print(df["country"].value_counts().head(15).to_string())
    print(f"\nCollection date range (sample):")
    print(df["collection_date"].value_counts().head(10).to_string())
    print(f"\nSaved:")
    print(f"  {OUTPUT_DIR}/mpox_genomes.fasta")
    print(f"  {OUTPUT_DIR}/mpox_metadata.csv")


if __name__ == "__main__":
    main()
