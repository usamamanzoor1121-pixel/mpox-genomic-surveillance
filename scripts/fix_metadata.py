"""Quick fix: re-extract metadata checking both country and geo_loc_name fields"""
from Bio import Entrez, SeqIO
import pandas as pd
import time, os

Entrez.email = "usamamanzoor1121@gmail.com"

df_old = pd.read_csv("data/raw/mpox_metadata.csv")
ids = df_old["accession"].tolist()

records_out = []
for i in range(0, len(ids), 20):
    batch = ids[i:i+20]
    print(f"  {i+1}-{i+len(batch)} of {len(ids)}")
    handle = Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text")
    for rec in SeqIO.parse(handle, "genbank"):
        country, coldate, host, isolate = "Unknown","Unknown","Unknown","Unknown"
        for f in rec.features:
            if f.type == "source":
                q = f.qualifiers
                country = q.get("geo_loc_name", q.get("country", ["Unknown"]))[0]
                coldate = q.get("collection_date", ["Unknown"])[0]
                host    = q.get("host", ["Unknown"])[0]
                isolate = q.get("isolate", ["Unknown"])[0]
        records_out.append({"accession": rec.id, "country": country,
                             "collection_date": coldate, "host": host, "isolate": isolate})
    handle.close()
    time.sleep(0.5)

df_new = pd.DataFrame(records_out)
df_new.to_csv("data/raw/mpox_metadata.csv", index=False)
print("Done. Country distribution:")
print(df_new["country"].value_counts().head(15))
