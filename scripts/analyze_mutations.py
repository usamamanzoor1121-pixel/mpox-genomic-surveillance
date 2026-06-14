"""
Script 3: Mutation hotspot and APOBEC3 signature analysis
- Merge Nextclade results with metadata
- Identify APOBEC3-consistent mutations (TC->TT and GA->AA signatures)
- Compare APOBEC3 burden across clades
- Genome-wide mutation hotspot visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from Bio import SeqIO

OUTPUT_DIR = "figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load reference sequence ────────────────────────────────────────────────────

print("Loading reference genome...")
ref_record = next(SeqIO.parse("data/raw/mpox_dataset/reference.fasta", "fasta"))
ref_seq = str(ref_record.seq).upper()
print(f"Reference length: {len(ref_seq):,} bp")

# ── Load Nextclade results and metadata ────────────────────────────────────────

print("Loading Nextclade results...")
nc = pd.read_csv("nextclade_output/nextclade.tsv", sep="\t", low_memory=False)
meta = pd.read_csv("data/raw/mpox_metadata.csv")

# Merge on accession (seqName matches accession)
nc["accession"] = nc["seqName"].str.replace(r"\s.*", "", regex=True)
df = nc.merge(meta, on="accession", how="left")
print(f"Merged dataset: {len(df)} genomes")

# ── Parse substitutions and classify APOBEC3 signature ─────────────────────────

def parse_substitutions(sub_str):
    """Parse 'T128C,G475A,...' into list of (ref, pos, alt) tuples. Position is 1-indexed."""
    if pd.isna(sub_str) or sub_str == "":
        return []
    subs = []
    for s in sub_str.split(","):
        if len(s) < 3:
            continue
        ref_base = s[0]
        alt_base = s[-1]
        try:
            pos = int(s[1:-1])
        except ValueError:
            continue
        subs.append((ref_base, pos, alt_base))
    return subs

def is_apobec3_signature(ref_base, pos, alt_base, ref_seq):
    """
    APOBEC3 signature: C->T in TC context (TC->TT), or G->A in GA context (GA->AA)
    pos is 1-indexed
    """
    idx = pos - 1  # convert to 0-indexed
    if ref_base == "C" and alt_base == "T":
        # check preceding base
        if idx > 0 and ref_seq[idx-1] == "T":
            return True
    elif ref_base == "G" and alt_base == "A":
        # check following base
        if idx < len(ref_seq)-1 and ref_seq[idx+1] == "A":
            return True
    return False

print("Analyzing mutations for APOBEC3 signatures...")

results = []
all_positions = []  # for hotspot plot

for _, row in df.iterrows():
    subs = parse_substitutions(row.get("substitutions", ""))
    n_total = len(subs)
    n_apobec = 0
    
    for ref_base, pos, alt_base in subs:
        all_positions.append({
            "accession": row["accession"],
            "clade": row["clade"],
            "position": pos,
            "ref": ref_base,
            "alt": alt_base,
        })
        if is_apobec3_signature(ref_base, pos, alt_base, ref_seq):
            n_apobec += 1
    
    apobec_fraction = n_apobec / n_total if n_total > 0 else 0
    
    results.append({
        "accession": row["accession"],
        "clade": row["clade"],
        "country": row["country"],
        "collection_date": row["collection_date"],
        "n_substitutions": n_total,
        "n_apobec3": n_apobec,
        "apobec3_fraction": apobec_fraction,
    })

results_df = pd.DataFrame(results)
mutations_df = pd.DataFrame(all_positions)

# ── Summary by clade ───────────────────────────────────────────────────────────

print("\n" + "="*55)
print("APOBEC3 SIGNATURE BURDEN BY CLADE")
print("="*55)
summary = results_df.groupby("clade").agg(
    n_genomes=("accession", "count"),
    mean_substitutions=("n_substitutions", "mean"),
    mean_apobec3_fraction=("apobec3_fraction", "mean"),
).round(3)
print(summary.to_string())

# ── Highlight the Ib/IIb recombinants ──────────────────────────────────────────

print("\n" + "="*55)
print("Ib/IIb RECOMBINANT GENOMES - DETAIL")
print("="*55)
recomb = results_df[results_df["clade"]=="Ib/IIb"]
print(recomb.to_string(index=False))

# ── Save processed data ─────────────────────────────────────────────────────────

os.makedirs("data/processed", exist_ok=True)
results_df.to_csv("data/processed/apobec3_summary.csv", index=False)
mutations_df.to_csv("data/processed/all_mutations.csv", index=False)
df.to_csv("data/processed/mpox_merged.csv", index=False)

# ── Figure 1: APOBEC3 fraction by clade (boxplot) ──────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sns.boxplot(data=results_df, x="clade", y="apobec3_fraction", ax=axes[0])
axes[0].set_title("APOBEC3-consistent Mutation Fraction by Clade")
axes[0].set_ylabel("Fraction of substitutions matching APOBEC3 signature")
axes[0].set_xlabel("Clade")

sns.boxplot(data=results_df, x="clade", y="n_substitutions", ax=axes[1])
axes[1].set_title("Total Substitutions by Clade")
axes[1].set_ylabel("Number of substitutions vs reference")
axes[1].set_xlabel("Clade")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/apobec3_by_clade.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {OUTPUT_DIR}/apobec3_by_clade.png")

# ── Figure 2: Genome-wide mutation hotspot ─────────────────────────────────────

fig, ax = plt.subplots(figsize=(16, 5))

bins = np.arange(0, len(ref_seq), 1000)
for clade, color in zip(["Ia","Ib","IIb","Ib/IIb"], ["#2a9d8f","#e76f51","#457b9d","#e63946"]):
    sub = mutations_df[mutations_df["clade"]==clade]
    if len(sub) == 0:
        continue
    ax.hist(sub["position"], bins=bins, alpha=0.5, label=f"Clade {clade}", color=color)

ax.set_xlabel("Genome position (bp)")
ax.set_ylabel("Mutation count (1kb bins)")
ax.set_title("Genome-wide Mutation Hotspots by Clade")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/mutation_hotspots.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUTPUT_DIR}/mutation_hotspots.png")

print("\nDone.")
