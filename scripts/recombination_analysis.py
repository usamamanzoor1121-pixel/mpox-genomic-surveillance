"""
Script 4: Recombination breakpoint analysis for Ib/IIb genomes
+ Refined APOBEC3 signature analysis using private (recently acquired) mutations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
from Bio import SeqIO

OUTPUT_DIR = "figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

GENOME_LENGTH = 197209

# ── Load reference sequence ────────────────────────────────────────────────────

ref_record = next(SeqIO.parse("data/raw/mpox_dataset/reference.fasta", "fasta"))
ref_seq = str(ref_record.seq).upper()

# ── Load gene annotations ──────────────────────────────────────────────────────

def load_gene_annotations(gff_path):
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            start, end = int(parts[3]), int(parts[4])
            attrs = parts[8]
            m = re.search(r'(?:gene_name|Name|gene)=([^;]+)', attrs)
            name = m.group(1) if m else "unknown"
            genes.append({"name": name, "start": start, "end": end})
    return pd.DataFrame(genes)

genes_df = load_gene_annotations("data/raw/mpox_dataset/genome_annotation.gff3")
print(f"Loaded {len(genes_df)} gene annotations")

# Genes near the genome terminus (last 20kb)
terminal_genes = genes_df[genes_df["end"] >= GENOME_LENGTH - 20000].sort_values("start")
print(f"\nGenes in terminal 20kb:")
print(terminal_genes.to_string(index=False))

# ── Load Nextclade results ──────────────────────────────────────────────────────

nc = pd.read_csv("nextclade_output/nextclade.tsv", sep="\t", low_memory=False)
nc["accession"] = nc["seqName"].str.replace(r"\s.*", "", regex=True)

# ── Parse substitution strings ──────────────────────────────────────────────────

def parse_subs(sub_str):
    if pd.isna(sub_str) or str(sub_str) in ("", "nan"):
        return []
    subs = []
    for s in str(sub_str).split(","):
        if len(s) < 3:
            continue
        ref_base, alt_base = s[0], s[-1]
        try:
            pos = int(s[1:-1])
        except ValueError:
            continue
        subs.append((ref_base, pos, alt_base))
    return subs

def get_private_mutations(row):
    all_subs = []
    for col in ["privateNucMutations.reversionSubstitutions",
                 "privateNucMutations.labeledSubstitutions",
                 "privateNucMutations.unlabeledSubstitutions"]:
        all_subs.extend(parse_subs(row.get(col, "")))
    return all_subs

def is_apobec3(ref_base, pos, alt_base, ref_seq):
    idx = pos - 1
    if ref_base == "C" and alt_base == "T" and idx > 0:
        return ref_seq[idx-1] == "T"
    if ref_base == "G" and alt_base == "A" and idx < len(ref_seq)-1:
        return ref_seq[idx+1] == "A"
    return False

# ── Refined APOBEC3 analysis on PRIVATE mutations only ──────────────────────────

print("\nComputing private-mutation APOBEC3 fractions...")

priv_results = []
for _, row in nc.iterrows():
    priv = get_private_mutations(row)
    n_total = len(priv)
    n_apobec = sum(is_apobec3(r, p, a, ref_seq) for r, p, a in priv)
    priv_results.append({
        "accession": row["accession"],
        "clade": row["clade"],
        "n_private": n_total,
        "n_apobec3_private": n_apobec,
        "apobec3_fraction_private": n_apobec / n_total if n_total > 0 else np.nan,
        "private_positions": [p for _, p, _ in priv],
    })

priv_df = pd.DataFrame(priv_results)

print("\n" + "="*55)
print("PRIVATE MUTATION APOBEC3 FRACTION BY CLADE")
print("="*55)
summary = priv_df.groupby("clade").agg(
    n_genomes=("accession","count"),
    mean_n_private=("n_private","mean"),
    mean_apobec3_fraction=("apobec3_fraction_private","mean"),
).round(3)
print(summary.to_string())

# ── Figure: recombination breakpoint map ────────────────────────────────────────

recomb_rows = priv_df[priv_df["clade"]=="Ib/IIb"]
ib_sample = priv_df[priv_df["clade"]=="Ib"].sample(min(20, (priv_df["clade"]=="Ib").sum()), random_state=42)

fig, axes = plt.subplots(2, 1, figsize=(16, 8), gridspec_kw={"height_ratios":[3,1]})

# Top: full genome view
y = 0
for _, row in recomb_rows.iterrows():
    axes[0].scatter(row["private_positions"], [y]*len(row["private_positions"]),
                     color="red", s=30, label="Ib/IIb recombinant" if y==0 else "")
    y += 1

for _, row in ib_sample.iterrows():
    axes[0].scatter(row["private_positions"], [y]*len(row["private_positions"]),
                     color="gray", s=10, alpha=0.5,
                     label="Typical Ib genome" if y==2 else "")
    y += 1

axes[0].set_xlim(0, GENOME_LENGTH)
axes[0].set_xlabel("Genome position (bp)")
axes[0].set_ylabel("Genome (each row = one isolate)")
axes[0].set_title("Private Mutation Positions: Ib/IIb Recombinants vs Typical Ib Genomes")
axes[0].legend(loc="upper left")
axes[0].axvspan(185000, GENOME_LENGTH, color="yellow", alpha=0.2)

# Bottom: zoomed terminal region with gene annotations
zoom_start = 180000
for _, row in recomb_rows.iterrows():
    positions_zoom = [p for p in row["private_positions"] if p >= zoom_start]
    axes[1].scatter(positions_zoom, [1]*len(positions_zoom), color="red", s=50, zorder=3)

for _, gene in terminal_genes.iterrows():
    if gene["end"] >= zoom_start:
        axes[1].barh(0, gene["end"]-max(gene["start"],zoom_start),
                      left=max(gene["start"],zoom_start), height=0.4,
                      color="steelblue", alpha=0.6)
        axes[1].text(max(gene["start"],zoom_start), 0.3, gene["name"],
                      fontsize=7, rotation=45)

axes[1].set_xlim(zoom_start, GENOME_LENGTH)
axes[1].set_ylim(-0.5, 1.5)
axes[1].set_xlabel(f"Genome position (bp) — zoomed to last {GENOME_LENGTH-zoom_start}bp")
axes[1].set_yticks([])
axes[1].set_title("Recombination Breakpoint Region — Gene Context")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/recombination_breakpoint.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {OUTPUT_DIR}/recombination_breakpoint.png")

# ── Save data ────────────────────────────────────────────────────────────────────

priv_df.drop(columns=["private_positions"]).to_csv("data/processed/private_mutations_summary.csv", index=False)
print("\nDone.")
