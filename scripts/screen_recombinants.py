"""
Script 6: Screen all genomes for terminal-region private mutation enrichment
- Tests whether each genome's private mutations are anomalously concentrated
  in the terminal ITR region (>=185,000 bp), as seen in the known Ib/IIb recombinants
- Uses binomial test against expected uniform distribution
"""

import pandas as pd
import numpy as np
from scipy.stats import binomtest
import matplotlib.pyplot as plt
import os

OUTPUT_DIR = "figures"
GENOME_LENGTH = 197209
TERMINAL_START = 185000
TERMINAL_LENGTH = GENOME_LENGTH - TERMINAL_START + 1
P_TERMINAL = TERMINAL_LENGTH / GENOME_LENGTH  # baseline probability if uniform

print(f"Terminal region: {TERMINAL_START}-{GENOME_LENGTH} ({TERMINAL_LENGTH} bp)")
print(f"Baseline probability (uniform): {P_TERMINAL:.4f}")

# ── Load data ───────────────────────────────────────────────────────────────────

nc = pd.read_csv("nextclade_output/nextclade.tsv", sep="\t", low_memory=False)
meta = pd.read_csv("data/raw/mpox_metadata.csv")
nc["accession"] = nc["seqName"].str.replace(r"\s.*", "", regex=True)
df = nc.merge(meta, on="accession", how="left")

# ── Parse all private mutations ──────────────────────────────────────────────────

def parse_subs(sub_str):
    if pd.isna(sub_str) or str(sub_str) in ("", "nan"):
        return []
    subs = []
    for s in str(sub_str).split(","):
        if len(s) < 3:
            continue
        try:
            pos = int(s[1:-1])
        except ValueError:
            continue
        subs.append(pos)
    return subs

def get_private_positions(row):
    positions = []
    for col in ["privateNucMutations.reversionSubstitutions",
                 "privateNucMutations.labeledSubstitutions",
                 "privateNucMutations.unlabeledSubstitutions"]:
        positions.extend(parse_subs(row.get(col, "")))
    return positions

print("Computing terminal enrichment for all genomes...")

results = []
for _, row in df.iterrows():
    positions = get_private_positions(row)
    k = len(positions)
    m = sum(1 for p in positions if p >= TERMINAL_START)
    
    if k >= 3:  # need minimum mutations for meaningful test
        pval = binomtest(m, k, P_TERMINAL, alternative="greater").pvalue
    else:
        pval = np.nan
    
    results.append({
        "accession": row["accession"],
        "clade": row["clade"],
        "country": row["country"],
        "collection_date": row["collection_date"],
        "qc_status": row["qc.overallStatus"],
        "n_private": k,
        "n_terminal": m,
        "terminal_fraction": m/k if k > 0 else np.nan,
        "pvalue": pval,
    })

results_df = pd.DataFrame(results)

# ── Identify candidates ─────────────────────────────────────────────────────────

# Exclude already-known recombinants
known_recomb = ["OZ478275.1", "OZ478273.1"]
candidates = results_df[
    (~results_df["accession"].isin(known_recomb)) &
    (results_df["n_private"] >= 3) &
    (results_df["pvalue"] < 0.05)
].sort_values("pvalue")

print("\n" + "="*60)
print(f"Genomes with significant terminal-region enrichment (p<0.05)")
print(f"excluding the 2 known recombinants")
print("="*60)
print(f"Total candidates: {len(candidates)}")
if len(candidates) > 0:
    print(candidates.to_string(index=False))
else:
    print("No additional candidates found at p<0.05")

# Show the known recombinants' values for reference
print("\n" + "="*60)
print("Known recombinants (for reference)")
print("="*60)
known = results_df[results_df["accession"].isin(known_recomb)]
print(known.to_string(index=False))

# ── Distribution summary ──────────────────────────────────────────────────────

print("\n" + "="*60)
print("Overall distribution of terminal_fraction (genomes with n_private>=3)")
print("="*60)
valid = results_df[results_df["n_private"]>=3]
print(valid["terminal_fraction"].describe().to_string())

# ── Figure ─────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(10,6))
valid_other = valid[~valid["accession"].isin(known_recomb)]
ax.scatter(valid_other["n_private"], valid_other["terminal_fraction"],
           alpha=0.4, s=20, color="gray", label="Other genomes")
known_plot = valid[valid["accession"].isin(known_recomb)]
ax.scatter(known_plot["n_private"], known_plot["terminal_fraction"],
           color="red", s=100, label="Known Ib/IIb recombinants", zorder=5)
if len(candidates) > 0:
    cand_plot = valid[valid["accession"].isin(candidates["accession"])]
    ax.scatter(cand_plot["n_private"], cand_plot["terminal_fraction"],
               color="orange", s=80, label="New candidates (p<0.05)", zorder=4)

ax.axhline(P_TERMINAL, color="black", linestyle="--", alpha=0.5, label="Expected (uniform)")
ax.set_xlabel("Number of private mutations")
ax.set_ylabel("Fraction in terminal region (>=185,000 bp)")
ax.set_title("Terminal-region Private Mutation Enrichment Screen")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/recombinant_screen.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {OUTPUT_DIR}/recombinant_screen.png")

results_df.to_csv("data/processed/recombinant_screen.csv", index=False)
