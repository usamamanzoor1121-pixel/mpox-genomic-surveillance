"""
Script 5: Final figures
- Refined APOBEC3 fraction (private mutations only) - corrected comparison
- Phylogenetic placement of Ib/IIb recombinants
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
import os

OUTPUT_DIR = "figures"

# ── Refined APOBEC3 figure ──────────────────────────────────────────────────────

priv_df = pd.read_csv("data/processed/private_mutations_summary.csv")

fig, ax = plt.subplots(figsize=(8, 6))
order = ["Ia", "Ib", "Ib/IIb", "IIb"]
sns.boxplot(data=priv_df, x="clade", y="apobec3_fraction_private", order=order, ax=ax)
sns.stripplot(data=priv_df, x="clade", y="apobec3_fraction_private", order=order,
               color="black", alpha=0.4, size=4, ax=ax)
ax.set_title("APOBEC3 Signature Fraction — Private (Recently Acquired) Mutations Only")
ax.set_ylabel("Fraction of private mutations matching APOBEC3 signature")
ax.set_xlabel("Clade")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/apobec3_private_refined.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUTPUT_DIR}/apobec3_private_refined.png")

# ── Phylogenetic placement of recombinants ─────────────────────────────────────

print("Loading tree...")
tree = Phylo.read("nextclade_output/nextclade.nwk", "newick")

nc = pd.read_csv("nextclade_output/nextclade.tsv", sep="\t", low_memory=False)
nc["accession"] = nc["seqName"].str.replace(r"\s.*", "", regex=True)
clade_map = dict(zip(nc["accession"], nc["clade"]))

recomb_names = nc[nc["clade"]=="Ib/IIb"]["seqName"].tolist()
print(f"Recombinant tip names: {recomb_names}")

# Find terminals matching recombinant names
terminals = tree.get_terminals()
recomb_terminals = [t for t in terminals if any(r in t.name for r in recomb_names)]
print(f"Found {len(recomb_terminals)} recombinant terminals in tree")

if recomb_terminals:
    # Get common ancestor and its clade for context
    target = recomb_terminals[0]
    path = tree.get_path(target)
    print(f"Path depth from root: {len(path)}")
    
    # Find a subtree: go up 2-3 levels from the target tip
    if len(path) >= 3:
        subtree_root = path[-3]
    else:
        subtree_root = tree.root
    
    subtree_terminals = subtree_root.get_terminals()
    print(f"Subtree contains {len(subtree_terminals)} tips")
    
    # Limit to reasonable size for plotting
    if len(subtree_terminals) > 30:
        # go down one more level
        subtree_root = path[-2]
        subtree_terminals = subtree_root.get_terminals()
        print(f"Reduced subtree contains {len(subtree_terminals)} tips")
    
    # Build a small tree for plotting
    import copy
    plot_tree = copy.deepcopy(subtree_root)
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(subtree_terminals)*0.3)))
    
    def get_label_color(clade):
        colors = {"Ia":"#2a9d8f","Ib":"#e76f51","IIb":"#457b9d","Ib/IIb":"#e63946"}
        return colors.get(clade, "gray")
    
    def label_func(clade_obj):
        if clade_obj.is_terminal():
            acc = clade_obj.name.split()[0] if clade_obj.name else ""
            cl = clade_map.get(acc, "?")
            return f"{acc} [{cl}]"
        return ""
    
    Phylo.draw(plot_tree, axes=ax, label_func=label_func, do_show=False)
    
    # Color recombinant tip labels red
    for text in ax.texts:
        if "Ib/IIb" in text.get_text():
            text.set_color("red")
            text.set_fontweight("bold")
    
    ax.set_title("Phylogenetic Placement — Ib/IIb Recombinant Genomes (red)")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/phylo_placement_recombinant.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUTPUT_DIR}/phylo_placement_recombinant.png")
else:
    print("Could not locate recombinant tips in tree - skipping placement plot")

print("\nDone.")
