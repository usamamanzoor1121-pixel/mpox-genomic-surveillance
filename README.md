# Mpox Genomic Surveillance Pipeline

A genomic surveillance pipeline for Monkeypox virus (MPXV) covering clade classification, mutation hotspot analysis, APOBEC3 signature characterization, and recombination breakpoint mapping using Nextclade and NCBI genomic data.

## Key Findings

### 1. Clade Distribution (300 recent MPXV genomes, Aug 2024 - Feb 2025)

- Clade Ib: 253 genomes (84%) - dominant, reflecting the ongoing DRC/Uganda outbreak
- Clade Ia: 31 genomes
- Clade IIb: 14 genomes
- Clade Ib/IIb (inter-clade recombinant): 2 genomes

### 2. Inter-clade Recombinant Identification

Two genomes from Puducherry, India (collected September 2025) were classified as the ad-hoc "Ib/IIb" recombinant clade in Nextclade's reference tree. This clade was added following the description of an inter-clade recombinant MPXV detected in a UK traveller returning from Asia (Pullan et al., 2025, virological.org).

Phylogenetic placement shows these two India genomes are the closest relatives to the UK reference recombinant strain (PP_004DYJ3) in the global tree, independently linking circulating genomes in India to the UK-reported recombinant lineage.

### 3. Recombination Breakpoint Mapping

Analysis of private (recently acquired) mutations in the recombinant genomes shows approximately 85% are concentrated in the terminal ~12kb of the genome (positions 185,000-197,209), overlapping the inverted terminal repeat (ITR) region and genes OPG205 through OPG001_dup. This is consistent with a block recombination event in the structurally variable ITR region rather than gradual point-mutation accumulation.

### 4. APOBEC3 Signature Analysis (private mutations only)

| Clade | Mean APOBEC3 fraction | Interpretation |
|---|---|---|
| IIb | ~0.84 | Heavy APOBEC3 editing, consistent with 2022+ sustained human transmission |
| Ib | ~0.46 | Moderate editing, consistent with ongoing DRC/Uganda transmission chains |
| Ia | ~0.11 | Low editing, consistent with more zoonotic spillover |
| Ib/IIb recombinant | ~0.13 | Low - similar to Ia, NOT Ib. Supports block transfer rather than accumulated edits |

## Pipeline

1. `fetch_mpox_genomes.py` - fetches 300 recent complete MPXV genomes from NCBI via Entrez, with metadata (country, collection date, host)
2. Nextclade CLI - clade assignment, mutation calling, phylogenetic placement against the Nextstrain Mpox all-clades dataset
3. `analyze_mutations.py` - genome-wide mutation hotspot analysis, initial APOBEC3 screening
4. `recombination_analysis.py` - private mutation analysis, recombination breakpoint mapping with gene annotations
5. `finalize_analysis.py` - refined APOBEC3 comparison, phylogenetic placement of recombinant genomes

## Figures

- `apobec3_by_clade.png` - total substitution burden by clade
- `mutation_hotspots.png` - genome-wide mutation density by clade
- `recombination_breakpoint.png` - private mutation positions and gene context at the breakpoint
- `apobec3_private_refined.png` - corrected APOBEC3 fraction using private mutations only
- `phylo_placement_recombinant.png` - phylogenetic placement of recombinant genomes

## Tech Stack

- Data: NCBI Entrez (Biopython)
- Clade classification: Nextclade CLI v3.21, Nextstrain Mpox all-clades dataset
- Analysis: Python, pandas, numpy, Biopython
- Visualization: matplotlib, seaborn

## Local Setup

git clone https://github.com/usamamanzoor1121-pixel/mpox-genomic-surveillance.git
cd mpox-genomic-surveillance
conda create -n mpox-surveillance python=3.11 -y
conda activate mpox-surveillance
pip install -r requirements.txt
python scripts/fetch_mpox_genomes.py
./nextclade dataset get --name nextstrain/mpox/all-clades --output-dir data/raw/mpox_dataset
./nextclade run --input-dataset data/raw/mpox_dataset --output-all nextclade_output/ data/raw/mpox_genomes.fasta
python scripts/analyze_mutations.py
python scripts/recombination_analysis.py
python scripts/finalize_analysis.py

## Limitations

- Genome-wide substitution counts vs reference conflate ancient clade divergence with recent mutation accumulation; the refined analysis uses private (clade-relative) mutations to address this
- QC flags ("mediocre"/"bad") in this dataset primarily reflect missing data (Ns) from low viral load clinical samples, not pipeline errors
- Sample size for clade IIb (n=14) and the recombinant clade (n=2) limits statistical power; findings should be treated as exploratory

## Author

Usama Manzoor
JSMU Diagnostic Laboratory and Blood Bank, Karachi, Pakistan
GitHub: https://github.com/usamamanzoor1121-pixel
LinkedIn: https://www.linkedin.com/in/usama-manzoor-042595182/
