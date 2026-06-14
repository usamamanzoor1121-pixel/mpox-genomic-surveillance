# Mpox Genomic Surveillance Pipeline

A genomic surveillance pipeline for Monkeypox virus (MPXV) covering clade classification, mutation hotspot analysis, APOBEC3 signature characterization, and recombination breakpoint screening using Nextclade and NCBI genomic data.

## Overview

This pipeline fetches recent MPXV genomes from NCBI, runs clade classification via Nextclade, and performs downstream analysis of mutation patterns including a statistical screen for inter-clade recombination signatures.

## Key Findings

### 1. Clade Distribution (300 recent MPXV genomes, Aug 2024 - Feb 2025)

- Clade Ib: 253 genomes (84%) - dominant, reflecting the ongoing DRC/Uganda outbreak
- Clade Ia: 31 genomes
- Clade IIb: 14 genomes
- Clade Ib/IIb (inter-clade recombinant): 2 genomes

### 2. Validation Against the WHO-Confirmed Ib/IIb Recombinant (Feb 2026)

The two genomes classified as "Ib/IIb" (accessions OZ478275.1 and OZ478273.1, Puducherry, India, September 2025) correspond to the recombinant MPXV case confirmed by WHO in February 2026 - the earliest known detection of an inter-clade Ib/IIb recombinant, alongside a UK case in a traveller returning from the Asia Pacific region (Pullan et al., 2025). This pipeline independently reproduced the correct clade classification for this case using public NCBI data and the Nextclade reference dataset, confirming the tooling works correctly on a real, high-profile surveillance event.

### 3. Recombination Breakpoint Mapping

Analysis of private (recently acquired) mutations in the recombinant genomes shows approximately 82% are concentrated in the terminal ~12kb of the genome (positions 185,000-197,209), overlapping the inverted terminal repeat (ITR) region and genes OPG205 through OPG001_dup. This is consistent with a block recombination event in the structurally variable ITR region rather than gradual point-mutation accumulation.

### 4. Screening for Additional Undetected Recombinants

WHO noted the UK and India cases fell ill several weeks apart with the same recombinant strain, raising the possibility of additional undetected cases. Using a binomial test for terminal-region private mutation enrichment, this pipeline screened all 300 genomes:

- The two known recombinants: p = 2.17e-32 (82% of 39 private mutations in a region representing 6.2% of the genome)
- All other 294 genomes: no candidates reached p < 0.05; median terminal fraction = 0

No additional candidates were detected in this sample. This sample is dominated by Uganda/DRC clade Ib genomes from Aug 2024 - Feb 2025, geographically and temporally distant from the suspected origin of the recombinant (South/Southeast Asia or Arabian Peninsula, mid-to-late 2025). Broader sampling targeting that region and timeframe would be needed for a more conclusive screen.

### 5. APOBEC3 Signature Analysis (private mutations only)

| Clade | Mean APOBEC3 fraction | Interpretation |
|---|---|---|
| IIb | ~0.84 | Heavy APOBEC3 editing, consistent with 2022+ sustained human transmission |
| Ib | ~0.46 | Moderate editing, consistent with ongoing DRC/Uganda transmission chains |
| Ia | ~0.11 | Low editing, consistent with more zoonotic spillover |
| Ib/IIb recombinant | ~0.13 | Low - similar to Ia, not Ib. Supports block transfer rather than accumulated edits |

## Pipeline

1. `fetch_mpox_genomes.py` - fetches recent complete MPXV genomes from NCBI via Entrez, with metadata (country, collection date, host)
2. Nextclade CLI - clade assignment, mutation calling, phylogenetic placement against the Nextstrain Mpox all-clades dataset
3. `analyze_mutations.py` - genome-wide mutation hotspot analysis, initial APOBEC3 screening
4. `recombination_analysis.py` - private mutation analysis, recombination breakpoint mapping with gene annotations
5. `finalize_analysis.py` - refined APOBEC3 comparison, phylogenetic placement of recombinant genomes
6. `screen_recombinants.py` - statistical screen for additional genomes with terminal-region recombination signatures

## Figures

- `apobec3_by_clade.png` - total substitution burden by clade
- `mutation_hotspots.png` - genome-wide mutation density by clade
- `recombination_breakpoint.png` - private mutation positions and gene context at the breakpoint
- `apobec3_private_refined.png` - corrected APOBEC3 fraction using private mutations only
- `phylo_placement_recombinant.png` - phylogenetic placement of recombinant genomes
- `recombinant_screen.png` - genome-wide screen for terminal-region recombination signatures

## Tech Stack

- Data: NCBI Entrez (Biopython)
- Clade classification: Nextclade CLI v3.21, Nextstrain Mpox all-clades dataset
- Analysis: Python, pandas, numpy, scipy, Biopython
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
python scripts/screen_recombinants.py

## Limitations

- Genome-wide substitution counts vs reference conflate ancient clade divergence with recent mutation accumulation; the refined analysis uses private (clade-relative) mutations to address this
- QC flags ("mediocre"/"bad") in this dataset primarily reflect missing data (Ns) from low viral load clinical samples, not pipeline errors
- Sample size for clade IIb (n=14) and the recombinant clade (n=2) limits statistical power
- The recombinant screening sample is geographically and temporally biased toward the Uganda/DRC Ib outbreak; a targeted search of South/Southeast Asian and Middle Eastern genomes from mid-to-late 2025 would be needed for a more conclusive screen for additional cases

## Author

Usama Manzoor
JSMU Diagnostic Laboratory and Blood Bank, Karachi, Pakistan
GitHub: https://github.com/usamamanzoor1121-pixel
LinkedIn: https://www.linkedin.com/in/usama-manzoor-042595182/
