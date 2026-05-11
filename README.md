# m6A GMUCT Cleavage Analysis

This repository contains analyses performed for BIOL 3999 investigating the relationship between m6A methylation and subcellular mRNA cleavage bias in *Arabidopsis thaliana*.

## Project Overview

Two complementary analyses were performed:

1. Coordinate-level overlap analysis
   - Tested whether cleavage sites from NAR Supplementary Table 2 overlap high-confidence Col-0 m6A peak coordinates.
   - Included observed vs expected overlap calculations and enrichment testing.

2. Gene-level flower CS m6A analysis
   - Compared flower control sample (flower CS) m6A peak annotations between nucleoplasmic and cytoplasmic mRNA cleavage-associated genes.
   - Included Fisher’s exact tests and visualization of methylation frequencies.

## Main Findings

- mRNA cleavage sites showed significant enrichment within Col-0 m6A peak regions relative to random expectation.
- Direct coordinate overlaps were strongly nucleoplasmic-enriched.
- Flower CS m6A peaks were more common among cytoplasmic cleavage-associated genes at the gene level.

## Repository Structure

- `src/`
  - Python analysis scripts
- `figures/`
  - Generated figures used in the writeup
- `data/`
  - Input and processed datasets
- `results/`
  - Tables and output summaries

## Software

Analyses were performed in Python using:
- pandas
- PyRanges
- NumPy
- SciPy
- matplotlib
