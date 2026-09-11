# Keloid-scRNAseq

This repository contains the analysis code for the single-cell RNA sequencing (scRNA-seq) and spatial transcriptomics study of keloid.

## Overview

The study investigates whether keloid fibroblasts retain a persistent late-remodeling state (Day 30-like) and lose Skin-like homeostasis, compared with normal scar and other chronic wound controls.

## Repository Structure

- `code/` – Analysis scripts ordered by stage (Stage 1 to Stage 13), plus the final publication figure script.
- `README.md` – This file.

## Analysis Pipeline

Scripts are numbered by stage and should be run in order:

1. Stage 1–2 – Raw data audit, matrix structure and sample mapping.
2. Stage 3 – Normal-wound reference QC and pseudobulk construction (GSE241132).
3. Stage 4 – Donor-held-out wound-state model and frozen signatures.
4. Stage 5 – Discovery cohort mapping (GSE163973).
5. Stage 6A–6B – Independent replication cohort (GSE181316).
6. Stage 7 – External keloid cohort support (GSE220300).
7. Stage 8 – Normal-wound spatial validation (GSE241124).
8. Stage 9 – Venous-ulcer disease control (GSE265972).
9. Stage 10 – Cross-modal exploratory support (GSE181297).
10. Stage 11 – Final stratified exact integration of core cohorts.
11. Stage 12 – Cross-cohort consensus genes and pathway interpretation.
12. Stage 13 – Wound30 signature robustness analysis.
13. Figure script – Final publication figure assembly.

## Requirements

- R (version 4.3.0 or higher)
- Key R packages:
  - Matrix
  - edgeR
  - DropletUtils
  - clusterProfiler
  - msigdbr
  - ggplot2, dplyr, tidyr, patchwork, ggrepel, ragg

## Data Availability

All primary datasets analyzed in this study are publicly available from the NCBI Gene Expression Omnibus (GEO) under accession numbers:

- GSE241132
- GSE163973
- GSE181316
- GSE220300
- GSE241124
- GSE265972
- GSE181297

## Citation

If you use this code, please cite:

> [Author names]. [Manuscript title]. *PLOS ONE* (under review).

## Contact

- GitHub: [@Wyh0519](https://github.com/Wyh0519)
