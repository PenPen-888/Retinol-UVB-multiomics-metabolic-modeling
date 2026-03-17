# Retinol-UVB Metabolic Reprogramming

This repository contains the code and processed data used in the study:

**"Cell-type-specific metabolic reprogramming underlies retinol’s efficacy against UVB-induced photoaging"**

## Overview

This study investigates cell-type-specific metabolic responses to retinol treatment under UVB stress using a multi-layered systems biology framework, including:

- Single-omics analysis (transcriptomics, proteomics, metabolomics)
- Cross-omics integration
- Multi-omics integration using DIABLO (mixOmics)
- Genome-scale metabolic modeling (ftINIT, RIPTiDe)
- Flux analysis and machine learning (random forest, Cliff’s delta) to identify key reactions

## Repository Structure
data/ # processed datasets used in the study
scripts/ # single-omics and cross-omics analysis scripts
model/ # genome-scale metabolic modeling code and results
figures/ # visualization scripts

## Requirements

The analyses were performed using the following environments:

- R (≥ 4.2)
- Python (≥ 3.9)
- MATLAB (R2024b)
- Gurobi Optimizer (compatible with MATLAB R2024b)
