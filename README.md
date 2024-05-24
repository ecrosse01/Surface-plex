# Surface-plex

## Surface-Plex Assay Analysis

Surface-Plex is a high-throughput assay that quanitfies the surface expression of therapeutic proteins of interest following treatment with thousands of compounds.
This repository contains the scripts and data files used for the analysis of sequencing reads from a Surface-Plex assay. 
The analysis includes read preprocessing, UMI extraction, alignment, batch correction, centered log ratio (CLR) normalization and differential expression testing.

## Relevant Paper and Abstract

**Title:** High-throughput quantification of cell-surface proteins for accelerated cancer therapeutic discovery

**Abstract:**
Cell-surface proteins on tumor cells serve as both key targets for immunotherapeutic strategies and markers of malignant cell phenotype and differentiation status. Pharmacologic agents that modulate surface expression of immunotherapeutic targets or induce terminal differentiation of cancer initiating cells therefore offer an attractive strategy for more robust and precise targeting of cancer. Both conventional and newer methods for characterizing cell surface protein expression on tumor cells, including flow cytometry and CITE-Seq, are invaluable for characterizing phenotypes but do not readily scale for high-throughput drug screening. Here we present Surface-plex, a technology which enables multiplexed, high-throughput sequencing readout of cell-surface proteins on tumor cells following treatment with thousands of compounds. We applied this technology to acute myeloid leukemia and identified potent chemical inducers of the therapeutically targetable antigen CD47 and several compounds as novel druggable promoters of blast differentiation. Surface-plex thus offers a scalable platform enabling rapid identification of compounds with therapeutic potential, applicable across cancer types, to accelerate targeted cancer therapeutic development.

**Link to Paper:** [URL to the paper once it is published]

## Repository Contents

- `Read_pre-processing.sh`: Shell script for preprocessing sequencing reads, including trimming, UMI extraction, alignment, and deduplication.
- `Surface-Plex_analysis.R`: R script for creating a read count matrix, performing batch correction, and generating plots.
- `custom_index.fa`: Custom .fasta file for generating the genome index used in the analysis.
- `Conditions.txt`: File containing the conditions and well indices for the experiment.

## Prerequisites

### Software Requirements

- **Shell Script**:
  - `STAR` (>= 2.7.9a)
  - `cutadapt` (>= 4.1)
  - `UMI-tools` (>= 1.0.1)
  - `SAMtools` (>= 1.16.1)

- **R Script**:
  - R (>= 4.0.0)
  - `tidyverse`
  - `RColorBrewer`
  - `magrittr`
  - `gtools`
  - `sva`
  - `purrr`

## Instructions

### 1. Preprocessing Reads

The preprocessing of sequencing reads involves generating an index, trimming reads, extracting UMIs, aligning and deduplicating reads.

a. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/Surface-plex.git
   cd Surface-plex

b. Use script Read_pre-processing.sh - the script is designed to be run in chunks in an interactive session but can be modified to be executed in full.

### 2. Analysis in R

  
