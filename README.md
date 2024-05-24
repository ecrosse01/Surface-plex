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
- `Surface-Plex_analysis.R`: R script for creating a read count matrix, performing batch correction, centered log ratio (CLR) normalization and differential expression testing.
- `custom_index_1.fa`: Custom .fasta file for generating the genome index used in the analysis for Assay 1 (see manuscript).
- `custom_index_2.fa`: Custom .fasta file for generating the genome index used in the analysis for Assay 2 (see manuscript).
- `Conditions_1.txt`: File containing the conditions and well indices for Assay 1.
- `Conditions_2.txt`: File containing the conditions and well indices for Assay 2.

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

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/Surface-plex.git
   cd Surface-plex

2. Use script Read_pre-processing.sh - the script is designed to be run in chunks in an interactive session but can be modified to be executed in full.

   For reanalysis of the published data use either 'custom_index_1.fa' or 'custom_index_2.fa' for Assays 1 and 2 from the mansucript, respectively.

### 2. Analysis in R

The R script generates a read count matrix, performs batch correction, centered log ratio (CLR) normalization and differential expression testing.

1. Ensure all R packages listed in the prerequisites are installed.

2. Run Surface-Plex_analysis.R in chunks. Code must be adapted for your specific experimental set up.

The conditions and well indices list for Assays 1 and 2 from the manuscript are in files 'Conditions_1.txt' and 'Conditions_2.txt' respectively.

### Data Availability
The raw sequencing data used in this analysis have been deposited in the Gene Expression Omnibus (GEO) under accession number GSE268052.

### Citation
If you use these scripts or data in your research, please cite our manuscript:

[TBD]

### License


