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

- `data/`
  - `custom_index_1.fa`: Custom .fasta file for generating the genome index used in the analysis for Assay 1 (see manuscript).
  - `custom_index_2.fa`: Custom .fasta file for generating the genome index used in the analysis for Assay 2 (see manuscript).
  - `Conditions_1.txt`: File containing the conditions and well indices for Assay 1.
  - `Conditions_2a.txt`: File containing the conditions and well indices for Assay 2 - day 1.
  - `Conditions_2b.txt`: File containing the conditions and well indices for Assay 2 - day 2.
  - `Assay1_raw_reads_matrix.csv`: Matrix of raw read counts for Assay 1 (see manuscript for details) - as input for tutorial.

- `scripts/`
  - `Read_pre-processing.sh`: Shell script for preprocessing sequencing reads, including trimming, UMI extraction, alignment, and deduplication.
  - `Generate_Reads_Matrix.R`: R script for creating a read count matrix from individual count files generated by 'Read_pre-processing.sh'
  - `ComBat_batch_correction.R`: R script for batch correcting the read count matrix if samples were processed on different days.
  - `Main_Analysis.R`: R script for CLR normalization, Mann Whitney statistical testing, and results table output with raw or batch corrected counts matrices as input.
    
- `README.md`: This file.

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
  - `magrittr`
  - `gtools`
  - `sva`
  - `purrr`

## TUTORIAL

This tutorial will take you through all the steps of the Surface-plex pipeline analysis. The initial steps require download of raw FastQ data from NCBI Gene Expression Omnibus with the accession number GSE268052. However, you can instead use the pre-generated read matrix for Assay1 from the manuscript, found in data/Assay1_raw_reads_matrix.csv, and **proceed from step 5 of the analysis**, to recreate the results table for this Assay.

### 1. Clone the repository 

  ```bash
git clone https://github.com/ecrosse01/Surface-plex.git
cd Surface-plex
```

### 2. Preprocessing Reads

The preprocessing of sequencing reads from raw FastQ involves generating an index, trimming reads, extracting UMIs, aligning and deduplicating reads.

1. Use script Read_pre-processing.sh - the script is designed to be run in chunks in an interactive session but can be modified to be executed in full.

   For reanalysis of the published data use either 'custom_index_1.fa' or 'custom_index_2.fa' for Assays 1 and 2 from the mansucript, respectively. The raw data can be downloaded from GEO with the accession number GSE268052.

   The output of this pipeline are directories containing individual .txt files with counts per antigen for each sample.

### 3. Analysis in R - Generate a matrix of read counts across all samples and antigens.

This R script generates directly utilizes the count file outputs from Step 2 to create a read count matrix for downstream analysis. Please modify code with correct directory paths.
The conditions and well indices list for Assays 1 and 2 from the manuscript are in files 'data/Conditions_1.txt' and 'data/Conditions_2a/2b.txt' respectively.

1. Ensure all R packages listed in the prerequisites are installed.

2. Run code
  ```r
  source("scripts/Generate_Reads_Matrix.R")
```
### 4. Analysis in R - Batch correction of matrix 

**This step is optional** depending on whether your data was processed across different days / sequenced on multiple lanes etc. Please modify code with correct directory paths.

1. Run code
   ```r
   source("scripts/ComBat_batch_correction.R")
   ```

### 5. Analysis in R - Main analysis - normalization, statistical testing and results output 

A sample matrix for running this main part of the analysis is provided in 'data/Assay1_raw_reads_matrix.csv'. This is a raw read counts matrix from Assay 1 from the manuscript (5x replicates of 84 compounds).

The steps in this code comprise:
1. Centered Log Ratio normalization of the matrix.
2. Mann-Whitney U statistical testing of differences between each compound and the control (DMSO).
3. Calculation of CLR expression differences between compound and the control.
4. Plotting of AnnexinV vs antigen expression variability for each compound (to determine threshold cutoffs for AnnexinV expression and remove compounds with likely highly cytotoxic properties).
   _See "output/Relationship_between_Compound_Variability_and_Cell_Death.png"_
6. Output of results table
   _See "output/Assay1_results.csv"_ 

Run code
   ```r
   source("scripts/Main_Analysis.R")
   ```

### Data Availability
The raw sequencing data used in this analysis have been deposited in the Gene Expression Omnibus (GEO) under accession number GSE268052.

### Citation
If you use these scripts or data in your research, please cite our manuscript:

[TBD]

### License


