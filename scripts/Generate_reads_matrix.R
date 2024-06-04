# Code for the analysis of sequencing reads from a Surface-Plex assay
# Including the following steps:
# 1. Load necessary packages
# 2. Create a matrix of read counts from the individual counts files
# 3. Reorder and rename columns based on conditions 

# 1. INITIALIZATION ----------------------------------------------------------

# Set directories for input files and output plots
counts_dir = "path/to/counts/dir"

# load necessary packages
library(tidyverse)
library(magrittr)
library(gtools)
library(sva)
library(purrr)

# 2. Make a matrix of read counts --------------------------------------------

# Get a list of all tab files in the directory
files <- list.files(counts_dir, pattern = "*.txt", full.names = TRUE)

## Load and left-join the tab files into a single tibble

tibble_list <- lapply(files, function(file) {
  df <- read_tsv(file, col_names = c("reads", "Antigen"))
  
  # Get the basename for the file
  file_base_name <- gsub("\\.txt$", "", basename(file))
  
  # Rename the 'reads' column to the basename of the file
  colnames(df)[1] <- file_base_name
  
  return(df)
})

tibble_merged <- purrr::reduce(tibble_list, full_join, by = "Antigen")

tibble_merged <- tibble_merged %>% select(Antigen, everything())

# Reorder the tibble based on natural sorting
new_order <- mixedorder(names(tibble_merged))
tibble_merged <- tibble_merged[, new_order]
tibble_merged <- tibble_merged %>% select(Antigen, everything())

# 3. Reorder and rename columns based on conditions --------------------------
## NOTE: Take care to make sure column names match order of well indices in tibble

Conds = read_tsv ("path/to/Conditions.txt")

new_col_names <- Conds$Conditions
colnames(tibble_merged)[2:4482] <- new_col_names

# save matrix
write.csv(tibble_merged, file = "path/to/saved/raw/matrix.csv")


