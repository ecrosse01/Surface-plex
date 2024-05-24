#########################################################################
#########################################################################
#########################################################################
#########################################################################

# Code for the analysis of sequencing reads from a Surface-Plex assay
# Including the following steps:
# 1. Load necessary packages
# 2. Create a matrix of read counts from the individual counts files
# 3. Reorder and rename columns based on conditions 
# 4. Perform batch correction using ComBat
# 5. CLR normalize the data
# 6. Mann Whitney Testing 
# 7. Plot AnnexinV diff vs. antigen expression variability index 
# 8. Get the final results table

#########################################################################
#########################################################################
#########################################################################
#########################################################################

# Set directories for input files and output plots
counts_dir = "path/to/counts/dir"


# 1. INITIALIZATION ----------------------------------------------------------

# load necessary packages
library(tidyverse)
library(RColorBrewer)
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

# 4. ComBat batch correction -------------------------------------------------

df = tibble_merged

# Remove columns with less than 500 counts

temp_df <- df %>%
  select(-Antigen) %>%
  select(where(~ is.numeric(.) && (sum(., na.rm = TRUE) >= 500)))

# Adding the "Antigen" column back to the filtered tibble
df <- bind_cols(df %>% select(Antigen), temp_df)

# Convert tibble to long format
df_long <- df %>%
  pivot_longer(cols = -Antigen, names_to = "condition", values_to = "expression") %>%
  separate(condition, into = c("plate", "rest_condition"), sep = "_", extra = "merge") %>%
  separate(rest_condition, into = c("replicate_day", "temp_condition"), sep = "_", extra = "merge") %>%
  separate(temp_condition, into = c("condition", "replicate"), sep = "_(?=[^_]+$)") %>%
  unite("replicate_id", replicate_day, replicate, sep = "_", remove = FALSE)

# Create sampleID for batch correction
df_long <- df_long %>%
  mutate(sampleID = paste(replicate_day, condition, replicate, sep = "_"))

# Pivot to wide format
df_wide <- df_long %>%
  select(Antigen, sampleID, expression) %>%
  pivot_wider(names_from = sampleID, values_from = expression)

# Preparing for ComBat
# The first column is Antigen, so we remove it for the expression matrix
expression_matrix <- as.matrix(df_wide[,-1])

# Extract batch information from column names (skip the first column which is 'Antigen')
batch <- gsub("^(R[0-9]+)_.*", "\\1", colnames(df_wide)[-1])

#ComBat normalization
expression_matrix_corrected <- ComBat(dat = expression_matrix, batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)

# Perform PCA on the original and corrected data
pca_original <- prcomp(t(expression_matrix), center = TRUE, scale. = TRUE)
pca_corrected <- prcomp(t(expression_matrix_corrected), center = TRUE, scale. = TRUE)

# Prepare data frames for ggplot
df_original <- data.frame(PC1 = pca_original$x[,1], PC2 = pca_original$x[,2], Batch = batch)
df_corrected <- data.frame(PC1 = pca_corrected$x[,1], PC2 = pca_corrected$x[,2], Batch = batch)

# Plotting
p1 <- ggplot(df_original, aes(x = PC1, y = PC2, color = Batch)) + 
  geom_point() + 
  ggtitle("Before Batch Correction") +
  theme_minimal()

p2 <- ggplot(df_corrected, aes(x = PC1, y = PC2, color = Batch)) + 
  geom_point() + 
  ggtitle("After Batch Correction") +
  theme_minimal()

df_corrected = as_tibble(expression_matrix_corrected)

# Add the Antigen column from the original df_wide
df_corrected <- tibble(Antigen = df_wide$Antigen, df_corrected)

# Create long verion of batch corrected tibble
df_corrected_long <- df_corrected %>%
  pivot_longer(cols = -Antigen, names_to = "condition", values_to = "expression") %>%
  separate(condition, into = c("replicate_day", "temp_condition"), sep = "_", extra = "merge") %>%
  separate(temp_condition, into = c("condition", "replicate"), sep = "_(?=[^_]+$)") %>%
  unite("replicate_id", replicate_day, replicate, sep = "_", remove = FALSE)



# 5. CLR normalize the data -----------------------------------------------

df_normalized <- df_corrected_long %>%
  # Group by Sample
  group_by(replicate_day, condition, replicate) %>%
  # Calculate the total counts per sample
  mutate(TotalCount = sum(expression)) %>%
  # Convert counts to proportions within each sample
  mutate(Proportion = expression / TotalCount) %>%
  # Log transformation of proportions (adding a small constant to avoid log(0))
  mutate(LogProportion = log(Proportion + 1e-6)) %>%
  # Centering: Subtract the average log value of all antigens in a sample
  mutate(CLR = LogProportion - mean(LogProportion, na.rm = TRUE)) 


# 6. Mann Whitney Testing ----------------------------------------------------

# Define control data
control_data <- df_normalized %>%
  filter(condition == "CTRL")

# Define mann_whitney test function
perform_mann_whitney <- function(antigen, test_data, control_data) {
  test_values <- test_data$CLR
  control_values <- control_data %>% 
    filter(Antigen == antigen) %>% 
    pull(CLR)
  
  # Check if there are enough non-missing values in both test and control groups
  if(length(na.omit(test_values)) < 3 || length(na.omit(control_values)) < 3) {
    return(tibble(statistic = NA, p.value = NA, method = "Insufficient data"))
  }
  
  test_result <- wilcox.test(test_values, control_values, exact = FALSE)
  broom::tidy(test_result)
}

# Run mann-whitney U test
results <- df_normalized %>%
  filter(condition != "CTRL") %>%
  group_by(Antigen, condition) %>%
  nest() %>%
  rowwise() %>%
  mutate(test_results = list(perform_mann_whitney(Antigen, data, control_data))) %>%
  unnest(cols = c(test_results))

# Calculate the mean for each compound/antigen combo and calculate the difference between that and control

Ave = df_normalized %>%
  group_by(Antigen, condition) %>%
  summarise(mean_expression = mean(CLR, na.rm = TRUE))

condition_c_avg <- Ave %>%
  filter(condition == "CTRL") 

# Calculate Diff
Diff <- Ave %>%
  left_join(condition_c_avg, by = c("Antigen")) %>%
  mutate(Diff = mean_expression.x - mean_expression.y) %>%
  rename(condition = "condition.x") %>%
  select(Antigen, condition, Diff)

results = results %>%
  left_join(Diff, by = c("Antigen", "condition"))



# 7. Plot AnnexinV diff vs. antigen expression variability index -------

# Convert Diff.x to absolute values in the original dataset
results$Abs_Diff_x <- abs(results$Diff)

# Calculate the variability index (mean of absolute Diff.x) for each compound
variability_index <- results %>%
  group_by(condition) %>%
  summarise(Variability_Index = mean(Abs_Diff_x))

# Calculate the mean AnnexinV Diff.x for each condition
annexinV_agg <- results %>%
  filter(Antigen == "AnnexinV") %>%
  group_by(condition) %>%
  summarise(AnnexinV_avg_diff = mean(Diff))

# Merge the variability index with the AnnexinV data
merged_data <- left_join(variability_index, annexinV_agg, by = "condition")

# Create the plot
ggplot(merged_data, aes(x = Variability_Index, y = AnnexinV_avg_diff)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Variability Index (Mean of Absolute Diff.x)",
       y = "Average AnnexinV Diff.x",
       title = "Relationship between Compound Variability and Cell Death")


# 8. FINAL RESULTS OUTPUT TABLE ----------------------------------------------

results = results %>% left_join (annexinV_agg, by = "condition") %>%
  select(Antigen, condition, statistic, p.value, method, Diff, AnnexinV_avg_diff)








