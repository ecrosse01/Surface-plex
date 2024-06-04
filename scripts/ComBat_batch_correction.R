# Code for the analysis of sequencing reads from a Surface-Plex assay
# Including the following steps:
# 4. Perform batch correction using ComBat


# load necessary packages
library(tidyverse)
library(magrittr)
library(gtools)
library(sva)
library(purrr)


# 4. ComBat batch correction -------------------------------------------------

df = read_csv(file = "path/to/raw/reads/matrix.csv")

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

# Save .csv
write.csv(df_corrected, file = "Batch_corrected_matrix.csv")



