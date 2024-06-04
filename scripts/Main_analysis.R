# Code for the analysis of sequencing reads from a Surface-Plex assay
# Including the following steps:
# 5. CLR normalize the data
# 6. Mann Whitney Testing 
# 7. Plot AnnexinV diff vs. antigen expression variability index 
# 8. Get the final results table

# load necessary packages
library(tidyverse)
library(magrittr)
library(gtools)
library(sva)
library(purrr)


# load raw or batch corrected reads matrix

df <- read_csv(file = "data/Assay1_raw_reads_matrix.csv")
df[,1] = NULL

# 5. CLR normalize the data -----------------------------------------------

# Convert tibble to long format
df_long <- df %>%
  pivot_longer(cols = -Antigen, names_to = "condition", values_to = "expression") %>%
  separate(condition, into = c("condition", "replicate"), sep = "_(?=[^_]+$)") %>% # make separate condition and replicate columns
  mutate(condition = if_else(condition == "C", "CTRL", condition)) # rename C to CTRL

# CLR normalization
df_normalized <- df_long %>%
  # Group by Sample
  group_by(condition, replicate) %>%
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


write.csv(results, file = "Assay1_results.csv")


