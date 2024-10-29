## ----------------------- Question 1 ------------------------------------------
# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Read the data
data <- read.csv("data/biomarker-raw.csv")

# Select a subset of columns for proteins (skip the first two columns)
protein_columns <- colnames(data)[3:length(colnames(data))]

# Randomly select 5 protein columns
set.seed(123) # Set seed for reproducibility
selected_proteins <- sample(protein_columns, 5)

# Convert selected columns to numeric and handle any issues with non-numeric data
data[selected_proteins] <- lapply(data[selected_proteins], function(x) as.numeric(as.character(x)))

# Filter to keep only rows without NA values in selected proteins
protein_data <- data %>% select(all_of(selected_proteins)) %>% drop_na()

# Reshape data for ggplot
protein_data_long <- protein_data %>% 
  pivot_longer(cols = everything(), names_to = "Protein", values_to = "Raw_Value") %>%
  mutate(Raw_Value = as.numeric(Raw_Value))  # Ensure Raw_Value is numeric

# Plot distributions using ggplot2
ggplot(protein_data_long, aes(x = Raw_Value)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  facet_wrap(~ Protein, scales = "free_x") +
  labs(title = "Distribution of Raw Values for Selected Proteins",
       x = "Raw Value", y = "Frequency") +
  theme_minimal()

# Cleaned data distribution after being log transformed

# Load required libraries
library(tidyverse)

# Load the cleaned biomarker data
load('data/biomarker-clean.RData')

sampled_proteins <- biomarker_clean %>%
  select(sample(names(.), 5)) %>%
  pivot_longer(cols = everything(), names_to = "Protein", values_to = "Transformed_Value")

# Plot histograms for each selected protein using facet_wrap (Transformed Values)
ggplot(sampled_proteins, aes(x = Transformed_Value)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_wrap(~ Protein, scales = "free") +
  labs(title = "Histograms of Transformed Values for Selected Proteins",
       x = "Transformed Value",
       y = "Frequency") +
  theme_minimal()

## ----------------------- Question 2 ------------------------------------------

## ----------------------- Question 3 ------------------------------------------

## ----------------------- Question 4 ------------------------------------------


