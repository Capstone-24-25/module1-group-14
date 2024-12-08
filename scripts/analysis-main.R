## ----------------------- Question 1 ------------------------------------------
# Load necessary libraries
library(tidyverse)
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
# exploratory analysis of outlying values

library(tidyverse)
library(ggplot2)

# function to find outliers using IQR
find_outliers <- function(x) {
  # IQR and upper/lower bounds
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # returns TRUE if outlier
  (x < lower_bound | x > upper_bound)
}

# testing function on data
outliers <- biomarker_clean2 %>%
  mutate(across(.cols = -c(group, ados), ~ find_outliers(.x))) %>%
  summarise(across(.cols = -c(group, ados), sum, na.rm = TRUE))

print(outliers)

# applying outlier function
# not 'group' and 'ados' : categorical
outlier_total <- biomarker_clean2 %>%
  mutate(across(.cols = -c(group, ados), ~ find_outliers(.x))) %>%
  rowwise() %>%
  mutate(outlier_count = sum(c_across(-c(group, ados)))) %>%
  ungroup() %>%
  select(group, ados, outlier_count)

# subjects with the most outliers
top_outliers <- outlier_total %>%
  arrange(desc(outlier_count))

print(top_outliers)

# which group of subjects has higher outlier counts
group_outliers <- outlier_total %>%
  group_by(group) %>%
  summarise(mean_outliers = mean(outlier_count),
            max_outliers = max(outlier_count),
            min_outliers = min(outlier_count),
            count = n())

print(group_outliers)


# boxplot to visualizing outlier count by subject (TD or ASD)
ggplot(outlier_total, aes(x = group, y = outlier_count)) +
  geom_boxplot() +
  labs(title = "Outlier Counts by Group", 
       x = "Group", 
       y = "Outlier Count") +
  theme_minimal()


# histogram for outlier counts by group, separated by color
ggplot(outlier_total, aes(x = outlier_count, fill = group)) +
  geom_histogram(binwidth = 5, position = "dodge", alpha = 0.7) +
  labs(title = "Histogram of Outlier Counts by Group",
       x = "Outlier Count",
       y = "Frequency") +
  theme_minimal() +
  scale_fill_manual(values = c("ASD" = "blue", "TD" = "orange"))  
## ----------------------- Question 3 ------------------------------------------

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

# Modification 1: Training Partition
set.seed(27)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8)
biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_train %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_train %>%
  select(-c(group, ados))

response <- biomarker_train %>% pull(group) %>% factor()

# fit RF
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_test <- biomarker_test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model
fit <- glm(class ~ ., 
           data = biomarker_sstar, 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

mod1_results <- biomarker_sstar_test %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

# Modification 2: Greater Than 10 Predictive Proteins (20 proteins)

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

# Modification 3: Fuzzy intersection to combine selected proteins

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- union(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

## ----------------------- Question 4 ------------------------------------------

data <- read.csv("biomarker-raw.csv")
var_names <- read_csv('biomarker-raw.csv', 
                      col_names = F,
                      show_col_types = F,
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}
biomarker_clean <- read_csv('biomarker-raw.csv', 
                            skip = 2,
                            show_col_types = F,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(.cols = -c(group, ados), 
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  select(group, ados, everything())








# Get all protein column names (excluding 'group' and 'ados' columns)
all_proteins <- biomarker_clean %>%
  select(-c(group, ados)) %>%
  colnames()

# Create all possible subsets of protein columns for grid search
protein_combinations <- list(all_proteins) # Start with all proteins as one set

# Modify the grid search function to search on all combinations of proteins
grid_search_auc <- function(biomarker_data, protein_combinations, response_col, n_trees = 1000) {
  results <- data.frame(Panel = character(), AUC = double(), stringsAsFactors = FALSE)
  
  for (combo in protein_combinations) {
    # Predictors
    predictors <- biomarker_data %>%
      select(all_of(combo))
    
    # Response variable
    response <- biomarker_data %>% pull(!!sym(response_col)) %>% factor()
    
    # Train/test split for cross-validation
    set.seed(123)
    train_index <- createDataPartition(response, p = 0.8, list = FALSE)
    train_data <- predictors[train_index, ]
    train_response <- response[train_index]
    test_data <- predictors[-train_index, ]
    test_response <- response[-train_index]
    
    # Random Forest model
    rf_model <- randomForest(x = train_data,
                             y = train_response,
                             ntree = n_trees,
                             importance = TRUE)
    
    # Calculate AUC
    test_pred <- predict(rf_model, test_data, type = "prob")[, 2]
    auc <- roc(test_response, test_pred, levels = rev(levels(test_response)))$auc
    
    results <- rbind(results, data.frame(Panel = paste(combo, collapse = ", "), AUC = auc))
  }
  
  best_panel <- results %>% arrange(desc(AUC)) %>% slice(1)
  return(best_panel)
}

# Run the grid search on all proteins
best_protein_panel <- grid_search_auc(biomarker_clean, list(all_proteins), response_col = "group", n_trees = 1000)
print(best_protein_panel)

## Lets see if we can do better





# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)


protein_panels <- list(proteins_s1, proteins_s2, intersect(proteins_s1, proteins_s2))

best_protein_panel <- grid_search_auc(biomarker_clean, protein_panels, response_col = "group", n_trees = 1000)

print(best_protein_panel)
