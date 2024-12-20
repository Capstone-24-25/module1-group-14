library(caret)
library(pROC)
library(randomForest)
library(dplyr)
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

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
