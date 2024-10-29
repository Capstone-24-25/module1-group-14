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

load('data/biomarker-clean.RData')

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






# function to perform grid search on protein panels
grid_search_auc <- function(biomarker_data, protein_panels, response_col, n_trees = 1000) {
  results <- data.frame(Panel = character(), AUC = double(), stringsAsFactors = FALSE)
  
  for (panel in protein_panels) {
    #predictors
    predictors <- biomarker_data %>%
      select(all_of(panel))
    
    #response variable 
    response <- biomarker_data %>% pull(!!sym(response_col)) %>% factor()
    
    #train/test split for cross-validation
    set.seed(123)
    train_index <- createDataPartition(response, p = 0.8, list = FALSE)
    train_data <- predictors[train_index, ]
    train_response <- response[train_index]
    test_data <- predictors[-train_index, ]
    test_response <- response[-train_index]
    
    #Random Forest model
    rf_model <- randomForest(x = train_data, 
                             y = train_response, 
                             ntree = n_trees, 
                             importance = TRUE)
    
    test_pred <- predict(rf_model, test_data, type = "prob")[, 2]
    auc <- roc(test_response, test_pred, levels = rev(levels(test_response)))$auc
    
    results <- rbind(results, data.frame(Panel = paste(panel, collapse = ", "), AUC = auc))
  }
  
  best_panel <- results %>% arrange(desc(AUC)) %>% slice(1)
  return(best_panel)
}

protein_panels <- list(proteins_s1, proteins_s2, intersect(proteins_s1, proteins_s2))

best_protein_panel <- grid_search_auc(biomarker_clean, protein_panels, response_col = "group", n_trees = 1000)

print(best_protein_panel)