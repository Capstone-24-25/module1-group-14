# temporarily removing the outlier trimming from preprocessing.R

library(tidyverse)


# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean2 <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  
  
  
  # log transform, center and scale, and trim
  #mutate(across(.cols = -c(group, ados), 
  #~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  
  # log transform, center and scale, without trimming
  # dont need .at=3 bc not part of scale() parameters
  mutate(across(.cols = -c(group, ados),
                ~ scale(log10(.x))[,1])) %>%
  
  
  
  # reorder columns
  select(group, ados, everything())

# export as r binary
save(list = 'biomarker_clean2', 
     file = 'data/biomarker-clean2.RData')