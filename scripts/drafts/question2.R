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


# boxplot to visualizing outlier count by group
ggplot(outlier_total, aes(x = group, y = outlier_count)) +
  geom_boxplot() +
  labs(title = "Outlier Counts by Group", 
       x = "Group", 
       y = "Outlier Count") +
  theme_minimal()


# histogram for outlier counts by group
ggplot(outlier_total, aes(x = outlier_count, fill = group)) +
  geom_histogram(binwidth = 5, position = "dodge", alpha = 0.7) +
  labs(title = "Histogram of Outlier Counts by Group",
       x = "Outlier Count",
       y = "Frequency") +
  theme_minimal() +
  scale_fill_manual(values = c("ASD" = "blue", "TD" = "orange"))  

'''
Both groups have similar mean outliers and total outliers.
There does not appear to be more outliers in group than the other, 
when looking at group_outliers.
As shown by the boxplot we can see that there is more variation in the TD group
than the ASD group when it comes to outliers.
'''


