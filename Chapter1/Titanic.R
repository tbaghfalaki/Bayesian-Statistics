# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # for combining plots

# Load data
Data <- read.csv("C:\\Users\\p80744tb\\Desktop\\new_book\\titanic\\train.csv", stringsAsFactors = FALSE)
names(Data)
# Convert columns to factors as appropriate
Data <- Data %>%
  mutate(
    Survived = factor(Survived, labels = c("No", "Yes")),
    Pclass = factor(Pclass),
    Sex = factor(Sex),
    SibSp = factor(SibSp),
    Parch = factor(Parch),
    Embarked = factor(Embarked)
  )

# Select numeric and categorical variables
num_vars <- c("Age", "Fare")
cat_vars <- c("Survived", "Pclass", "Sex", "SibSp", "Parch", "Embarked")

# Prepare data for plotting
df_num_long <- Data %>%
  select(all_of(num_vars)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value")

df_cat_long <- Data %>%
  select(all_of(cat_vars)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value")

# Numeric variable histograms
p_num <- ggplot(df_num_long, aes(x = Value)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 30, alpha = 0.7, na.rm = TRUE) +
  facet_wrap(~ Variable, scales = "free_x") +
  labs(title = "Distribution of Numeric Variables") +
  theme_minimal()

# Categorical variable barplots
p_cat <- ggplot(df_cat_long, aes(x = Value)) +
  geom_bar(fill = "steelblue", na.rm = TRUE) +
  facet_wrap(~ Variable, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Categorical Variables") +
  theme_minimal()

# Combine numeric and categorical plots vertically
p_num / p_cat
