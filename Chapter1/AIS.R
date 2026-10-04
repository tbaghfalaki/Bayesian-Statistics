rm(list=ls())
# Load required packages
library(sn)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load AIS data
data(ais)


# Categorical variables (sex and sport)
ggplot(ais, aes(x = sex)) +
  geom_bar(fill = "mediumseagreen") +
  labs(title = "Distribution of Sex", x = "Sex", y = "Count") +
  theme_minimal()


dim(ais)

table(ais$sex)

# Sport name mapping (adjust names as needed)
sport_names <- c(
  BB = "Basketball",
  Gym = "Gymnastics",
  Net = "Netball",
  Row = "Rowing",
  Swim = "Swimming",
  T400m = "Track 400m",
  T800m = "Track 800m",
  Tenn = "Tennis",
  Wp = "Water Polo"
)

ais %>%
  mutate(sport_full = recode(sport, !!!sport_names)) %>%
  count(sport_full) %>%
  ggplot(aes(x = sport_full, y = n)) +
  geom_col(fill = "mediumpurple", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 4) +
  labs(title = "Distribution of Sport",
       x = "Sport",
       y = "Count") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Numeric variables — put them in long format for faceting
num_vars <- c("BMI", "RCC", "WCC", "Hc", "Hg", "SSF", "Bfat", "LBM", "Ht", "Wt")

ais_long <- ais %>%
  select(all_of(num_vars)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

# Histogram with density overlay
ggplot(ais_long, aes(x = Value)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "skyblue", alpha = 0.6) +
  geom_density(color = "indianred1", size = .51) +
  facet_wrap(~Variable, scales = "free", ncol = 4) +
  labs(title = "Histograms with Density Curves of AIS Numeric Variables",
       x = "Value", y = "Density") +
  theme_minimal(base_size = 14)
