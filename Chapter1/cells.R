# Load required packages
library(ggplot2)
library(dplyr)
library(patchwork)  # for combining plots
library(bellreg)

data(cells)
# Assuming your dataset is named 'cells'
# and contains variables: cells, smoker, weight

# Compute summary statistics
summary_data <- cells %>%
  group_by(weight, smoker) %>%
  summarise(
    mean_cells = mean(cells, na.rm = TRUE),
    prop_nonzero = mean(cells > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Upper panel: Mean of infected cells
p1 <- ggplot(summary_data, aes(x = weight, y = mean_cells, fill = factor(smoker))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    x = "Body Mass Category",
    y = "Mean Number of Infected Cells",
    fill = "Smoking Status"
  ) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "orange"),
                    labels = c("Smoker", "Non-smoker")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")

# Lower panel: Proportion of nonzero counts
p2 <- ggplot(summary_data, aes(x = weight, y = prop_nonzero, fill = factor(smoker))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    x = "Body Mass Category",
    y = "Proportion of Nonzero Counts",
    fill = "Smoking Status"
  ) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "orange"),
                    labels = c("Smoker", "Non-smoker")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

# Combine panels vertically
(p1 / p2) + plot_annotation(tag_levels = 'A')
