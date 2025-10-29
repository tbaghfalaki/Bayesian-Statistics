# Load libraries
library(R2jags)
library(bellreg)
library(dplyr)

# Load data
data(cells)

# Prepare dataset
cells <- cells %>%
  mutate(
    smoker = as.numeric(smoker),
    gender = as.numeric(gender),
    weightf = factor(weight, levels = c("normal", "over", "obese")),
    agef = factor(age, levels = c("young", "mid", "old"))
  )

# Model matrix: Intercept, Smoking, BMI(Overweight/Obese), Age(Middle/Old), Gender
X <- model.matrix(~ smoker + weightf + agef + gender, data = cells)
K <- ncol(X)
N <- nrow(X)
Y <- cells$cells   # response variable (count of infected cells)

# Tabulate counts of Y
y_counts <- table(Y)

# Tabulate counts of Y
y_counts <- table(Y)

# Basic improved barplot
barplot(
  y_counts,
  main = "Distribution of Infected Cells",
  xlab = "Number of Infected Cells",
  ylab = "Frequency",
  col = "skyblue",
  border = "white"
)
vcdExtra::zero.test(Y)
 