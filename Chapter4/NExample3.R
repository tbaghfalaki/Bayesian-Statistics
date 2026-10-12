# --- Load required package and AIS data ---
library(sn)

data(ais)       # load AIS dataset
bmi <- ais$BMI   # extract BMI variable

# --- Compute sample statistics ---
n     <- length(bmi)
x_bar <- mean(bmi)
s     <- sd(bmi)
s2    <- s^2

# --- Hypotheses and prior ---
mu0  <- 23       # Null hypothesis mean
tau2 <- 100      # Prior variance under H1

# --- Compute Bayes Factor (closed-form, plug-in sigma^2 = s^2) ---
# Log-scale computation for numerical stability
logBF10 <- 0.5 * (log(s2) - log(s2 + n * tau2)) +
  (n * tau2 * (x_bar - mu0)^2) / (2 * s2 * (s2 + n * tau2))

BF10 <- exp(logBF10)

# --- Print results ---
cat("Sample size n       =", n, "\n")
cat("Sample mean x_bar   =", x_bar, "\n")
cat("Sample sd s         =", s, "\n")
cat("Sample variance s^2 =", s2, "\n")
cat("Bayes Factor BF10   =", BF10, "\n")
