rm(list=ls())
# Target PMF
target_pmf <- c(0.15, 0.12, 0.09, 0.08, 0.12, 0.14, 0.09, 0.07, 0.04, 0.10)
sum(target_pmf)
j_vals <- 1:10
q_val <- 0.1
k <- max(target_pmf / q_val)

# Sampling function using acceptance-rejection
acDiscrete <- function(n) {
  y <- c()
  while (length(y) < n) {
    candidate <- sample(1:10, 1)
    u <- runif(1)
    if (u <= target_pmf[candidate] / (k * q_val)) {
      y <- c(y, candidate)
    }
  }
  return(y)
}

# Generate samples
set.seed(123)  # for reproducibility
n <- 10000
samples <- acDiscrete(n)

# Empirical PMF
empirical_pmf <- table(factor(samples, levels = 1:10)) / n

# Plot both PMFs
barplot(
  rbind(target_pmf, empirical_pmf),
  beside = TRUE,
  col = c("skyblue", "salmon"),
  names.arg = j_vals,
  ylim = c(0, 0.14),
  legend.text = c("Target PMF", "Sampled PMF"),
  args.legend = list(x = "topright", bty = "n"),
  main = "Comparison of Target and Sampled PMFs",
  ylab = "Probability",
  xlab = "j"
)
