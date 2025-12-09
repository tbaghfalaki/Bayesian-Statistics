# Define income variables for the three waves
income_vars <- c("income.1", "income.2", "income.3")

# Set up colors and labels
colors <- c("lightblue", "lightgreen", "lightpink")
wave_labels <- c("Wave 13", "Wave 14", "Wave 15")

# Plot raw income histograms
par(mfrow = c(2,3))  # 1 row, 3 columns for side-by-side plots
for (i in 1:3) {
  hist(Data[[income_vars[i]]],
       main = paste("Income", wave_labels[i]),
       xlab = "Income",
       col = colors[i],
       breaks = 30)
}

# Plot log-transformed income histograms
for (i in 1:3) {
  hist(log(Data[[income_vars[i]]]),
       main = paste("Log Income", wave_labels[i]),
       xlab = "Log(Income)",
       col = colors[i],
       breaks = 30)
}
