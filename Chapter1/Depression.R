# Load required packages
library(ggplot2)

Data <- read.csv("C:\\Users\\p80744tb\\Desktop\\Bayesian Statistics\\Rcodes\\Data\\student_depression_dataset.csv", stringsAsFactors = FALSE)
attach(Data)
#names(Data)



# 1. Gender (categorical)
ggplot(Data, aes(x = Gender)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Distribution of Gender", x = "Gender", y = "Count") +
  theme_minimal()

# 2. Age (numeric)
ggplot(Data, aes(x = Age)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.6) +
  labs(title = "Age Distribution with Density Curve", x = "Age", y = "Density") +
  theme_minimal()

# 3. City (categorical)
ggplot(Data, aes(x = City)) +
  geom_bar(fill = "coral") +
  labs(title = "Distribution of City", x = "City", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 4. Profession (categorical)
ggplot(Data, aes(x = Profession)) +
  geom_bar(fill = "darkseagreen") +
  labs(title = "Distribution of Profession", x = "Profession", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. Academic Pressure (numeric/ordinal)
ggplot(Data, aes(x = Academic.Pressure)) +
  geom_histogram(binwidth = 1, fill = "orchid", color = "black", alpha = 0.7) +
  labs(title = "Academic Pressure Distribution", x = "Academic Pressure", y = "Count") +
  theme_minimal()

# 6. Work Pressure (numeric/ordinal)
ggplot(Data, aes(x = Work.Pressure)) +
  geom_histogram(binwidth = 1, fill = "plum", color = "black", alpha = 0.7) +
  labs(title = "Work Pressure Distribution", x = "Work Pressure", y = "Count") +
  theme_minimal()

# 7. CGPA (numeric)
ggplot(Data, aes(x = CGPA)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgreen", color = "black", alpha = 0.6) +
  #geom_density(color = "darkgreen", size = 1) +
  labs(title = "CGPA Distribution", x = "CGPA", y = "Density") +
  theme_minimal()

# 8. Study Satisfaction (numeric/ordinal)
ggplot(Data, aes(x = Study.Satisfaction)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Study Satisfaction Distribution", x = "Study Satisfaction", y = "Count") +
  theme_minimal()

# 9. Job Satisfaction (numeric/ordinal)
ggplot(Data, aes(x = Job.Satisfaction)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Job Satisfaction Distribution", x = "Job Satisfaction", y = "Count") +
  theme_minimal()

# 10. Sleep Duration (numeric)
# Create a frequency table
sleep_table <- as.data.frame(table(Data$Sleep.Duration))

# Rename columns for clarity
colnames(sleep_table) <- c("Sleep.Duration", "Count")

# Plot bar chart
ggplot(sleep_table, aes(x = Sleep.Duration, y = Count)) +
  geom_bar(stat = "identity", fill = "mediumpurple") +
  labs(title = "Frequency of Sleep Duration", x = "Sleep Duration (hours)", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 11. Dietary Habits (categorical)
ggplot(Data, aes(x = Dietary.Habits)) +
  geom_bar(fill = "goldenrod") +
  labs(title = "Distribution of Dietary Habits", x = "Dietary Habits", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 12. Degree (categorical)
ggplot(Data, aes(x = Degree)) +
  geom_bar(fill = "deepskyblue") +
  labs(title = "Distribution of Degree", x = "Degree", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 13. Suicidal Thoughts (binary)
ggplot(Data, aes(x = Have.you.ever.had.suicidal.thoughts..)) +
  geom_bar(fill = "tomato") +
  labs(title = "Suicidal Thoughts (Yes/No)", x = "Suicidal Thoughts", y = "Count") +
  theme_minimal()

# 14. Work Study Hours (numeric)
workstudy_table <- as.data.frame(table(Data$Work.Study.Hours))
barplot(workstudy_table$Freq, names.arg = workstudy_table$Var1,
        las = 2, col = "steelblue",
        main = "Frequency of Work and Study Hours",
        xlab = "Work and Study Hours",
        ylab = "Count")


# 15. Financial Stress (numeric/ordinal)
# Remove rows where Financial.Stress is "?" or NA
clean_data <- Data[Data$Financial.Stress != "?" & !is.na(Data$Financial.Stress), ]

# Create frequency table without "?"
Financial.Stress_table <- as.data.frame(table(clean_data$Financial.Stress))

# Plot barplot
barplot(Financial.Stress_table$Freq, names.arg = Financial.Stress_table$Var1,
        las = 2, col = "steelblue",
        main = "Frequency of Financial Stress Levels",
        xlab = "Financial Stress",
        ylab = "Count")



# 16. Family History of Mental Illness (binary)
ggplot(Data, aes(x = Family.History.of.Mental.Illness)) +
  geom_bar(fill = "seagreen") +
  labs(title = "Family History of Mental Illness", x = "History", y = "Count") +
  theme_minimal()

# 17. Depression (numeric score)
Depression_table <- as.data.frame(table(Data$Depression))

barplot(Depression_table$Freq, names.arg = Depression_table$Var1,
        las = 2, col = "goldenrod",
        main = "Frequency of Depression Status",
        xlab = "Depression",
        ylab = "Count")



# 13. Suicidal Thoughts (binary)
ggplot(Data, aes(x = Have.you.ever.had.suicidal.thoughts..)) +
  geom_bar(fill = "tomato") +
  labs(title = "Suicidal Thoughts (Yes/No)", x = "Suicidal Thoughts", y = "Count") +
  theme_minimal()
# Create frequency table for suicidal thoughts variable
suicide_table <- as.data.frame(table(Data$Have.you.ever.had.suicidal.thoughts..))

# Plot barplot
barplot(suicide_table$Freq, names.arg = suicide_table$Var1,
        las = 2, col = "tomato",
        main = "Frequency of Suicidal Thoughts",
        xlab = "Suicidal Thoughts (Yes/No)",
        ylab = "Count")
