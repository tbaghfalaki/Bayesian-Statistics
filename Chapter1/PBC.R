rm(list = ls())
library(JM) # for data pbc2.id
library(cmprsk) # for cuminc()
library(survminer) # for ggcompetingrisks()
library(cowplot) # for theme_cowplot()
library(ggsci) # for scale_fill_jco()
library(ggplot2)
library(dplyr)
library(tidyr)

data("pbc2")
data("pbc2.id")

library(ggplot2)

# Serum Bilirubin
ggplot(pbc2, aes(x = year, y = serBilir, group = id)) +
  geom_line(alpha = 0.3) +
  labs(title = "Serum Bilirubin",
       x = "Years",
       y = "Serum Bilirubin") +
  theme_minimal()

# Albumin
ggplot(pbc2, aes(x = year, y = albumin, group = id)) +
  geom_line(alpha = 0.3) +
  labs(title = "Albumin",
       x = "Years",
       y = "Albumin") +
  theme_minimal()

# Serum Cholesterol
ggplot(pbc2, aes(x = year, y = serChol, group = id)) +
  geom_line(alpha = 0.3) +
  labs(title = "Serum Cholesterol",
       x = "Years",
       y = "Serum Cholesterol") +
  theme_minimal()

# SGOT
ggplot(pbc2, aes(x = year, y = SGOT, group = id)) +
  geom_line(alpha = 0.3) +
  labs(title = "SGOT",
       x = "Years",
       y = "SGOT") +
  theme_minimal()

# Platelets
ggplot(pbc2, aes(x = year, y = platelets, group = id)) +
  geom_line(alpha = 0.3) +
  labs(title = "Platelets",
       x = "Years",
       y = "Platelets") +
  theme_minimal()









# Competing risks analysis


data("pbc2.id")
names(pbc2.id)
attach(pbc2.id)

pbc2.id$status1 <- as.numeric(pbc2.id$status) - 1
fit <- cmprsk::cuminc(
  ftime = pbc2.id$years, fstatus = pbc2.id$status,
  cencode = "alive", group = ""
)

ggcompetingrisks(fit,
                 conf.int = TRUE, multiple_panels = FALSE,
                 palette = c("#E95420", "#3CB371"), legend.labs = c("Dead", "Transplant"),
                 subset = "", legend.title = ""
) + ggtitle("") + xlab("Time (year)") +
  ylab("Cumulative incidence of event") + theme_cowplot() + scale_fill_jco()



ggcompetingrisks(fit,
                 conf.int = TRUE,
                 legend.title = "", cencode = "alive"
) +
  xlab("Time (year)") +
  theme_cowplot() + scale_fill_jco() + theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.direction = "vertical"
  ) +
  theme(
    rect = element_rect(fill = "gray80", colour = "gray80"),
    panel.background = element_rect(
      fill = "gray80",
      colour = NA
    ),
    plot.background = element_rect(colour = "gray80")
  )
