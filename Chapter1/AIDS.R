library(ggplot2)
library(JM)

data(aids)

# Plot CD4 trajectories by drug with one color per panel
ggplot(aids, aes(x = obstime, y = CD4, group = patient)) +
  geom_line(alpha = 0.5, aes(color = drug)) +  # map color to drug
  facet_wrap(~ drug) +
  scale_x_continuous(breaks = c(0, 2, 6, 12, 18)) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E")) +  # assign one color per drug level
  labs(title = "CD4 Count Trajectories Over Time by Drug",
       x = "Observation Time (weeks)",
       y = "CD4 Count",
       color = "Drug") +
  theme_minimal() +
  theme(legend.position = "none")  # hide legend if desired




library(survminer)
splots <- list()
Bfit1 <- survfit(Surv(Time, death) ~ drug, data = aids.id) 
ggsurvplot(Bfit1, conf.int = TRUE,
           risk.table = FALSE, risk.table.col = "strata",
           ggtheme = theme_bw() , legend.title = element_blank(), 
           legend.labs = c("ddC", "ddI"),
           palette = c("#1F77B4", "#FF7F0E") ) 
