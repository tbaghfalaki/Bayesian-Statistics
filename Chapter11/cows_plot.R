rm(list=ls())
library(ggplot2)
library(dplyr)
library(nlme)

# ── Load and prepare data ───────────────────────────────────────────────
data(Milk)
df <- as.data.frame(Milk)

names(df) <- tolower(names(df))

df$diet_lupine   <- ifelse(df$diet == "lupine",        1, 0)
df$diet_barlupin <- ifelse(df$diet == "barley+lupine", 1, 0)
df$cow_id        <- as.integer(factor(df$cow))
df$time_lin      <- df$time
df$time_ch       <- pmax(df$time - 3, 0)
df$time_ch2      <- pmax(df$time - 3, 0)^2

cat("N =", nrow(df), "\n")
cat("Cows =", length(unique(df$cow_id)), "\n")
cat("Weeks =", length(unique(df$time)), "\n")
cat("Diet levels:", levels(df$diet), "\n")
unique(df$diet)

# ── Colour palette ─────────────────────────────────────────────────────
diet_colours <- c("barley"        = "#2166AC",
                  "barley+lupins" = "#F4A582",
                  "lupins"        = "#D6604D")

# ── Mean trajectory ────────────────────────────────────────────────────
mean_traj <- df %>%
  group_by(diet, time) %>%
  summarise(mean_protein = mean(protein), .groups = "drop")

# ── Spaghetti plot ─────────────────────────────────────────────────────
p_spaghetti <- ggplot(df, aes(x = time, y = protein,
                              group = cow, colour = diet)) +
  geom_line(alpha = 0.30, linewidth = 0.4) +
  geom_vline(xintercept = 3,
             linetype   = "dashed",
             colour     = "grey30",
             linewidth  = 0.7) +
  geom_line(data        = mean_traj,
            aes(x       = time,
                y       = mean_protein,
                group   = diet,
                colour  = diet),
            linewidth   = 1.50,
            alpha       = 1.0,
            inherit.aes = FALSE) +
  facet_wrap(~ diet, ncol = 3) +
  scale_colour_manual(values = diet_colours) +
  scale_x_continuous(breaks = seq(1, 19, by = 3)) +
  labs(x = "Week",
       y = "Milk protein content") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(p_spaghetti)
 