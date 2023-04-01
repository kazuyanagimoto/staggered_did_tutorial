library(tidyverse)
library(here)
library(bacondecomp)
library(fixest)

load(here("output/r/simulation/1_gen_data/data.rds"))

## Diagnosis
# Goodman-Bacon Decomposition (https://github.com/evanjflack/bacondecomp.git)
df_bacon <- bacon(y ~ is_treated, data, "id", "t")

ggplot(df_bacon) +
  aes(x = weight, y = estimate, shape = type, color = type) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0) + 
  theme_bw() +
  labs(x = "Weight", y = "Estimate", shape = NULL, color = NULL) + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

ggsave(here("output/r/simulation/2_diagnosis/bacon_decomp.pdf"), width = 6, height = 4)

# Jakiela Diagnosis (https://pjakiela.github.io/TWFE/)
# Weight
model_tr_resid <- feglm(is_treated ~ 1 | id + t, data)
model_y_resid <- feglm(y ~ 1 | id + t, data)

df_jakiela <- data |>
  mutate(tr_resid = resid(model_tr_resid),
         y_resid = resid(model_y_resid))

df_jakiela |>
  mutate(label_wgt = case_when(
    t < tr_time ~ "Comparison observation",
    t >= tr_time & tr_resid > 0 ~ "Treatment observations - positive weight",
    t >= tr_time & tr_resid < 0 ~ "Treatment observations - negative weight"),
    label_wgt = factor(label_wgt, levels = c(
      "Comparison observation",
      "Treatment observations - positive weight",
      "Treatment observations - negative weight")),
    label_grp = recode_factor(group,
        `1` = "Treat in 1989",
        `2` = "Treat in 1998",
        `3` = "Treat in 2007")) |>
  ggplot(aes(x = t, y = label_grp, color = label_wgt)) +
  geom_point(size = 3, shape = "square") +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(values = c("grey", "darkgreen", "maroon")) +
  scale_x_continuous(breaks = seq(1980, 2015, 5)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank())

ggsave(here("output/r/simulation/2_diagnosis/jakiela_weight.pdf"), width = 8, height = 3)


# Heterogeneity

df_jakiela |>
  ggplot(aes(x = tr_resid, y = y_resid, color = is_treated)) +
  geom_point(shape = 1, size = 2.5, alpha = 0.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
  scale_color_manual(labels = c("Treatment observations", "Comparison observations"),
                     values = c("maroon", "darkgreen")) +
  labs(x = "D residual", y = "Y residual", color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.8, 0.9))

ggsave(here("output/r/simulation/2_diagnosis/jakiela_resid.pdf"), width = 6, height = 4)
