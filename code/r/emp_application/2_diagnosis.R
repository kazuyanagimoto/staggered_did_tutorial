library(tidyverse)
library(here)
library(fixest)

load(here("output/r/emp_application/1_overview/data.rds"))

data <- data |>
  filter(!is.na(y))

model_tr_resid <- feols(is_treated ~ 1 | year + pca_id^month + pca_id[log_load], data)
model_y_resid <- feols(y ~ 1 | year + pca_id^month + pca_id[log_load], data)

df_jakiela <- data |>
  mutate(tr_resid = resid(model_tr_resid),
         tr_resid2 = tr_resid^2,
         denom = sum(tr_resid2),
         w = tr_resid/denom,
         y_resid = resid(model_y_resid)) |>
  mutate(wm = mean(w), .by = c(group, time))
# Weight

df_jakiela |>
  mutate(label_wgt = case_when(
    time < tr_time ~ "Comparison observation",
    time >= tr_time & wm > 0 ~ "Treatment observations - positive weight",
    time >= tr_time & wm < 0 ~ "Treatment observations - negative weight"),
    label_wgt = factor(label_wgt, levels = c(
      "Comparison observation",
      "Treatment observations - positive weight",
      "Treatment observations - negative weight")),
    group = dense_rank(tr_time)) |>
  filter(!is.na(tr_time)) |>
  ggplot(aes(x = time, y = group, color = label_wgt)) +
  geom_point(size = 3, shape = "square") +
  labs(x = NULL, y = "Timing group", color = NULL) +
  scale_y_continuous(breaks = 1:17) +
  scale_color_manual(values = c("grey", "darkgreen", "maroon")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank())

ggsave(here("output/r/emp_application/2_diagnosis/jakiela_weight.png"), width = 6, height = 4)

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
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.82, 0.95))

ggsave(here("output/r/emp_application/2_diagnosis/jakiela_resid.png"), width = 6, height = 4)
