library(dplyr)
library(ggplot2)
library(here)

set.seed(123)
n_obs <- 1000
T <- 36 # 1980~2015

df_ind <- tibble(
  id = 1:n_obs,
  ai = rnorm(n_obs, 0, 0.5),
  group = sample(1:3, n_obs, replace = TRUE, prob = c(17/50, 18/50, 15/50)),
  tr_time = recode(group, `1` = 1989, `2` = 1998, `3` = 2007)
)

data <- df_ind |>
 slice(rep(1:n(), each = T)) |>
 mutate(t = rep(1980:2015, time = n_obs),
        lt = rep(rnorm(T, 0, 0.5), time = n_obs),
        is_treated = t >= tr_time,
        y = ai + lt + rnorm(n_obs * T, 0, 0.5) + case_when(
          group == 1 ~ rnorm(n_obs * T, 0.5, 0.2) * (t - 1989) * is_treated,
          group == 2 ~ rnorm(n_obs * T, 0.3, 0.2) * (t - 1998) * is_treated,
          group == 3 ~ rnorm(n_obs * T, 0.1, 0.2) * (t - 2007) * is_treated),
        rel_time = t - tr_time,
        group = factor(group))

save(data, file = here("output/r/simulation/1_gen_data/data.rds"))

## Plot

data |>
  group_by(t, group) |>
  summarize(y_mean = mean(y), y_min = min(y), y_max = max(y), .groups = "drop") |>
  mutate(group = recode_factor(group,
                    `1` = "1989",
                    `2` = "1998",
                    `3` = "2007")) |>
  ggplot(aes(x = t, y = y_mean, ymin = y_min, ymax = y_max, color = group, fill = group)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(alpha = 0.1, color = NA) +
  geom_vline(xintercept = 1989.5, linetype = "dashed") +
  geom_vline(xintercept = 1997.5, linetype = "dashed") +
  geom_vline(xintercept = 2007.5, linetype = "dashed") +
  labs(x = NULL, y = NULL, color = NULL, fill = NULL,
       title = "Case 2. Heterogeneous/Dynamic Effect",
       subtitle = "Simulation 6 in Baker et al.") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.1, 0.9),
    plot.title.position = "plot")

ggsave(here("output/r/simulation/1_gen_data/baker22_sim6.pdf"), width = 6, height = 6)


## TWFE regression
feols(y ~ is_treated | id + t, data = data, cluster = "id")