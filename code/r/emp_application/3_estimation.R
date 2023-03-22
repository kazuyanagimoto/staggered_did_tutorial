library(here)
library(tidyverse)
library(patchwork)
library(microbenchmark)

library(fixest)

load(here("output/r/emp_application/1_overview/data.rds"))
INF <- as.integer(1000)

## 0. Benchmark TWFE
est_naive <- function() {
  data |>
    mutate(rel_time = replace_na(rel_time, INF)) |>
    feols(
      y ~ i(rel_time, ref = c(min(rel_time):-25, -1, 25:max(rel_time))) |
          year + pca_id^month + pca_id[log_load],
        data = _, weights = ~wgt, cluster = "pca_modate")
}

## 1. Fully-saturated TWFE

# (1-B) Sun and Abraham


est_sunab <- function() {
  data |>
    mutate(tr_time = replace_na(tr_time, INF)) |>
    feols(y ~ sunab(tr_time, time, ref.p = c(-INF:-25, -1, 25:INF)) |
              year + pca_id^month + pca_id[log_load],
          data = _, weights = ~wgt, cluster = "pca_modate")
}


## Plot
plot_ctb <- function(ctb, title) {
  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0) |>
    mutate(lci = coef - 1.96 * se,
          rci = coef + 1.96 * se) |>
  filter(between(rel_time, -12, 18)) |>
  ggplot(aes(x = rel_time, y = coef, ymin = lci, ymax = rci)) +
  geom_point(size = 2, color = "skyblue") +
  geom_line(color = "skyblue") +
  geom_vline(xintercept = 0.0, linetype = "longdash") +
  geom_hline(yintercept = 0.0) +
  geom_ribbon(alpha = 0.2, fill = "skyblue") +
  scale_x_continuous(breaks = seq(-12, 18, by = 3)) +
  labs(x = "Months since the event",
       y = "Estimates",
       title = title) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot")
}

p1 <- est_naive() |> coeftable() |> plot_ctb("Naive TWFE Event Study")
p2 <- est_sunab() |> aggregate(agg = "t") |> plot_ctb("Sun and Abraham")
p3 <- p1 # Replace it with Callaway and Sant'Anna
p4 <- p2 # Replace it with de Chaisemartin and D'Haultfoeuille
p5 <- p1 # Replace it with Borusyk, Jaravel, and Spiess
p6 <- p2 # Replace it with Gardener

p2 <- res_sunab |> aggregate(agg = "t") |> plot_ctb("Sun and Abraham")
(p1 + p2) / (p3 + p4)
ggsave(here("output/r/emp_application/3_estimation/alt_est1.pdf"), width = 6, height = 6)

(p1 + p2) / (p5 + p6)
ggsave(here("output/r/emp_application/3_estimation/alt_est2.pdf"), width = 6, height = 6)

## Benchmark
# I may not use microbenchmark here, since it might take too much time for some methods
mbm <- microbenchmark(
    sunab = est_sunab(),
    times = 10
)
mbm
