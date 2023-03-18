library(here)
library(tidyverse)
library(patchwork)
library(microbenchmark)

# Estimator
library(fixest) # TWFE, Sun and Abraham

# Data
load(here("output/r/simulation/1_gen_data/data.rds"))
data <- data |>
  mutate(rel_time = t - tr_time)

## 0. Benchmark TWFE

est_naive <- function() {
  feols(y ~ i(rel_time, ref = c(-1, 18:27)) | id + t, data = data, cluster = "id")
}

## 1. Fully-saturated TWFE

# (1-A) Stacked regression

est_stack <- function() {
  data |>
    filter(
      between(t, 1989 - 5, 1989 + 5) |
      (between(t, 1998 - 5, 1998 + 5) & group == 2) |
      (between(t, 2007 - 5, 2007 + 5) & group == 3)) |>
    feols(y ~ i(rel_time, ref = c(-1, -27:-6, 15:27)) | id + t,
          data = _, cluster = "id")
}


# (1-B) Sun and Abraham

est_sunab <- function() {
   feols(y ~ sunab(tr_time, t, ref.c = max(tr_time)) | id + t,
             data = data, cluster = "id")
}


# (1-C) Wooldridge

## 2. Rolling methods

# (2-A) Callaway-Sant'Anna

# (2-B) de Chaisemartin and D'Haultfoeuille


## 3. Inputation methods

# (3-A) Borusyak et al.

# (3-B) Gardner


#-------------------------------------------
# Plot & Benchmark
#-------------------------------------------

## Plot
plot_ctb <- function(ctb, title) {
  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0) |>
    mutate(lci = coef - 1.96 * se,
          rci = coef + 1.96 * se) |>
  filter(between(rel_time, -5, 5)) |>
  mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4)) |>
  ggplot(aes(x = rel_time, y = coef, ymin = lci, ymax = rci)) +
  geom_point(size = 2, color = "skyblue") +
  geom_line(color = "skyblue") +
  geom_line(aes(x = rel_time, y = val_true), linetype = "longdash", color = "orange") +
  geom_vline(xintercept = 0.0, linetype = "longdash") +
  geom_hline(yintercept = 0.0) +
  geom_ribbon(alpha = 0.2, fill = "skyblue") +
  scale_x_continuous(breaks = -5:5) +
  labs(x = "Periods since the event",
       y = "Estimates",
       title = title) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot")
}

p1 <- est_naive() |> coeftable() |> plot_ctb("Naive TWFE Event Study")
p2 <- est_stack() |> coeftable() |> plot_ctb("Stacked Regression")
p3 <- est_sunab() |> aggregate(agg = "t") |> plot_ctb("Sun and Abraham")
#p4 <- WRITE HERE


(p1 + p2) / (p3 + p3)
ggsave(here("output/r/simulation/3_estimation/alt_est1.pdf"), width = 6, height = 6)

## Benchmark

mbm <- microbenchmark(
    sunab = est_sunab()
)

# Yanagimoto will write here some exporting function of the benchmarks



