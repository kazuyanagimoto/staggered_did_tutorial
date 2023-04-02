library(here)
library(tidyverse)
library(patchwork)
library(microbenchmark)
# Estimator
library(fixest)
library(did)
library(DIDmultiplegt)
library(didimputation)
library(did2s)

load(here("output/r/emp_application/1_overview/data.rds"))
INF <- as.integer(1000)

## 0. Benchmark TWFE
est_naive <- function() {
  ctb <- data |>
    mutate(rel_time = replace_na(rel_time, INF)) |>
    feols(y ~ i(rel_time, ref = c(min(rel_time):-25, -1, 25:max(rel_time))) |
          year + pca_id^month + pca_id[log_load],
        data = _, weights = ~wgt, cluster = "pca_modate") |> coeftable()
  
  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0)
}

## 1. Fully-saturated TWFE

# (1-B) Sun and Abraham

est_sunab <- function() {
  ctb <- data |>
    mutate(tr_time = replace_na(tr_time, INF)) |>
    feols(y ~ sunab(tr_time, time, ref.p = c(-INF:-25, -1, 25:INF)) |
              year + pca_id^month + pca_id[log_load],
          data = _, weights = ~wgt, cluster = "pca_modate") |>
    aggregate(agg = "t")
  
  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0)
}

## 2. Rolling methods

# (2-A) Callaway-Sant'Anna

est_clsa <- function() {
  # make sure to code "never treated group" = 0
  data_ca <- data |> mutate(tr_time = replace(tr_time, is.na(tr_time), 0))
  
  # Absorb FEs and time-varying Covariates a la Caetano et al. (2022)
  data_ca <- data_ca |> mutate(pca_id_month = paste(pca_id, month))
  mod <- feols(y ~ -1 | year + pca_id_month + pca_id[log_load],
               data = subset(data_ca, is_treated=="FALSE"))
  data_ca$resid <- data_ca$y - predict(mod, data_ca)
  
  mod <- did::att_gt(yname = "resid",
                     tname = "time",
                     gname = "tr_time",
                     panel = FALSE,
                     clustervars = "pca_modate",
                     weightsname = "wgt",
                     control_group = "notyettreated",
                     data = data_ca) |> 
    aggte(type = "dynamic", na.rm = T)
  
  tibble(rel_time = mod$egt,
         coef = mod$att.egt,
         se = mod$se.egt)
}

# (2-B) de Chaisemartin and D'Haultfoeuille

est_dcdh <- function() {
  # Absorb FEs and time-varying Covariates a la Caetano et al. (2022)
  data_dcdh <- data |> mutate(pca_id_month = paste(pca_id, month))
  mod <- feols(y ~ 0 | year + pca_id_month + pca_id[log_load],
               data = subset(data_dcdh, is_treated=="FALSE"))
  data_dcdh$resid <- data_dcdh$y - predict(mod, data_dcdh)
  
  mod <- DIDmultiplegt::did_multiplegt(df = data_dcdh,
                                       Y = "resid",
                                       G = "pca_id",
                                       T = "time",
                                       D = "is_treated",
                                       dynamic = 18, placebo = 12, brep = 20) 
  # No firstdiff_placebo option as in Stata
  # No weight option
  # Error if we use cluster option (cluster = pca_id_month)
  
  ests <- mod[grepl("^placebo_|^effect|^dynamic_", names(mod))]
  
  tibble(rel_time = names(ests),
         coef = as.numeric(ests),
         se = as.numeric(mod[grepl("^se_*", names(mod))])) |>
    mutate(rel_time = str_replace(rel_time, "placebo_", "-"),
           rel_time = str_replace(rel_time, "effect", "0"),
           rel_time = str_replace(rel_time, "dynamic_", ""),
           rel_time = as.integer(rel_time))
  
}

## 3. Imputation methods

# (3-A) Borusyak et al.

est_bjs <- function() {
  mod <- data |>
    rename(dep_var = y) |> # BUG: We cannot use "y" as "yname"
    mutate(pca_id_month = paste(pca_id, month)) |>
    didimputation::did_imputation(yname = "dep_var",
                                  gname = "tr_time",
                                  tname = "time",
                                  idname = "pca_id",
                                  first_stage = ~ 0 | year + month + pca_id_month,
                                  wname = "wgt",
                                  horizon = TRUE,
                                  pretrends = -12:-1,
                                  cluster_var = "pca_modate")
  
  tibble(rel_time = as.integer(mod$term),
         coef = mod$estimate,
         se = mod$std.error)
}


# (3-B) Gardner

est_gard <- function() {
  ctb <- did2s::did2s(data,
                      yname = "y",
                      treatment = "is_treated",
                      cluster_var = "pca_modate",
                      first_stage = ~ 0 | year + pca_id^month + pca_id[log_load],
                      second_stage = ~ i(rel_time, ref = c(-1)),
                      weights = "wgt") |> coeftable()
    # Error if not using bootstrap

  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0)
}

## Plot
plot_did <- function(ctb, title) {
  ctb |>
    filter(between(rel_time, -12, 18)) |>
    mutate(lci = coef - 1.96 * se,
           rci = coef + 1.96 * se) |>
    ggplot(aes(x = rel_time, y = coef, ymin = lci, ymax = rci)) +
    geom_point(size = 2, color = "skyblue") +
    geom_line(color = "skyblue") +
    geom_line(aes(x = rel_time, y = coef), linetype = "longdash", color = "skyblue") +
    geom_vline(xintercept = 0.0, linetype = "longdash") +
    geom_hline(yintercept = 0.0) +
    geom_ribbon(alpha = 0.2, fill = "skyblue") +
    scale_x_continuous(breaks = seq(-12,18,6)) +
    labs(x = "Periods since the event",
         y = "Estimates",
         title = title) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title.position = "plot")
}

p1 <- est_naive() |> plot_did("Naive TWFE Event Study")
p2 <- est_sunab() |> plot_did("Sun and Abraham")
p3 <- est_clsa() |> plot_did("Callaway and Sant'Anna")
p4 <- est_dcdh() |> plot_did("de Chaisemartin and D'Haultfoeuille")
p5 <- est_bjs() |> plot_did("Borusyak, Jaravel, Spiess")
p6 <- est_gard() |> plot_did("Gardner")

(p1 + p2) / (p3 + p4)
ggsave(here("output/r/emp_application/3_estimation/alt_est1.pdf"), width = 6, height = 6)

(p1 + p2) / (p5 + p6)
ggsave(here("output/r/emp_application/3_estimation/alt_est2.pdf"), width = 6, height = 6)

## Benchmark
# I may not use microbenchmark here, since it might take too much time for some methods
mbm <- microbenchmark(
  sunab = est_sunab(),
  wdrg = est_wdrg(),
  clsa = est_clsa(),
  bjs = est_bjs(),
  gard = est_gard(),
  times = 100
)
mbm

mbm_dcdh <- microbenchmark(
  dcdh = est_dcdh(),
  times = 1
)
mbm_dcdh

