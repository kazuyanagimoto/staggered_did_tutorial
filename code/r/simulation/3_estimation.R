library(here)
library(tidyverse)
library(patchwork)
library(microbenchmark)
library(multcomp) # May not be necessary
library(broom) # May not be necessary
# Estimator
library(fixest) # TWFE, Sun and Abraham, Wooldridge (Manual implementation)
library(did) # Callaway-Sant'Anna
library(DIDmultiplegt) # de Chaisemartin and D'Haultfoeuille
library(didimputation) # Borusyak, Jaravel, Spiess
library(did2s) # Gardner

# Data
load(here("output/r/simulation/1_gen_data/data.rds"))

## 0. Benchmark TWFE

est_naive <- function() {
  ctb <- feols(y ~ i(rel_time, ref = c(-1, 18:27)) | id + t,
               data = data, cluster = "id") |> coeftable()
  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0)
}

## 1. Fully-saturated TWFE

# (1-A) Stacked regression

est_stack <- function() {
  ctb <- data |>
    filter(
      between(t, 1989 - 5, 1989 + 5) |
      (between(t, 1998 - 5, 1998 + 5) & group == 2) |
      (between(t, 2007 - 5, 2007 + 5) & group == 3)) |>
    feols(y ~ i(rel_time, ref = c(-1, -27:-6, 15:27)) | id + t,
          data = _, cluster = "id") |> coeftable()
  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0)
}


# (1-B) Sun and Abraham

est_sunab <- function() {
   ctb <- feols(y ~ sunab(tr_time, t, ref.c = max(tr_time)) | id + t,
                data = data, cluster = "id") |> aggregate(agg = "t")
   
   tibble(name = rownames(ctb),
          coef = ctb[,"Estimate"],
          se = ctb[,"Std. Error"]) |>
     separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
     add_row(rel_time = as.integer(-1), coef = 0, se = 0)
}

# (1-C) Wooldridge

est_wdrg <- function() {

  ctb <- feols(y ~ i(rel_time, group, ref = c(-1, 27), ref2 = 3) | t + group,
              data = data, cluster = "id") |>
        coeftable()

  wtb <- tibble(
    group = factor(1:2),
    wgt = c(mean(data$group == 1), mean(data$group == 2))
  )

  tibble(name = rownames(ctb),
          coef = ctb[,"Estimate"],
          se = ctb[,"Std. Error"]) |>
     separate(name, into = c(NA, "rel_time", "group"), sep = "::") |>
     separate(rel_time, into = c("rel_time", NA), sep = ":", convert = TRUE) |>
     filter(between(rel_time, -5, 5)) |>
     add_row(rel_time = as.integer(-1), group = "1", coef = 0, se = 0) |>
     add_row(rel_time = as.integer(-1), group = "2", coef = 0, se = 0) |>
     left_join(wtb, by = "group") |>
     group_by(rel_time) |>
     summarize(coef = weighted.mean(coef, wgt),
               se = sqrt(weighted.mean(se^2, wgt^2)))

}


est_wdrg_original <- function() {

  data_wdrg <- data

  for (i in 0:27) {
    v_name <- paste("lag", i, sep ='')
    data_wdrg <- data_wdrg |>
      mutate(!!v_name := ifelse(rel_time == i, 1, 0))
  }

  for (i in 1:27) {
    v_name <- paste("lead", i, sep ='')
    data_wdrg <- data_wdrg |>
      mutate(!!v_name := ifelse(rel_time == -i, 1, 0))
  }

  data_wdrg <- data_wdrg |>
    mutate(lead1 = 0)


  data_wdrg <- data_wdrg |> 
    mutate(across(lead1:lead27, list( g2 = \(x) x * (group == 2))),
           across(lead1:lead27, list( g3 = \(x) x * (group == 3))),
           across(lag0:lag27, list( g2 = \(x) x * (group == 2))),
           across(lag0:lag27, list( g3 = \(x) x * (group == 3))))
  xvars_lag <- names(data)[grep("^lag", names(data))]
  xvars_lead <- names(data)[grep("^lead", names(data))]
  formula <- as.formula(paste0("y", " ~ ", paste(xvars_lag, collapse = "+"), " + ",
                               paste(xvars_lead, collapse = "+"), " | t + group"))
  mod <- feols(formula, data = data_wdrg, cluster = "id")
  
  list <- list()
  for (i in 0:5) {
    coefs <- mod$coefficient
    var1 <- ifelse(length(grep(paste0("^lag", i, "TRUE"), names(coefs))) != 0,
                   names(coefs)[grep(paste0("^lag", i, "TRUE"), names(coefs))], 0)
    var2 <- ifelse(length(grep(paste0("^lag", i, "_g2"), names(coefs))) != 0,
                   names(coefs)[grep(paste0("^lag", i, "_g2"), names(coefs))], 0)
    var3 <- ifelse(length(grep(paste0("^lag", i, "_g3"), names(coefs))) != 0,
                   names(coefs)[grep(paste0("^lag", i, "_g3"), names(coefs))], 0)
    hypo <- paste0(var1, " + ", var2, " * (2/3) + ", var3, " * (1/3) = 0")
    mod.lh <- glht(mod, linfct = c(hypo))
    mod.conf <- confint(mod.lh)
    list <- c(list, mod.conf$confint, as.integer(i))
  }
  
  list <- c(list, 0, 0, 0, as.integer(-1))
  for (i in 2:5) {
    coefs <- mod$coefficient
    var1 <- ifelse(length(grep(paste0("^lead", i, "TRUE"), names(coefs))) != 0,
                   names(coefs)[grep(paste0("^lead", i, "TRUE"), names(coefs))], 0)
    var2 <- ifelse(length(grep(paste0("^lead", i, "_g2"), names(coefs))) != 0,
                   names(coefs)[grep(paste0("^lead", i, "_g2"), names(coefs))], 0)
    var3 <- ifelse(length(grep(paste0("^lead", i, "_g3"), names(coefs))) != 0,
                   names(coefs)[grep(paste0("^lead", i, "_g3"), names(coefs))], 0)
    hypo <- paste0(var1, " + ", var2, " * (2/3) + ", var3, " * (1/3) = 0")
    mod.lh <- glht(mod, linfct = c(hypo))
    mod.conf <- confint(mod.lh)
    list <- c(list, mod.conf$confint, as.integer(-i))
  }
  
  unlist(list) |> matrix(ncol = 4, byrow = T) |> data.frame() |>
      rename(coef = X1, rci = X2, lci = X3, rel_time = X4) |>
      mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4))
}


## 2. Rolling methods

# (2-A) Callaway-Sant'Anna

est_clsa <- function() {
   mod <- did::att_gt(yname = "y",
                      tname = "t",
                      idname = "id",
                      gname = "tr_time",
                      control_group = "notyettreated",
                      data = data) |> aggte(type = "dynamic")
   tibble(rel_time = mod$egt,
          coef = mod$att.egt,
          se = mod$se.egt)
}


# (2-B) de Chaisemartin and D'Haultfoeuille

est_dcdh <- function() {
  mod <- DIDmultiplegt::did_multiplegt(df = data,
                                       Y = "y", G = "id", T = "t", D = "is_treated",
                                       dynamic = 5, placebo = 5, brep = 100,
                                       cluster = "id", parallel = TRUE)
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
          didimputation::did_imputation(yname = "dep_var",
                                       gname = "tr_time",
                                       tname = "t",
                                       idname = "id",
                                       horizon = TRUE,
                                       pretrends = -5:-1)
  tibble(rel_time = as.integer(mod$term),
         coef = mod$estimate,
         se = mod$std.error)
}

# (3-B) Gardner

est_gard <- function() {
  ctb <- did2s::did2s(data, yname = "y", treatment = "is_treated", cluster_var = "id",
             first_stage = ~ 0 | id + t,
             second_stage = ~i(rel_time, ref = c(-1))) |> coeftable()
  tibble(name = rownames(ctb),
          coef = ctb[,"Estimate"],
          se = ctb[,"Std. Error"]) |>
      separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
      add_row(rel_time = as.integer(-1), coef = 0, se = 0)
}


#-------------------------------------------
# Plot & Benchmark
#-------------------------------------------

## Plot
plot_did <- function(ctb, title) {
  ctb |>
    filter(between(rel_time, -5, 5)) |>
    mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4),
           lci = coef - 1.96 * se,
           rci = coef + 1.96 * se) |>
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

p1 <- est_naive() |> plot_did("Naive TWFE Event Study")
p2 <- est_stack() |> plot_did("Stacked Regression")
p3 <- est_sunab() |> plot_did("Sun and Abraham")
p4 <- est_wdrg() |> plot_did("Wooldridge")
p5 <- est_clsa() |> plot_did("Callaway and Sant'Anna")
p6 <- est_dcdh() |> plot_did("de Chaisemartin and D'Haultfoeuille")
p7 <- est_bjs() |> plot_did("Borusyak, Jaravel, Spiess")
p8 <- est_gard() |> plot_did("Gardner")

(p1 + p2) / (p3 + p4)
ggsave(here("output/r/simulation/3_estimation/alt_est1.pdf"), width = 6, height = 6)
(p5 + p6) / (p7 + p8)
ggsave(here("output/r/simulation/3_estimation/alt_est2.pdf"), width = 6, height = 6)



## Benchmark

mbm <- microbenchmark(
    sunab = est_sunab(),
    clsa = est_clsa(),
    bjs = est_bjs(),
    gard = est_gard(),
    times = 100
)

mbm_dcdh <- microbenchmark(
  dcdh = est_dcdh(),
  times = 10
)

tb_mbm <- summary(mbm) |>
  dplyr::select(method = expr, time = median, num_eval = neval) |>
  mutate(time = time / 1000)

tb_mbm_dcdh <- summary(mbm_dcdh) |>
  dplyr::select(method = expr, time = median, num_eval = neval)

tb_mbm |>
  bind_rows(tb_mbm_dcdh) |>
  mutate(method = recode_factor(method,
    sunab = "Sun and Abraham",
    clsa = "Callway and Sant'Anna",
    dcdh = "de Chaisemartin and D'Haultfoeuille",
    bjs = "Borusyak, Jaravel, Spiess",
    gard = "Gardener"),
    time = sprintf("%02.f:%02.f:%03.f",
                   time %/% 60,
                   round(time) %% 60,
                   round(time * 1000) %% 1000)) |>
  arrange(method) |>
  write_tsv(here("output/r/simulation/3_estimation/bench_sim.tsv"))
