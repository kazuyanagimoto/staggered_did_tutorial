library(here)
library(tidyverse)
library(patchwork)
library(microbenchmark)
library(multcomp)
library(broom)
# Estimator
library(fixest) # TWFE, Sun and Abraham
library(did) # Callaway-Sant'Anna
library(DIDmultiplegt) # de Chaisemartin and D'Haultfoeuille
library(didimputation) # Borusyak, Jaravel, Spiess

# Data
load(here("output/r/simulation/1_gen_data/data.rds"))
data <- data |>
  mutate(rel_time = t - tr_time)

for (i in 0:27) {
  v_name <- paste("lag", i, sep ='')
  data <- data |>
    mutate(!!v_name := ifelse(rel_time == i, 1, 0))
}

for (i in 1:27) {
  v_name <- paste("lead", i, sep ='')
  data <- data |>
    mutate(!!v_name := ifelse(rel_time == -i, 1, 0))
}

data <- data |>
  mutate(lead1 = 0)

## 0. Benchmark TWFE

est_naive <- function() {
  ctb <- feols(y ~ i(rel_time, ref = c(-1, 18:27)) | id + t,
               data = data, cluster = "id") |> coeftable()
  tibble(name = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
    add_row(rel_time = as.integer(-1), coef = 0, se = 0) |>
    mutate(lci = coef - 1.96 * se,
           rci = coef + 1.96 * se) |>
    filter(between(rel_time, -5, 5)) |>
    mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4))
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
    add_row(rel_time = as.integer(-1), coef = 0, se = 0) |>
    mutate(lci = coef - 1.96 * se,
           rci = coef + 1.96 * se) |>
    filter(between(rel_time, -5, 5)) |>
    mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4))
}


# (1-B) Sun and Abraham

est_sunab <- function() {
   ctb <- feols(y ~ sunab(tr_time, t, ref.c = max(tr_time)) | id + t,
                data = data, cluster = "id") |> aggregate(agg = "t")
   
   tibble(name = rownames(ctb),
          coef = ctb[,"Estimate"],
          se = ctb[,"Std. Error"]) |>
     separate(name, into = c(NA, "rel_time"), sep = "::", convert = TRUE) |>
     add_row(rel_time = as.integer(-1), coef = 0, se = 0) |>
     mutate(lci = coef - 1.96 * se,
            rci = coef + 1.96 * se) |>
     filter(between(rel_time, -5, 5)) |>
     mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4))
}

# (1-C) Wooldridge

est_wdrg <- function() {
  data <- data |> 
    mutate(across(lead1:lead27, list( g2 = \(x) x * (group == 2))),
           across(lead1:lead27, list( g3 = \(x) x * (group == 3))),
           across(lag0:lag27, list( g2 = \(x) x * (group == 2))),
           across(lag0:lag27, list( g3 = \(x) x * (group == 3))))
  xvars_lag <- names(data)[grep("^lag", names(data))]
  xvars_lead <- names(data)[grep("^lead", names(data))]
  formula <- as.formula(paste0("y", " ~ ", paste(xvars_lag, collapse = "+"), " + ",
                               paste(xvars_lead, collapse = "+"), " | t + group"))
  mod <- feols(formula, data = data, cluster = "id")
  
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
          se = mod$se.egt) |>
     filter(rel_time < 6, rel_time > -6) |>
     mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4),
            lci = coef - 1.96 * se,
            rci = coef + 1.96 * se)
}

# (2-B) de Chaisemartin and D'Haultfoeuille

est_dcdh <- function() {
  mod <- DIDmultiplegt::did_multiplegt(df = data,
                                       Y = "y", G = "id", T = "t", D = "is_treated",
                                       dynamic = 5, placebo = 5, brep = 100,
                                       cluster = "id")
  ests <- mod[grepl("^placebo_|^effect|^dynamic_", names(mod))]
  data.frame(rel_time = names(ests),
             coef = as.numeric(ests),
             se = as.numeric(mod[grepl("^se_placebo|^se_effect|^se_dynamic", names(mod))])) |> 
    mutate(rel_time = sub("effect", "0", rel_time)) |>
    mutate(rel_time = sub("placebo_", "-", rel_time)) |>
    mutate(rel_time = as.integer(sub("dynamic_", "", rel_time))) |>
    mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4),
             lci  = coef - se*1.96,
             rci = coef + se*1.96)
}

## 3. Imputation methods

# (3-A) Borusyak et al.

est_bjs <- function() {
  data <- data |> rename(dep_var = y)
  mod <- didimputation::did_imputation(data = data,
                                       yname = "dep_var",
                                       gname = "tr_time",
                                       tname = "t",
                                       idname = "id",
                                       horizon = TRUE,
                                       pretrends = -5:-1)
  tibble(rel_time = as.integer(mod$term),
                coef = mod$estimate,
                se = mod$std.error) |>
    filter(rel_time < 6, rel_time > -6) |>
    mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4),
           lci = coef - 1.96 * se,
           rci = coef + 1.96 * se)
}

# (3-B) Gardner

est_gard <- function() {
  t_d <- as.data.frame(i(data$t))
  colnames(t_d) <- paste0("t",names(t_d))
  id_d <- as.data.frame(i(data$id))
  colnames(id_d) <- paste0("id",names(id_d))
  data <- bind_cols(data, t_d, id_d)
  fe_vars <- names(data)[grep("^t[0-9]+|^id[0-9]+", names(data))]
  temp <- data[,c("y","is_treated",fe_vars)]
  mod_first <- lm(y ~ . - 1 - is_treated - t1980 - id1, data = subset(temp, is_treated == FALSE))
  data$pred <- predict(mod_first, newdata = data)
  data$resid <- data$y - data$pred
  xvars_lead <- names(data)[grep("^lead", names(data))]
  xvars_lag <- names(data)[grep("^lag", names(data))]
  formula <- as.formula(paste0("resid", " ~ ", paste(xvars_lead, collapse = "+"), " + ",
                               paste(xvars_lag, collapse = "+"), "-1"))
  ctb <- feols(formula, data = data, cluster = "id") |> coeftable()
  lead1 <- c(0,0,0,0)
  ctb <- rbind(ctb, lead1)
  tibble(rel_time = rownames(ctb),
         coef = ctb[,"Estimate"],
         se = ctb[,"Std. Error"]) |>
    mutate(rel_time = sub("lead", "-", rel_time)) |>
    mutate(rel_time = as.integer(sub("lag", "", rel_time))) |>
    filter(rel_time < 6, rel_time > -6) |>
    mutate(val_true = if_else(between(rel_time, -5, 0), 0, rel_time * 0.4),
           lci = coef - 1.96 * se,
           rci = coef + 1.96 * se)

}

#-------------------------------------------
# Plot & Benchmark
#-------------------------------------------

## Plot
plot_did <- function(ctb, title) {
  ctb |> ggplot(aes(x = rel_time, y = coef, ymin = lci, ymax = rci)) +
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

p1 <- est_naive()

|> plot_did("Naive TWFE Event Study")
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
    wdrg = est_wdrg(),
    clsa = est_clsa(),
    dcdh = est_dcdh(),
    bjs = est_bjs(),
    gard = est_gard()
)

# Yanagimoto will write here some exporting function of the benchmarks

