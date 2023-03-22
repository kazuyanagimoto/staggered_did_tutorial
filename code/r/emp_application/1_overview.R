library(here)
library(tidyverse)

cicala <- haven::read_dta(here("data/cicala_aer_2022_ready.dta"))

stata_ym <- function(year, month) {
    as.integer(12 * (year - 1960) + month - 1)
}

data <- cicala |>
  rename(y = log_ideal_trade_surplus,
         pca_id = pca_abbrev99) |>
  mutate(is_treated = market_operation == 1,
         tr_time = mkt_ym,
         time = stata_ym(year, month),
         date = lubridate::ym(str_c(year, "-", month)),
         rel_time = time - tr_time)

save(data, file = here("output/r/emp_application/1_overview/data.rds"))

# Check distribution of treatement over time

tibble(
  name = c("Num. of PCAs", "Num. of periods (months)", "Num. of timing groups",
           "Num. of obs. (daily)", "Fraction treated"),
  value = c(length(unique(data$pca_id)), length(unique(data$time)),
            length(unique(data$mkt_ym)), nrow(data),
            round(mean(data$is_treated), digits = 3))) |>
  mutate(value = as.character(value)) |>
  write_tsv(here("output/r/emp_application/1_overview/sum_stats1.tsv"))

data |>
  mutate(treat = if_else(is_treated, "Treat = 1", "Treat = 0")) |>
  select(imperf_netgen_cost, ideal_trade_surplus, treat) |>
  group_by(treat) |>
  summarize(across(everything(), list(mean = mean, sd = sd), .names = "{.fn}.{.col}")) |>
    write_tsv(here("output/r/emp_application/1_overview/sum_stats2.tsv"))




# Check how many timing groups?

data |>
  group_by(month = lubridate::floor_date(date, "month")) |>
  summarize(obs_by_month = n(),
            mkt_count = sum(is_treated),
            mkt_frac = mkt_count / obs_by_month) |>
  ggplot(aes(x = month, y = mkt_frac)) +
  geom_line() +
  labs(x = NULL, y = "Fraction of Market Dispatch") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(here("output/r/emp_application/1_overview/treate_timing.pdf"), width = 6, height = 6)

  
