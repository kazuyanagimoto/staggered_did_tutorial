library(ggplot2)
library(here)
library(bacondecomp)
library(fixest)

load(here("output/r/simulation/1_gen_data/data.rds"))

## Diagnosis
# Goodman-Bacon Decomposition (https://github.com/evanjflack/bacondecomp.git)
df_bacon <- bacon(y ~ is_treated, data, "id", "t")

ggplot(df_bacon) +
  aes(x = weight, y = estimate, shape = factor(type)) +
  geom_point() +
  geom_hline(yintercept = 0) + 
  theme_minimal() +
  labs(x = "Weight", y = "Estimate", shape = "Type")+
  theme(legend.title = element_blank())

# Jakiela Diagnosis (https://pjakiela.github.io/TWFE/)
# Weight
model_tr_resid = fixest::feglm(is_treated ~ 1 | id + t, data)

data <- data |>
  mutate(tr_resid = resid(model_tr_resid),
         tr_resid2 = resid(model_tr_resid)^2,
         denom = sum(tr_resid2),
         w = tr_resid/denom)

ggplot() + 
  # Treatment observations
  geom_point(data = data[data$t < data$tr_time,], aes(x = t, y = group), shape = 15, size = 4, color = "grey") +
  # Positive weight
  geom_point(data = data[data$t >= data$tr_time & data$tr_resid > 0,], aes(x = t, y = group), shape = 15, size = 4, color = "darkgreen") +
  # Negative weight
  geom_point(data = data[data$t >= data$tr_time & data$tr_resid < 0,], aes(x = t, y = group), shape = 15, size = 4, color = "maroon") +
  
  scale_y_continuous(breaks = c(1, 2, 3), labels = c("Treat in 1989", "Treat in 1998", "Treat in 2007"), expand = expansion(mult = 0.2)) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limits = c(1980, 2015), breaks = seq(1980, 2015, 5)) +
  
  theme(aspect.ratio = 0.1,
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 8),
        axis.ticks = element_blank())


# Heterogeneity
model_y_resid = fixest::feglm(y ~ 1 | id + t, data)

data <- data |>
  mutate(y_resid = resid(model_y_resid))

ggplot(data, aes(x = tr_resid, y = y_resid)) +
  geom_point(aes(color = factor(is_treated)), shape = 1, size = 2.5) +
  geom_smooth(aes(color = factor(is_treated)), method = "lm", se = FALSE, linewidth = 1, formula = y ~ x, linetype = 1) +
  scale_color_manual(labels = c("Treatment observations", "Comparison observations"), values = c("maroon", "darkgreen")) +
  
  xlab("D residual") +
  ylab("Y residual") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box = "horizontal",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))
