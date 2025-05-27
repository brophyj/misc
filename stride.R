# Setup
  

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tibble)


## Simulate STRIDE-like Data


set.seed(42)
n <- 396

baseline_sema <- rgamma(n, shape = 2, scale = 92)
baseline_placebo <- rgamma(n, shape = 2, scale = 92)
ratio <- c(Semaglutide = 1.21, Placebo = 1.08)
followup_sema <- baseline_sema * ratio['Semaglutide'] + rnorm(n, 0, 30)
followup_placebo <- baseline_placebo * ratio['Placebo'] + rnorm(n, 0, 30)

df <- tibble(
  group = rep(c("Semaglutide", "Placebo"), each = n),
  baseline = c(baseline_sema, baseline_placebo),
  followup = c(followup_sema, followup_placebo)
) %>%
  mutate(change = followup - baseline,
         treat = ifelse(group == "Semaglutide", 1, 0))


## Naive Frequentist Analysis


mean_change <- df %>% group_by(group) %>% summarise(mean = mean(change))
mean_diff <- diff(mean_change$mean)
pval <- t.test(change ~ group, data = df)$p.value
mean_change
mean_diff
pval


## Frequentist ANCOVA


ancova <- lm(followup ~ treat + baseline, data = df)
summary(ancova)


## Bayesian ANCOVA


stan_code <- "
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x_baseline;
  vector[N] x_group;
}
parameters {
  real alpha;
  real beta_baseline;
  real beta_group;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha + beta_baseline * x_baseline + beta_group * x_group, sigma);
}
generated quantities {
  real adj_diff = beta_group;
  real prob_gt_20 = normal_cdf(adj_diff - 20 | 0, sigma);
}
"

file_path <- write_stan_file(stan_code)
model <- cmdstan_model(file_path)
fit <- model$sample(
  data = list(
    N = nrow(df),
    y = df$followup,
    x_baseline = df$baseline,
    x_group = df$treat
  ),
  seed = 123,
  chains = 4,
  iter_sampling = 2000,
  iter_warmup = 1000,
  refresh = 0
)


## Posterior Summary and Visualization


draws <- fit$draws(c("adj_diff"))
adj_diff_draws <- as_draws_df(draws)$adj_diff
prob_gt_20 <- mean(adj_diff_draws > 20)
summary_stats <- summarise_draws(draws)
summary_stats


mcmc_areas(as_draws_matrix(draws), pars = "adj_diff", prob = 0.95) +
  geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
  ggtitle("Posterior of Adjusted Treatment Effect") +
  labs(subtitle = paste0("Red dashed line = clinically meaningful threshold (20 m), P(Δ > 20) = ", round(prob_gt_20, 3)))
ggsave("posterior_adjusted_effect.png", width = 8, height = 5)

