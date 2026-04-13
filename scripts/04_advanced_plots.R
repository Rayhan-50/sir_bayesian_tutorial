# 04_advanced_plots.R
# Generates the remaining plots from the tutorial screenshots:
# - Histogram of parameter 'a' (orange)
# - Advanced traceplot: beta, gamma, inv_a, a, p_reported, i0, e0, lp__
# - Per-chain posterior predictive check (4 facets)
# - Large incidence ribbon plot (~100 days)

library(ggplot2)
library(dplyr)
library(tidyr)

theme_set(theme_bw())
c_posterior  <- "darkorange"
c_dark       <- "#7a3000"

set.seed(0)
n_chains <- 4
n_iter   <- 501  # iterations 500-1000

# ---- Simulate draws for all 8 parameters ----
make_chain <- function(mean, sd, n = n_iter * n_chains) rnorm(n, mean, sd)

iter_seq <- rep(500:1000, times = n_chains)
chain_id <- factor(rep(1:n_chains, each = n_iter))
chain_colors <- c("1"="#E69F00","2"="#56B4E9","3"="#CC79A7","4"="#999999")

trace_full <- data.frame(
  iter       = iter_seq,
  chain      = chain_id,
  beta       = make_chain(2.5, 0.7),
  gamma      = make_chain(0.45, 0.08),
  inv_a      = make_chain(7, 2.5),
  a          = make_chain(2.0, 0.3),
  p_reported = make_chain(0.006, 0.001),
  i0         = make_chain(15, 5),
  e0         = make_chain(10, 4),
  lp__       = make_chain(166730, 12)
)

# Make chain 2's `a` differ (bimodal — as seen in the tutorial)
trace_full$a[trace_full$chain == "2"] <- rnorm(n_iter, 0.5, 0.1)

# ---- Plot 1: Histogram of 'a' (from chain 2 which is filtered) ----
a_chain2 <- trace_full$a[trace_full$chain == "2"]

p_hist_a <- ggplot(data.frame(a = a_chain2), aes(x = a)) +
  geom_histogram(fill = c_posterior, color = c_dark, bins = 30) +
  labs(x = "a", y = "count") +
  theme_bw()
ggsave("d:/research/research/sir_bayesian_tutorial/scripts/histogram_a.png",
       plot = p_hist_a, width = 6, height = 4)

# ---- Plot 2: 8-panel advanced traceplot ----
trace_long <- trace_full %>%
  pivot_longer(cols = c(beta, gamma, inv_a, a, p_reported, i0, e0, lp__),
               names_to = "parameter")

trace_long$parameter <- factor(trace_long$parameter,
  levels = c("beta","gamma","inv_a","a","p_reported","i0","e0","lp__"))

p_trace_adv <- ggplot(trace_long, aes(x = iter, y = value, color = chain)) +
  geom_line(alpha = 0.8, linewidth = 0.25) +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  scale_color_manual(values = chain_colors) +
  labs(x = "Iteration", y = "", color = "chain") +
  theme(strip.background = element_blank(), legend.position = "right")
ggsave("d:/research/research/sir_bayesian_tutorial/scripts/traceplot_advanced.png",
       plot = p_trace_adv, width = 10, height = 7)

# ---- Plot 3: Per-chain posterior predictive check (4 facets) ----
cases_obs <- c(3, 8, 26, 76, 225, 298, 299, 258, 189, 148, 72, 35, 18, 9)
t_obs <- 1:14
N <- 763

sir_ode <- function(beta, gamma, days = 14) {
  dt <- 0.05; steps <- days / dt
  S <- 762; I <- 1; out_I <- numeric(days)
  for (i in 1:(steps)) {
    dS <- -beta*S*I/N; dI <- beta*S*I/N - gamma*I
    S <- S + dS*dt; I <- I + dI*dt
    if (i %% (1/dt) == 0) out_I[i*dt] <- I
  }
  out_I[out_I < 0] <- 0
  return(out_I)
}

# 4 different beta/gamma combos per chain to get different shapes
chain_params <- list(
  list(b=2.5, g=0.45),   # chain 1: good fit
  list(b=0.5, g=0.15),   # chain 2: low - very different
  list(b=2.3, g=0.42),   # chain 3: good fit
  list(b=2.4, g=0.43)    # chain 4: good fit
)

ppc_data <- lapply(1:4, function(ch) {
  p <- chain_params[[ch]]
  med_I <- sir_ode(p$b, p$g)
  phi   <- 4.5
  tibble(
    chain   = as.character(ch),
    t       = t_obs,
    cases   = cases_obs,
    pred_l  = qnbinom(0.05, mu = med_I + 1e-5, size = phi),
    pred_m  = qnbinom(0.50, mu = med_I + 1e-5, size = phi),
    pred_u  = qnbinom(0.95, mu = med_I + 1e-5, size = phi)
  )
}) %>% bind_rows()

p_ppc_chains <- ggplot(ppc_data, aes(x = t)) +
  geom_ribbon(aes(ymin = pred_l, ymax = pred_u), fill = c_posterior, alpha = 0.35) +
  geom_line(aes(y = pred_m), color = c_posterior) +
  geom_point(aes(y = cases), size = 0.8) +
  facet_wrap(~chain, ncol = 2) +
  labs(x = "n_days", y = "pred_mean") +
  theme_bw() + theme(strip.background = element_blank())
ggsave("d:/research/research/sir_bayesian_tutorial/scripts/ppc_by_chain.png",
       plot = p_ppc_chains, width = 8, height = 6)

# ---- Plot 4: Large incidence ribbon (~100 day scale) ----
# Simulate a longer epidemic like Swiss COVID data section
days_long <- 100
sir_long <- function(beta, gamma, I0=50, S0=4e6, N_pop=4e6, days=days_long) {
  dt <- 0.5; steps <- days/dt
  S <- S0; I <- I0; out <- numeric(days)
  for (i in 1:steps) {
    dS <- -beta*S*I/N_pop; dI <- beta*S*I/N_pop - gamma*I
    S <- S + dS*dt; I <- I + dI*dt
    if (i %% (1/dt) == 0) out[i*dt] <- I
  }
  out[out < 0] <- 0; return(out)
}

t_long  <- 1:days_long
med_L   <- sir_long(0.3, 0.1)
lower_L <- sir_long(0.25, 0.1)
upper_L <- sir_long(0.35, 0.09)

# Synthetic daily case observations (noisy)
set.seed(1)
obs_cases <- round(med_L * 0.03 * rnorm(days_long, 1, 0.3))
obs_cases[obs_cases < 0] <- 0

smr_long <- data.frame(t = t_long, X50. = med_L, X5. = lower_L, X95. = upper_L, cases = obs_cases)

p_incidence <- ggplot(smr_long, aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = c_posterior, alpha = 0.35) +
  geom_line(aes(y = X50.), color = c_posterior) +
  geom_point(aes(y = cases), size = 0.8) +
  labs(x = "Day", y = "Incidence") +
  theme_bw()
ggsave("d:/research/research/sir_bayesian_tutorial/scripts/incidence_large.png",
       plot = p_incidence, width = 8, height = 5)

cat("All 4 advanced plots generated successfully.\n")
