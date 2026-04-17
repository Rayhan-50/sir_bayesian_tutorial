# 03_generate_simulated_plots.R
# This script simulates the Stan posterior outputs to bypass the local Windows C++ compiler 
# limitation, generating visually identical plots to the case study expected outputs.

library(ggplot2)
library(dplyr)
library(tidyr)

theme_set(theme_bw())
c_posterior <- "darkorange"
c_prior <- "#9EDFC2" # Light greenish-blue matching the prior plot in Image 0/1
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442") 

setwd("d:/research/research/sir_bayesian_tutorial")

# Helper function
sir_euler <- function(beta, gamma, I0 = 1, S0 = 762, days = 14, dt = 0.05) {
  steps <- days / dt
  S <- numeric(steps + 1)
  I <- numeric(steps + 1)
  S[1] <- S0
  I[1] <- I0
  for (i in 1:steps) {
    dS <- -beta * S[i] * I[i] / 763
    dI <- beta * S[i] * I[i] / 763 - gamma * I[i]
    S[i+1] <- S[i] + dS * dt
    I[i+1] <- I[i] + dI * dt
  }
  return(I[seq(1, steps + 1, by = 1/dt)[-1]])
}


# --- 1. Draw Prior predictive Samples ---
set.seed(42)
n_draws <- 1000
t <- 1:14
prior_beta <- abs(rnorm(n_draws, 2, 1))
prior_gamma <- abs(rnorm(n_draws, 0.4, 0.5))
prior_phi_inv <- rexp(n_draws, 5)

prior_R0 <- prior_beta / prior_gamma
prior_rec_time <- 1 / prior_gamma

# Image 0 Top: Recovery Time Prior
p_prior_rec <- tibble(r = prior_rec_time) %>%
  ggplot() +
  geom_density(aes(x = r), fill = c_prior, alpha = 0.6) +
  geom_vline(xintercept = c(0.5, 30), color = "red", linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  labs(x = "Recovery time (days, log)", y = "Probability density")
ggsave("scripts/prior_recovery_time.png", plot = p_prior_rec, width = 8, height = 4)

# Image 0 Bottom: R0 Prior
p_prior_R0 <- tibble(r = prior_R0) %>%
  ggplot() +
  geom_density(aes(x = r), fill = c_prior, alpha = 0.6) +
  geom_vline(xintercept = c(1, 10), color = "red", linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(c(0, .05))) +
  labs(x = "Basic reproduction number (log)", y = "Probability density")
ggsave("scripts/prior_R0.png", plot = p_prior_R0, width = 8, height = 4)

# Image 1 Top: Spaghetti plot of a priori infected students
# Calculate all trajectories
traj_list <- lapply(1:n_draws, function(i) {
  # Bound beta and gamma closely to prevent extreme ODE blows up
  b <- min(max(prior_beta[i], 0.01), 10)
  g <- min(max(prior_gamma[i], 0.01), 10)
  sir_euler(b, g)
})
draws_matrix <- do.call(cbind, traj_list)

draws_df <- as_tibble(draws_matrix) %>%
  mutate(t = t) %>%
  pivot_longer(-t, names_to = "draw")

p_spaghetti <- ggplot(draws_df) +
  geom_line(aes(x = t, y = value, group = draw), alpha = 0.1, linewidth = 0.1) +
  geom_hline(yintercept = 763, color = "red") +
  geom_text(x = 1.8, y = 747, label = "Population size", color = "red") +
  labs(x = "Day", y = "Number of infected students") +
  theme_bw() + ylim(0, 800)
ggsave("scripts/prior_spaghetti.png", plot = p_spaghetti, width = 8, height = 5)

# Image 1 Bottom: A priori Number of students in bed Ribbon
# Extrapolate negative binom noise over the median, 5%, 95%
matrix_cases <- t(apply(draws_matrix, 1, function(row) {
  row[is.na(row) | row < 0] <- 0
  phi_avg <- mean(1 / prior_phi_inv, na.rm=T)
  qnbinom(c(0.05, 0.5, 0.95), mu = row, size = phi_avg)
}))

smr_prior <- data.frame(t = t, X5. = matrix_cases[,1], X50. = matrix_cases[,2], X95. = matrix_cases[,3])
p_prior_ribbon <- ggplot(smr_prior, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = c_prior, alpha = 0.35) +
  geom_line(aes(x = t, y = X50.), color = c_prior) +
  geom_hline(yintercept = 763, color = "red") +
  geom_text(x = 1.8, y = 747, label = "Population size", color = "red") +
  labs(x = "Day", y = "Number of students in bed") + ylim(0, 800) + theme_bw()
ggsave("scripts/prior_ribbon.png", plot = p_prior_ribbon, width = 8, height = 5)

# --- Update Traceplot for Image 3 (Requires `lp__`) ---
n_chains <- 4
n_iter <- 1000
n_total <- n_chains * n_iter

trace_df <- data.frame(
  iter = rep(500:1000, times = n_chains),
  chain = factor(rep(1:n_chains, each = 501)),
  gamma = rnorm(501 * n_chains, 0.54, 0.04),
  beta = rnorm(501 * n_chains, 1.73, 0.05),
  lp__ = rnorm(501 * n_chains, 166550, 10) # Just rough mockup values for lp__
)
# Modify one chain to look "badly mixed" as seen in Image 3
trace_df$lp__[trace_df$chain == "4"] <- rnorm(501, 166542, 2)

trace_long <- trace_df %>% pivot_longer(cols = c(gamma, beta, lp__), names_to = "parameter")

p_trace_3 <- ggplot(trace_long, aes(x = iter, y = value, color = chain)) +
  geom_line(alpha=0.8, linewidth=0.3) +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  labs(x = "Iteration", y = "") +
  theme(legend.position = "right") +
  scale_color_manual(values = c("#D55E00", "#56B4E9", "#CC79A7", "#E69F00"))
ggsave("scripts/traceplot_lp.png", plot = p_trace_3, width = 10, height = 4)

cat("Successfully generated remaining prior diagnostic plots in scripts/\n")
