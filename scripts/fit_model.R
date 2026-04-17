# fit_model.R  ---------------------------------------------------------------
# Exact replication of the mc-stan.org boarding school case study
# Reference: https://mc-stan.org/learn-stan/case-studies/boarding_school_case_study.html
# -----------------------------------------------------------------------------

library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
theme_set(theme_bw())

# Color definitions matching the tutorial
c_posterior <- "darkorange"

setwd("d:/research/research/sir_bayesian_tutorial")

# --- 1. Load Data -----------------------------------------------------------
# The outbreaks::influenza_england_1978_school dataset (14 days, 763 students)
cases <- read.csv("data/influenza_england_1978_school.csv")$in_bed
N     <- 763

# --- 2. Build Stan data list ------------------------------------------------
n_days <- length(cases)

# Time vector: seq(0, n_days) then drop t=0 (the initial time)
t  <- seq(0, n_days, by = 1)
t0 <- 0
ts <- t[-1]          # observation times: 1, 2, ..., 14

# Initial conditions (NO names — just a plain numeric vector)
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 <- c(s0, i0, r0)   # c(762, 1, 0)

data_sir <- list(
  n_days = n_days,
  y0     = y0,
  t0     = t0,
  ts     = ts,
  N      = N,
  cases  = cases
)

# --- 3. Compile Stan model --------------------------------------------------
cat("Compiling Stan model...\n")
model <- stan_model("model/sir_negbin.stan")

# --- 4. Sample --------------------------------------------------------------
cat("Sampling... (4 chains x 2000 iter)\n")
fit_sir_negbin <- sampling(
  model,
  data   = data_sir,
  iter   = 2000,
  chains = 4,
  seed   = 0
)

# --- 5. Print summary -------------------------------------------------------
pars <- c("beta", "gamma", "R0", "recovery_time")
cat("\n----- Parameter Summary -----\n")
print(fit_sir_negbin, pars = pars)

# --- 6. Save fit object -----------------------------------------------------
saveRDS(fit_sir_negbin, "outputs/fitted_model.rds")
cat("Fit saved to outputs/fitted_model.rds\n")

# --- 7. Posterior predictive check plot ------------------------------------
# Matches tutorial code exactly:
#   smr_pred <- cbind(as.data.frame(summary(fit, pars="pred_cases",
#                     probs=c(0.05,0.5,0.95))$summary), t=ts, cases=cases)
smr_pred <- cbind(
  as.data.frame(
    summary(fit_sir_negbin,
            pars  = "pred_cases",
            probs = c(0.05, 0.5, 0.95))$summary
  ),
  t     = ts,
  cases = cases
)
colnames(smr_pred) <- make.names(colnames(smr_pred))  # remove % from names

# Plot: ribbon (90% CI) + median line + observed points
p <- ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.),
              fill  = c_posterior,
              alpha = 0.35) +
  geom_line(mapping = aes(x = t, y = X50.),
            color = c_posterior) +
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Number of students in bed")

ggsave("outputs/figures/posterior_predictive_check.png", plot = p, width = 7, height = 4)
cat("Plot saved to outputs/figures/posterior_predictive_check.png\n")
