# 1978 Boarding School Influenza SIR Model

This project implements a classic Epidemiological SIR (Susceptible-Infected-Recovered) model using Bayesian inference in Stan and R. It analyzes the famous 1978 outbreak of influenza at a British boarding school.

## Project Structure

- `data/influenza_england_1978_school.csv`: The daily count of students confined to bed.
- `model/sir_negbin.stan`: The Bayesian SIR model written in Stan, using a Negative Binomial likelihood, and an ODE solver.
- `scripts/01_fit_sir.R`: The complete analysis script using `cmdstanr` to load the data, compile/run the model, and plot a posterior predictive check.

## Prerequisites

To run this project, make sure you have R installed and the following packages:

```R
install.packages(c("posterior", "dplyr", "tidyr", "ggplot2", "cmdstanr"))
```

You must also install the CmdStan backend. If you haven't yet, you can do this from R:

```R
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
install_cmdstan()
```

## How to Run the Project

### Using Antigravity / VS Code Terminal

1. Open this repository (`d:\research\research\sir_bayesian_tutorial`) in your IDE or terminal.
2. Ensure you have the proper R environment active.
3. Depending on your OS, simply execute the R script from the command line:

```bash
# In your terminal
Rscript scripts/01_fit_sir.R
```

The script will:
- Compile the Stan model (`sir_negbin.stan`).
- Sample the posterior distribution using 4 parallel chains.
- Print a summary of the model parameters (`beta`, `gamma`, `R0`, `recovery_time`).
- Save the results as an RDS file (`scripts/fit_sir.rds`).
- Generate a posterior predictive plot with an orange ribbon for the 90% credible interval and black dots for the observed data, saving it as `scripts/posterior_predictive_plot.png`.
