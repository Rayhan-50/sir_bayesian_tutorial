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

###  VS Code Terminal

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




Goal Description
The project will be heavily cleaned, refactored, and organized to focus solely on reproducing the Boarding School Influenza Bayesian SIR Stan case study. Files related to Swiss COVID SEIR models and redundant visualizations will be archived. Scripts will be renamed to professional standard names, with their internal code stripped of unrelated components, simplified, and nicely commented. Unimportant assets like installers and repetitive outputs will be relocated to an archive folder hierarchy.

User Review Required
WARNING

I intend to move the .docx and .pdf files to docs/ and name the main report final_report.docx. Please confirm if you want me to literally convert it to .pdf or if keeping it as .docx (or renaming to pdf without true conversion, which might corrupt it) is fine.

CAUTION

The scripts 00_master_reproduce.R and 06_summary_figure.R will have significant chunks of code removed (specifically, anything relating to the Swiss COVID-19 SEIR models) to meet the instruction to focus purely on the Boarding School dataset.

Proposed Changes
Structure & File Relocation
Create clean directories: outputs/figures, outputs/tables, docs, and archive/....
[DELETE] / [ARCHIVE]
We will move the following potentially unnecessary items to archive/:

rtools44_installer.exe (To archive/experimental/ to skip deleting)
create_docx.py (To archive/old_scripts/)
scripts/03_generate_simulated_plots.R, scripts/04_advanced_plots.R, scripts/05_final_plots.R (To archive/old_scripts/)
model/sir_negbin.rds (To archive/unused_models/)
Duplicate and repetitive PNG files in scripts tracking to archive/duplicate_plots/
[MODIFY] & [NEW] Cleaned Scripts
[NEW] scripts/fit_model.R
Renamed from scripts/01_fit_sir.R. Path dependencies updated (e.g., loading model as model/sir_negbin.stan and exporting to outputs/fitted_model.rds).
Saving figure to outputs/figures/posterior_predictive_check.png.
[NEW] scripts/diagnostics.R
Renamed and heavily whittled down from scripts/00_master_reproduce.R.
Removes COVID/SEIR content (Section 4, and poor-mixing traces from Section 3).
Modifies output paths to strictly write to outputs/figures/.
[NEW] scripts/generate_figures.R
Renamed from scripts/06_summary_figure.R.
Slices out COVID-19 plots (removing Panels H, J, K, L) and generates a clean, single summary plot for the SIR boarding school case study.
[NEW] README.md & cleanup_report.md
Generating fresh, comprehensive documentation outlining the project's layout, objectives, and running instructions.
Open Questions
Is it acceptable to keep the Stan files inside model/, and when rstan compiles them naturally it will create a .rds or compile it in memory, is that sufficient or do we explicitly need to handle compiled .rds output here?
Where should the summary table data (outputs/tables/summary_table.csv) originate from? I plan to append a short CSV export snippet in diagnostics.R that dumps print(fit_sir_negbin) statistical summaries to fulfill this requirement. Is this acceptable?
Verification Plan
Automated Tests
Run ls -R (via PowerShell equivalent) to confirm directory tree perfectly matches requested configuration.
Check file sizes/content of the refactored R scripts to ensure removal of Swiss COVID/SEIR code and correct file path references to outputs/.
Refactoring Bayesian SIR Project
Authentication Required
Please sign in.
