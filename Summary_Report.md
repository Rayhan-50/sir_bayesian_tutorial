# Summary of Bayesian Disease Transmission Modeling

This document outlines the step-by-step findings from our reproduction of the Bayesian disease transmission modeling case study (Grinsztajn et al.). Our results are visualized in the `summary_figure.png` composite chart which is divided into twelve specific analytical steps.

### Step A: Raw Observed Data
The first step maps the empirical data from the 1978 Boarding School influenza outbreak. Showing the number of students confined to bed across January and February gives us the baseline outbreak timeline that our modeling seeks to recreate.

### Step B: SIR Posterior Predictive Check
Here, we overlay our basic SIR model's 95% credible interval posteriors atop the actual observed raw data dots. The tight ribbon tracking the daily reported cases confirms that our deterministic SIR structure successfully captures the fundamental disease mechanics of the outbreak.

### Step C: Latent Infected Students
We project the true, unobserved number of infected students within the population—many of whom may not be confined to bed yet. Because our model accounts for the mathematical gap between being infectious and being reported, we provide a robust estimate of the true outbreak size hidden beneath the raw reporting.

### Step D: Prior Expectations for Recovery Time
Before feeding data to our model, we establish wide, physically plausible prior bounds for the infection recovery time using log scales. The sensitivity checks (red dashed lines) confirm that our mathematical priors mathematically restrict the model into a medically realistic range between half a day to thirty days.

### Step E: Prior Expectations for Reproduction Number ($R_0$)
We outline our initial assumption distribution for the basic reproduction number ($R_0$). By setting the density explicitly between sensible boundaries (between 1 and 10), we ensure that the MCMC sampler doesn't waste time exploring physically impossible disease exponential rates.

### Step F: Prior Predictive Ribbon (Spaghetti Plot)
We feed our initial priors back into the simulation to forecast the epidemic without any knowledge of the real-world dataset. Providing a widely varying span of possible trajectories ensures that our priors are sufficiently broad and non-restrictive before being clamped down by the actual data evidence.

### Step G: Traceplot - Well-Mixed Inference Space (SIR)
We conduct MCMC diagnostics for the simple SIR model by observing four independent sampling chains. Because all colored chains are tightly woven together—resembling a fuzzy, static caterpillar—we statistically prove a healthy convergence on the global optimal parameter values.

### Step H: Traceplot - Poorly-Mixed Diagnostic Warning (SEIR)
We introduce an explicit error case where an isolated MCMC sampling chain (Chain 2) becomes mathematically stuck in a localized probability valley. Visualizing this mismatch illustrates why tracing diagnostics is critically necessary when upgrading to more complex SEIR modeling architectures.

### Step I: Prior vs. Posterior Adjustments
This section cleanly visualizes how much the real-world data changed our initial parameter assumptions across Transmission, Recovery Rate, Overdispersion, and Detection matrices. Narrow, high-peaked posterior curves rising out of wide, flat priors demonstrate that our case-study data provided overwhelmingly strong statistical signals.

### Step J: Tracking the Effective Reproduction Number ($R(t)$)
Transitioning into the 2020 Switzerland COVID-19 dataset, we map the dynamic $R(t)$ curve across time as it incorporates public health forcing constraints. We distinctly observe the exponential virus spread collapsing immediately below the critical threshold ($R<1$) strictly corresponding with the March 15th government interventions.

### Step K: Final SEIR Posterior Predictive Check
This final plot combines our complex SEIR configurations, governmental policy forcing parameters, and underreporting estimates into one final fit against the real-world Switzerland observations. The close tracking guarantees that the Bayesian framework holds robustly even at the scale of public data.

### Step L: Posterior Statistical Formulations
The text table summarizes our key parameter estimates by providing Posterior Medians alongside 90% Confidence Interval brackets. These quantified final measurements give public policy researchers exact values for key benchmarks such as the viral recovery rate ($\gamma$), daily reporting rate ($p_{reported}$), and structural delay timelines ($\nu$).
