functions {
  // SIR ODE system: S = Susceptible, I = Infected, R = Recovered
  // beta: contact/transmission rate, gamma: recovery rate
  // N: total population (passed as integer data, not a parameter)
  real[] sir(real t,
             real[] y,
             real[] theta,
             real[] x_r,
             int[]  x_i) {

    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = x_i[1];

    real beta  = theta[1];
    real gamma = theta[2];

    // dS/dt = -beta * S * I / N
    real dS_dt = -beta * I * S / N;
    // dR/dt =  gamma * I
    real dR_dt =  gamma * I;
    // dI/dt derived from conservation (S+I+R = const)
    real dI_dt = -(dS_dt + dR_dt);

    return {dS_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days;     // number of observation days
  real y0[3];              // initial conditions [S0, I0, R0]
  real t0;                 // initial time
  real ts[n_days];         // observation times
  int N;                   // total population
  int cases[n_days];       // observed in-bed counts
}

transformed data {
  real x_r[0];             // no real fixed data needed by ODE
  int  x_i[1] = { N };    // pass N as integer fixed data
}

parameters {
  real<lower=0> gamma;     // recovery rate
  real<lower=0> beta;      // transmission rate
  real<lower=0> phi_inv;   // inverse overdispersion (phi = 1/phi_inv)
}

transformed parameters {
  // ODE solution: y[day, compartment]
  real y[n_days, 3];
  real phi = 1.0 / phi_inv;

  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;
    // Solve SIR ODE at each observation time
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}

model {
  // Priors (truncated at 0 by <lower=0> in parameters block)
  beta    ~ normal(2,   1);
  gamma   ~ normal(0.4, 0.5);
  phi_inv ~ exponential(5);

  // Negative Binomial likelihood:
  // col(to_matrix(y), 2) extracts the I compartment for all days
  cases ~ neg_binomial_2(col(to_matrix(y), 2), phi);
}

generated quantities {
  // Derived epidemiological quantities
  real R0            = beta / gamma;
  real recovery_time = 1.0 / gamma;

  // Posterior predictive samples (for posterior predictive check plot)
  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y), 2), phi);
}
