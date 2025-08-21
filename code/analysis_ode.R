library(brms)


# Model -------------------------------------------------------------------

# custom Stan code to deal with the dynamic simulation
# this will be fed into the brms model
# This represents the Zeide function z mentioned in Equations 1 and 2 of
# the main text
Zeide <- "
// Specify dynamical system (ODEs)
// Note: This will need to change in order to include mortality predictions too
vector ode_lnD_zeide(
  real t,
  vector lnD,
  real log_a,
  real log_b,
  real log_c
){
  vector[1] dlnDdt;
  dlnDdt[1] = exp(log_a + expm1(log_b) * lnD[1] - exp(log_c) * expm1(lnD[1]));
  return dlnDdt;
}
// Integrate ODE and prepare output
// Note: This code is written to be observation by observation
real predict_zeide(
  real D0,
  real delta_t,
  real log_a,
  real log_b,
  real log_c
){
  vector[1] y0 = [ log(D0) ]';
  real times[1];
  array[1] vector[1] lnD1;
  times[1] = delta_t;
  lnD1 = ode_rk45(
    ode_lnD_zeide,
    y0,
    0,
    times,
    log_a,
    log_b,
    log_c
  );
  vector[1] y1 = [ lnD1[1,1] ]';
  return(y1[1]);
}
"

# brms model formula
# D1 is final diameter
# D0 is initial diameter
# delta_t is time lapsed
# loga, logb and logc are the Zeide growth parameters
# AccName is species names in our use case
# eta is a placeholder for the linear predictor in brms syntax
frml <- bf(
    D1 ~ eta,
    nlf(eta ~ predict_zeide(D0, delta_t, loga, logb, logc)),
    loga ~ 0 + AccName,
    logb ~ 0 + AccName,
    logc ~ 0 + AccName,
    sigma ~ 0 + AccName,
    nl = TRUE
)

# priors
mypriors <- c(
    prior(normal(0, 0.5), nlpar = "loga"),
    prior(normal(-1, 0.5), nlpar = "logb"),
    prior(normal(-1, 0.5), nlpar = "logc")
)

# mcmc sampling
# we used multithreading which might be a better idea for a high-performance cluster
# you'll need to change the threads and cores settings depending on your computer
mod <- brm(
    frml,
    prior = mypriors,
    stanvars = stanvar(scode = Zeide, block = "functions"),
    data = dat_train,  # replace with your own data
    family = brmsfamily("lognormal", link_sigma = "log"),
    control = list(adapt_delta = 0.99),
    chains = 4,
    backend = "cmdstanr",
    threads = threading(24),
    cores = 4,
    init = 0
)


# AFTER RUNNING YOUR MODEL
# we have to recompile the model with rstan to be able to generate predictions
# (I don't understand it either...)
# modRstan <- brm(
#     frml,
#     prior = mypriors,
#     stanvars = stanvar(scode = Zeide, block = "functions"),
#     data = dat_train,     # replace with your own data
#     backend = "rstan",
#     chains=0,
#     family = brmsfamily("lognormal", link_sigma = "log")
# )