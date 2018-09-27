## ------------------------------------------------------------------------
library(mbd)

## ------------------------------------------------------------------------
set.seed(14)

## ------------------------------------------------------------------------
lambda <- 0.2 # sympatric speciation rate
mu <- 0.15 # extinction rate;
nu <- 2.0 # multiple allopatric speciation trigger rate
q <- 0.1 # single-lineage speciation probability
sim_pars <- c(lambda, mu, nu, q)

## ------------------------------------------------------------------------
soc <- 2 # Use a crown age
crown_age <- 1
cond <- 1 # Condition on non-extinction

sim <- mbd_sim(
  pars = sim_pars, 
  soc = soc, 
  age = crown_age, 
  cond = cond, 
  tips_interval = c(0, Inf)
)

## ------------------------------------------------------------------------
sim <- mbd_sim_checked(
  mbd_params = create_mbd_params(lambda = lambda, mu = mu, nu = nu, q = q),
  crown_age = crown_age,
  conditioned_on = "non_extinction"
)

## ---- fig.show='hold', fig.width=7, fig.height=7-------------------------
graphics::plot(sim$tas)

## ---- fig.show='hold', fig.width=7, fig.height=7-------------------------
graphics::plot(sim$tes)

## ------------------------------------------------------------------------
knitr::kable(head(sim$L))

## ------------------------------------------------------------------------
knitr::kable(head(sim$brts))

## ------------------------------------------------------------------------
mbd::mbd_loglik(
  pars = sim_pars, 
  brts = sim$brts, 
  soc = soc, 
  cond = cond, 
  missnumspec = 0
)

## ------------------------------------------------------------------------
phylogeny <- ape::read.tree(text = "((A:1, B:1):2, C:3);")  
ape::plot.phylo(phylogeny)
brts <- ape::branching.times(phylogeny)

## ------------------------------------------------------------------------
lambda <- 0.3 # sympatric speciation rate
mu <- 0.1 # extinction rate
nu <- 0.11 # multiple allopatric speciation trigger rate
q <- 0.15 # single-lineage speciation probability

## ------------------------------------------------------------------------
idparsopt <- 4 # Only optimize the fourth parameter, q
idparsfix <- c(1, 2, 3) # Fix all parameters except q
parsfix <- c(lambda, mu, nu) # Use the known values for the fixed parameters
initparsopt <- q # Set an initial guess for q
out <- mbd_ml(
  brts = brts, 
  initparsopt = initparsopt, 
  idparsopt = idparsopt, 
  parsfix = parsfix, 
  idparsfix = idparsfix, 
  soc = soc, 
  cond = cond
)
knitr::kable(out)

## ------------------------------------------------------------------------
init_param_values <- create_mbd_params(
  lambda = lambda, mu = mu, nu = nu, q = q
)

## ------------------------------------------------------------------------
fixed_params <- create_mbd_params_selector(
  lambda = TRUE, 
  mu = TRUE,
  nu = TRUE
)

## ------------------------------------------------------------------------
estimated_params <- create_mbd_params_selector(
  q = TRUE
)

## ------------------------------------------------------------------------
out <- mbd_calc_max_lik(
  branching_times = brts,
  init_param_values = init_param_values,
  fixed_params = fixed_params,
  estimated_params = estimated_params,
  init_n_species = 2,
  conditioned_on = "non_extinction"
)
knitr::kable(out)

