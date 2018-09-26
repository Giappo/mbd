## ------------------------------------------------------------------------
library(mbd)

## ------------------------------------------------------------------------
set.seed(14)

## ------------------------------------------------------------------------
lambda <- 0.2 # sympatric speciation rate
mu <- 0.15 # extinction rate;
nu <- 2.0 # multiple allopatric speciation trigger rate
q <- 0.1 # single-lineage speciation probability
sim_pars <- c(lambda, mu, nu, q); 

## ------------------------------------------------------------------------
soc <- 2 # Use a crown age
crown_age <- 10
cond <- 1 # Condition on non-extinction

sim <- mbd_sim(
  pars = sim_pars, 
  soc = soc, 
  age = crown_age, 
  cond = cond, 
  tips_interval = c(0, Inf)
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
sim <- mbd_sim(
  pars = sim_pars, 
  soc = soc, 
  age = crown_age / 10, 
  cond = cond, 
  tips_interval = c(0, Inf)
)
ape::plot.phylo(sim$tes)

## ------------------------------------------------------------------------
idparsopt <- 4 # Only optimize the fourth parameter, q
ids <- 1:4
idparsfix <- ids[-idparsopt] # Fix all parameters except q
parsfix <- sim_pars[idparsfix] # Use the known values for the fixed parameters
initparsopt <- 0.15 # Set an initial guess for q
mbd_ML(
  brts = sim$brts, 
  initparsopt = initparsopt, 
  idparsopt = idparsopt, 
  parsfix = parsfix, 
  idparsfix = idparsfix, 
  soc = soc, 
  cond = cond
)

