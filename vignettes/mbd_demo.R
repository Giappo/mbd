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

