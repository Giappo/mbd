## ------------------------------------------------------------------------
# devtools::install_github("Giappo/mbd", quiet = TRUE) # nolint
library(mbd)

## ------------------------------------------------------------------------
set.seed(2)

## ------------------------------------------------------------------------
lambda <- 0.2 # sympatric speciation rate
mu <- 0.15 # extinction rate;
nu <- 2.0 # multiple allopatric speciation trigger rate
q <- 0.1 # single-lineage speciation probability

## ------------------------------------------------------------------------
crown_age <- 1
sim_pars <- c(lambda, mu, nu, q)
sim <- mbd::mbd_sim(
  pars = sim_pars,
  n_0 = 2, # Use a crown age
  age = crown_age,
  cond = 1 # Condition on non-extinction
)

## ---- fig.show='hold', fig.width=7, fig.height=7-------------------------
graphics::plot(sim$full_tree)

## ---- fig.show='hold', fig.width=7, fig.height=7-------------------------
graphics::plot(sim$reconstructed_tree)

## ------------------------------------------------------------------------
knitr::kable(head(sim$l_matrix))

## ------------------------------------------------------------------------
knitr::kable(head(sim$brts))

## ------------------------------------------------------------------------
mbd::mbd_loglik(
  pars = c(lambda, mu, nu, q),
  brts = sim$brts,
  n_0 = 2, # Crown age
  cond = 1  # Non-extinction
)

## ------------------------------------------------------------------------
phylogeny <- ape::read.tree(text = "((A:1, B:1):2, C:3);")
ape::plot.phylo(phylogeny)
brts <- ape::branching.times(phylogeny)

## ------------------------------------------------------------------------
brts <- sim$brts
start_pars <- c(0.2, 0.15, 1, 0.15)
n_0 <- 2
cond <- 1
out <- mbd::mbd_ml(
  start_pars = start_pars,
  brts = brts,
  cond = cond,
  n_0 = n_0,
  optim_ids = c(FALSE, FALSE, FALSE, TRUE),
  verbose = TRUE
)
knitr::kable(out)
