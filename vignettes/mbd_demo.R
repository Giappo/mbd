## ------------------------------------------------------------------------
# devtools::install_github("Giappo/mbd", quiet = TRUE)
library(mbd)

## ------------------------------------------------------------------------
set.seed(2)

## ------------------------------------------------------------------------
lambda <- 0.2 # sympatric speciation rate
mu <- 0.15 # extinction rate;
nu <- 2.0 # multiple allopatric speciation trigger rate
q <- 0.1 # single-lineage speciation probability

## ------------------------------------------------------------------------
if (1 == 2) {
  crown_age <- 1
  sim_pars <- c(lambda, mu, nu, q)
  sim <- mbd_sim(
    pars = sim_pars, 
    n_0 = 2, # Use a crown age 
    age = crown_age,
    cond = 1 # Condition on non-extinction
  )
}

## ---- fig.show='hold', fig.width=7, fig.height=7-------------------------
if (1 == 2) {
  graphics::plot(sim$full_tree)
}

## ---- fig.show='hold', fig.width=7, fig.height=7-------------------------
if (1 == 2) {
  graphics::plot(sim$reconstructed_tree)
}

## ------------------------------------------------------------------------
if (1 == 2) {
  knitr::kable(head(sim$l_matrix))
}

## ------------------------------------------------------------------------
if (1 == 2) {
  knitr::kable(head(sim$brts))
}

## ------------------------------------------------------------------------
if (1 == 2) {
  mbd::mbd_loglik(
    pars = c(lambda, mu, nu, q), 
    brts = sim$brts, 
    n_0 = 2, # Crown age 
    cond = 1  # Non-extinction 
  )
}

## ------------------------------------------------------------------------
phylogeny <- ape::read.tree(text = "((A:1, B:1):2, C:3);")  
ape::plot.phylo(phylogeny)
brts <- ape::branching.times(phylogeny)

## ------------------------------------------------------------------------
if (1 == 2) {
  brts <- sim$brts
  start_pars <- c(0.2, 0.15, 1, 0.15)
  optim_ids <- c(FALSE, FALSE, FALSE, TRUE)
  n_0 <- 2
  cond <- 1
  out <- mbd::mbd_ml(
    start_pars = start_pars,
    true_pars = sim_pars,
    optim_ids = optim_ids,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    verbose = TRUE
  )
  knitr::kable(out)
}

