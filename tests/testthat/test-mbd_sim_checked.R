context("mbd_sim_checked")

test_that("same as classic interface", {
  
  lambda <- 0.2 # sympatric speciation rate
  mu <- 0.15 # extinction rate;
  nu <- 2.0 # multiple allopatric speciation trigger rate
  q <- 0.1 # single-lineage speciation probability
  sim_pars <- c(lambda, mu, nu, q)
  crown_age <- 1
  
  set.seed(42)
  classic_sim <- mbd_sim(
    pars = sim_pars, 
    soc = 2, 
    age = crown_age, 
    cond = 1
  )

  set.seed(42)
  new_sim <- mbd_sim(
    pars = sim_pars, 
    soc = 2, 
    age = crown_age, 
    cond = 1
  )
  expect_equal(names(classic_sim), names(new_sim))
  expect_equal(classic_sim$btrs, new_sim$btrs)
})
