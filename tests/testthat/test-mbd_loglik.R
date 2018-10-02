context("mbd_loglik")

test_that("silent if verbose is FALSE", {

  skip("fix")
  expect_silent(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1), 
      brts = c(1, 2, 3, 4), 
      soc = 2, # Crown age 
      cond = 1  # Non-extinction
    )
  )
})

test_that("uninformative error message", {

  skip("fix")
  expect_silent(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1), 
      brts = c(1, 2, 2, 2, 3), 
      soc = 2, # Crown age 
      cond = 1  # Non-extinction
    )
  )
})

test_that("likelihoods using expo and lsoda are identical", {
  soc <- 2
  cond <- 1
  lambda <- 0.12
  mu <- 0.23
  nu <- 0.34
  q <- 0.45
  simpars <- c(lambda, mu, nu, q)
  brts <- c(1, 2, 2)
  loglik_expo  <- mbd::mbd_loglik(
    pars = simpars, brts = brts, soc = soc, cond = cond, methode = "expo"
  )
  loglik_lsoda <- mbd::mbd_loglik(
    pars = simpars, brts = brts, soc = soc, cond = cond, methode = "lsoda"
  )
  expect_equal(loglik_expo, loglik_lsoda)
})

test_that("Convergence fail?", {
  lambda <- 0.1
  mu <- 0.2
  nu <- 0.3
  q <- 0.4
  pars <- c(lambda, mu, nu, q)
  crown <- 2
  # Here it works
  expect_silent(
    mbd::mbd_loglik(
      pars = pars, 
      brts = c(1, 2, 2), # OK
      soc = crown
    )
  )
  # Is this really what is intended?
  expect_error(
    mbd::mbd_loglik(
      pars = pars, 
      brts = c(1, 2, 3), # Not OK
      soc = crown
    ),
    "Maximum number of species exceeds R's memory limits"
  )
})
