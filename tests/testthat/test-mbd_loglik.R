context("mbd_loglik")

test_that("silent if verbose is FALSE", {

  skip("fix")
  expect_silent(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    )
  )
})

test_that(
  "number of items to replace is not a multiple of replacement length", {

  skip("number of items to replace is not a multiple of replacement length, #5")
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

test_that("q is zero", {

  skip("Cannot do ML estimation on BD trees, Issue #4, #4")
  # pars[4] is q
  loglik <- mbd::mbd_loglik(
    pars = c(0.2, 0.1, 2.0, 0.0),
    brts = c(1, 2, 3),
    soc = 2 # Crown age
  )
  expect_equal(loglik, -2.6008336107847660479)
})

test_that("nu is zero", {

  # pars[3] is nu
  loglik <- mbd::mbd_loglik(
    pars = c(0.2, 0.1, 0.0, 0.1),
    brts = c(1, 2, 3),
    soc = 2 # Crown age
  )
  expect_equal(loglik, -2.6008336107847660479)
})

test_that("nu and q are zero", {

  # pars[3] is nu
  # pars[4] is q
  loglik <- mbd::mbd_loglik(
    pars = c(0.2, 0.1, 0.0, 0.0),
    brts = c(1, 2, 3),
    soc = 2 # Crown age
  )
  expect_equal(loglik, -2.6008336107847660479)
})




test_that("abuse", {

  expect_error(
    mbd::mbd_loglik(
      pars = c(0.1), # Too few
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'pars' must have a length of four"
  )

  expect_error(
    mbd::mbd_loglik(
      pars = c(NaN, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'pars' cannot contain NaNs"
  )
  expect_error(
    mbd::mbd_loglik(
      pars = c(Inf, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'pars' cannot contain Infs"
  )
  expect_error(
    mbd::mbd_loglik(
      pars = c(-12.34, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'lambda' must be positive"
  )
  expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, -12.34, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'mu' must be positive"
  )
  expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, -12.34, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'nu' must be positive"
  )
  expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, -12.34),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'q' must be positive"
  )
  expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 12.34),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'q' must be less or equal to one"
  )
  expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1,  # Non-extinction
      minimum_multiple_births = -1234
    ),
    "'minimum_multiple_births' must be positive"
  )
})
