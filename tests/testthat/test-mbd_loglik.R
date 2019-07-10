context("mbd_loglik")

# silent if verbose is FALSE ----
test_that("silent if verbose is FALSE", {

  pars <- c(0.2, 0.15, 2, 0.1)
  n_0 <- 2
  age <- 5
  cond <- 0
  brts <- c(3, 2, 1)

  testthat::expect_silent(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0 = n_0,
      cond = cond
    )
  )
})

# likelihoods using expo and lsodes are identical ----
test_that("likelihoods using expo and lsodes are identical", {

  pars <- c(0.23, 0.12, 0.5, 0.24)
  brts <- c(3, 2, 1)
  n_0  <- 2
  cond <- 1
  loglik_expo  <- mbd::mbd_loglik(
    pars = pars, brts = brts, n_0 = n_0, cond = cond, methode = "expo"
  )
  loglik_lsodes <- mbd::mbd_loglik(
    pars = pars, brts = brts, n_0 = n_0, cond = cond, methode = "lsodes"
  )
  testthat::expect_equal(loglik_expo, loglik_lsodes)
})

# nu is zero ----
test_that("nu is zero", {

  # pars[3] is nu
  pars <- c(0.2, 0.1, 0, 0.1)
  brts <- c(3, 2, 1)
  n_0  <- 2
  cond <- 0

  pars2 <- c(0, cond, 0, 0, n_0)
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0 = n_0, # Crown age
      cond = cond
    ),
    DDD::bd_loglik(
      pars1 = pars[1:3],
      pars2 = pars2,
      brts = brts,
      missnumspec = 0
    )
  )
})

# q is zero ----
test_that("q is zero", {

  # pars[4] is q
  pars <- c(0.2, 0.1, 2.0, 0)
  brts <- c(3, 2, 1)
  n_0  <- 2
  cond <- 0

  pars2 <- c(0, cond, 0, 0, n_0)
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0 = n_0,
      cond = cond,
      safety_threshold = 0
    ),
    DDD::bd_loglik(
      pars1 = pars[1:3],
      pars2 = pars2,
      brts = brts,
      missnumspec = 0
    )
  )
})

# nu and q are zero ----
test_that("nu and q are zero", {

  # pars[3] is nu
  # pars[4] is q
  pars <- c(0.2, 0.1, 0.0, 0.0)
  brts <- c(1, 2, 3)
  n_0  <- 2
  cond <- 0

  pars2 <- c(0, cond, 0, 0, n_0)
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0 = n_0,
      cond = cond,
      safety_threshold = 0
    ),
    DDD::bd_loglik(
      pars1 = pars[1:3],
      pars2 = pars2,
      brts = brts,
      missnumspec = 0
    )
  )
})

# abuse ----
test_that("abuse", {

  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.1), # Too few
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    "'pars' must have a length of four"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(NaN, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    "'pars' cannot contain NaNs"
  )
  testthat::expect_true(
    mbd::mbd_loglik(
      pars = c(Inf, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ) == -Inf
  )
  testthat::expect_true(
    mbd::mbd_loglik(
      pars = c(-12.34, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ) == -Inf
  )
  testthat::expect_true(
    mbd::mbd_loglik(
      pars = c(0.2, -12.34, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ) == -Inf
  )
  testthat::expect_true(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, -12.34, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ) == -Inf
  )
  testthat::expect_true(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, -12.34),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ) == -Inf
  )
  testthat::expect_true(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 12.34),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ) == -Inf
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1),
      brts = c(1, 2, 2, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    "At any time you cannot have more speciations than number of species."
  )
})
