context("pmb_loglik")

test_that("use", {

  # pars[2] must be zero
  loglik <- mbd::pmb_loglik(
    pars = c(0.2, 0.0, 2.0, 0.1),
    brts = c(1, 2, 3),
    n_0 = 2 
  )
  testthat::expect_equal(loglik, -3.6017356241900162495)
})

test_that("q is zero", {

  # pars[2] must be zero
  # pars[4] is q
  loglik <- mbd::pmb_loglik(
    pars = c(0.2, 0.0, 2.0, 0.0),
    brts = c(1, 2, 3),
    n_0 = 2 
  )
  testthat::expect_equal(loglik, -3.2271163556401454287)
})

test_that("nu is zero", {

  # pars[2] must be zero
  # pars[3] is nu
  loglik <- mbd::pmb_loglik(
    pars = c(0.2, 0.0, 0.0, 0.1),
    brts = c(1, 2, 3),
    n_0 = 2 
  )
  testthat::expect_equal(loglik, -3.2271163556401454287)
})

test_that("nu and q are zero", {

  # pars[2] must be zero
  # pars[3] is nu
  loglik <- mbd::pmb_loglik(
    pars = c(0.2, 0.0, 0.0, 0.0),
    brts = c(1, 2, 3),
    n_0 = 2 
  )
  testthat::expect_equal(loglik, -3.2271163556401454287)
})

test_that("abuse", {

  lambda <- 0.2
  mu <- 0.0 # Must be zero
  nu <- 2.0
  q <- 0.1
  pars <- c(lambda, mu, nu, q)

  testthat::expect_silent(
    mbd::pmb_loglik(
      pars = pars,
      brts = c(1, 2, 3),
      n_0 = 2 
    )
  )
  testthat::expect_error(
    mbd::pmb_loglik(
      pars = c(NaN, 0.0, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2
    ),
    "'pars' cannot contain NaNs"
  )
  testthat::expect_equal(
    mbd::pmb_loglik(
      pars = c(Inf, 0.0, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2 
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::pmb_loglik(
      pars = c(-123, 0.0, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2 
    ),
    -Inf
  )
  testthat::expect_error(
    mbd::pmb_loglik(
      pars = c(0.2, 12.34, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2 
    ),
    "this function works only for mu = 0!"
  )
  testthat::expect_equal(
    mbd::pmb_loglik(
      pars = c(0.2, 0.0, -12.34, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2 
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::pmb_loglik(
      pars = c(0.2, 0.0, 2.0, -12.34),
      brts = c(1, 2, 3),
      n_0 = 2 
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::pmb_loglik(
      pars = c(0.2, 0.0, 2.0, 12.34),
      brts = c(1, 2, 3),
      n_0 = 2 
    ),
    -Inf
  )
})

test_that("pmb_loglik is called correctly by mbd_loglik", {

  pars <- c(0.2, 0, 1.5, 0.2)
  brts <- c(5, 4, 3, 3, 2, 2)
  n_0   <- 2
  
  testthat::expect_equal(
    pmb <- mbd::pmb_loglik(
      pars = pars,
      brts = brts,
      n_0 = n_0 
    ),
    mbd <- mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0 = n_0
    )
  )
  
})
