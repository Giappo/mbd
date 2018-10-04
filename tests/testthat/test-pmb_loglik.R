context("pmb_loglik")

test_that("use", {

  # pars[2] must be zero
  loglik <- mbd::pmb_loglik(
    pars = c(0.2, 0.0, 2.0, 0.1),
    brts = c(1, 2, 3),
    soc = 2 # Crown age
  )
  expect_equal(loglik, -3.6017356241900162495)
})

test_that("abuse", {

  lambda <- 0.2
  mu <- 0.0 # Must be zero
  nu <- 2.0
  q <- 0.1
  pars <- c(lambda, mu, nu, q)

  expect_silent(
    mbd::pmb_loglik(
      pars = pars,
      brts = c(1, 2, 3),
      soc = 2 # Crown age
    )
  )
  expect_error(
    mbd::pmb_loglik(
      pars = c(NaN, 0.0, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2 # Crown age
    ),
    "'pars' cannot contain NaNs"
  )
  expect_error(
    mbd::pmb_loglik(
      pars = c(Inf, 0.0, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2 # Crown age
    ),
    "'pars' cannot contain Infs"
  )
})
