context("nee_functions")

test_that("1 minus functions work", {

  lambda <- 0.25
  mu <- 0.1
  t <- 6
  testthat::expect_equal(
    1 - mbd::p_t(t = t, lambda = lambda, mu = mu),
    mbd::one_minus_pt(t = t, lambda = lambda, mu = mu)
  )
  testthat::expect_equal(
    1 - mbd::ut(t = t, lambda = lambda, mu = mu),
    mbd::one_minus_ut(t = t, lambda = lambda, mu = mu)
  )

  ###

  lambda <- 0.25
  mu <- 0.25
  t <- 7
  testthat::expect_equal(
    1 - mbd::p_t(t = t, lambda = lambda, mu = mu),
    mbd::one_minus_pt(t = t, lambda = lambda, mu = mu)
  )
  testthat::expect_equal(
    1 - mbd::ut(t = t, lambda = lambda, mu = mu),
    mbd::one_minus_ut(t = t, lambda = lambda, mu = mu)
  )

})
