context("pmb_loglik")

test_that("use", {

  # pars[2] must be zero
  loglik <- mbd::pmb_loglik(
    pars = c(0.2, 0.0, 2.0, 0.1),
    brts = c(1, 2, 3),
    soc = 2 # Crown age
  )
  expect_true(loglik > -Inf)
})
