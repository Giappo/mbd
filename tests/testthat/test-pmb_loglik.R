context("pmb_loglik")

test_that("use", {

  expect_silent(
    mbd::pmb_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2 # Crown age
    )
  )
})

test_that("Convergence fail", {
  lambda <- 0.1
  mu <- 0.2
  nu <- 0.3
  q <- 0.4
  pars <- c(lambda, mu, nu, q)
  crown <- 2
  # Here it works
  loglik <- mbd::pmb_loglik(
    pars = pars,
    brts = c(1, 2, 2),
    soc = crown
  )
  expect_equal(loglik, -Inf)
})
