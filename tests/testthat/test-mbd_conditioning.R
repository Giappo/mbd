context("mbd_conditioning")

# cond == 0 yields pc == 1; cond == 1 yields pc <= 1;
test_that("cond == 0 yields pc == 1; cond == 1 yields pc <= 1", {

  pars <- c(0.2, 0.15, 2.0, 0.1)
  brts <- c(5, 4, 3, 3, 1)
  n_0 <- 2

  testthat::expect_equal(
    cond_prob(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 0
    ),
    1
  )
  testthat::expect_lt(
    cond_prob(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 1,
      lx = 16
    ),
    1
  )
})

# conditioned likelihood is greater than unconditioned likelihood ----
test_that("conditioned likelihood is greater than unconditioned likelihood", {

  pars <- c(0.2, 0.15, 2.0, 0.1)
  brts <- c(5, 4, 3, 3, 1)
  n_0 <- 2

  testthat::expect_true(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 0
    ) <=
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 1
    )
  )
})

# right conditioning ----
test_that("right conditioning", {

  pars <- c(0.2, 0.15, 0, 0.1)
  brts <- c(5, 4, 1)
  n_0  <- 2
  cond <- 1

  pars2 <- c(0, cond, 0, 0, n_0)
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0 = n_0,
      cond = cond,
      lx = 300
    ),
    DDD::bd_loglik(
      pars1 = pars[1:3],
      pars2 = pars2,
      brts = brts,
      missnumspec = 0
    ),
    tolerance = 1e-5
  )
})

test_that("accurate and fast", {

  skip("Issue #77: https://github.com/Giappo/mbd/issues/77")

  lxs <- seq(from = 50, to = 500, by = 50)
  for (l in seq_along(lxs)) {
    lx <- lxs[l]
    time_cond <- system.time(prob_cond <- cond_prob(
      pars = pars,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      lx = lx
    ))[[3]]
    # unconditioned likelihood time for the same branching times
    time_likelihood <- system.time(
      mbd_loglik(pars = pars, brts = brts, n_0 = n_0, cond = 0, lx = lx)
    )
    if (time_cond < time_likelihood && (1 - prob_cond) <= 0.05) {
      break
    }
  }
  # we want the conditioning time to be smaller than the time
  #  for the full likelihood computation.
  expect_true(time_cond < time_likelihood)
  # we want the result to be accurate enough to have an error smaller than 5%
  expect_true((1 - prob_cond) <= 0.05)
})

# abuse ----
test_that("abuse", {

  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.1), # Too few
      brts = c(1, 2, 3),
      n_0 = 2, # Crown age
      cond = 1  # Non extinction of the phylogeny
    ),
    "'pars' must have a length of four"
  )

  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(NaN, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2, # Crown age
      cond = 1  # Non extinction of the phylogeny
    ),
    "'pars' cannot contain NaNs"
  )
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = c(Inf, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = c(-12.34, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, -12.34, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, -12.34, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, -12.34),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 12.34),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
})
