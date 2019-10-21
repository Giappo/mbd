context("cond_probs")

# cond == 0 yields pc == 1; cond == 1 yields pc <= 1;
test_that("cond == 0 yields pc == 1; cond == 1 yields pc <= 1", {

  pars <- c(0.2, 0.15, 2.0, 0.1)
  brts <- c(5, 4, 3, 3, 1)
  n_0 <- 2
  lx <- 25

  testthat::expect_equal(
    cond_prob_p(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 0,
      lx = lx
    ),
    1
  )
  testthat::expect_lt(
    cond_prob_p(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 1,
      lx = lx
    ),
    1
  )
  testthat::expect_equal(
    cond_prob_q(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 0,
      lx = lx
    ),
    1
  )
  testthat::expect_lt(
    cond_prob_q(
      pars = pars,
      brts = brts,
      n_0  = n_0,
      cond = 1,
      lx = lx
    ),
    1
  )
})

test_that("accurate and fast", {

  pars <- c(0.2, 0.10, 2.0, 0.15)
  brts <- c(10)
  n_0 <- 2
  prob_cond_sim <- cond_prob_sim(
    pars = pars,
    brts = brts,
    cond = cond,
    n_0 = n_0
  )
  lxs <- seq(from = 25, to = 33, by = 4)
  for (l in seq_along(lxs)) {
    lx <- lxs[l]
    time_cond_p <- system.time(
      prob_cond_p <- cond_prob_p(
        pars = pars,
        brts = brts,
        cond = cond,
        n_0 = n_0,
        lx = lx
      )
    )[[3]]
    time_cond_q <- system.time(
      prob_cond_q <- cond_prob_q(
        pars = pars,
        brts = brts,
        cond = cond,
        n_0 = n_0,
        lx = lx
      )
    )[[3]]
    # unconditioned likelihood time for the same branching times
    time_likelihood <- system.time(
      mbd_loglik(pars = pars, brts = brts, n_0 = n_0, cond = 0, lx = lx ^ 2)
    )[[3]]
    # we want the conditioning time to be smaller than the time
    #  for the full likelihood computation.
    expect_true(time_cond_p < time_likelihood)
    expect_true(time_cond_q < time_likelihood)
  }
})
