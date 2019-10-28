context("cond_probs")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("dps sum up to zero", {

  q <- 0.2
  lx <- 300
  pp <- matrix(0, ncol = lx, nrow = lx)
  pp[2:10, 2:10] <- 1
  pp <- pp / sum(pp)
  mm <- 2:(lx + 1)
  matrices <- cond_prob_p_matrices(q = q, lx = lx)
  nu_matrix <- matrices$nu_matrix
  pp2 <- matrices$empty_pp
  m1 <- matrices$m1
  m2 <- matrices$m2
  pp2[mm, mm] <- pp

  dp_lambda <-
    cond_prob_dp_lambda(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_mu <-
    cond_prob_dp_mu(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_nu <-
    cond_prob_dp_nu(pp = pp, nu_matrix = nu_matrix)

  dp_nu2 <- cond_prob_dp_nu2(pp = pp, nu_matrix = nu_matrix)

  expect_equal(
    sum(dp_lambda),
    0,
    tolerance = 10 * .Machine$double.eps
  )
  expect_equal(
    sum(dp_mu),
    0,
    tolerance = 10 * .Machine$double.eps
  )
  expect_equal(
    sum(dp_nu),
    0,
    tolerance = 10 * .Machine$double.eps
  )
  expect_equal(
    sum(dp_nu2),
    0,
    tolerance = 10 * .Machine$double.eps
  )
})

test_that("diagonal in p_nu_matrix is (1 - q) ^ m", {
  q <- 0.2
  lx <- 50
  nu_matrix <- cond_prob_p_matrices(q = q, lx = lx)$nu_matrix
  testthat::expect_equal(
    diag(nu_matrix),
    (1 - q) ^ (0:(lx - 1))
  )
})

test_that("p_m1_m2 sums up to one", {
  pars <- c(0.2, 0.1, 2.5, 0.2)
  lx <- 30
  brts <- c(2)
  matrices <- cond_prob_p_matrices(
    q = pars[4],
    lx = lx
  )

  p_m1_m2 <- prob_cond_get_p_m1_m2(
    pars = pars,
    brts = brts,
    matrices = matrices
  )

  testthat::expect_gt(
    sum(p_m1_m2),
    0.99
  )
})

test_that("cond == 0 yields pc == 1; cond == 1 yields pc <= 1", {

  pars <- c(0.2, 0.15, 2.0, 0.1)
  brts <- c(5, 4, 3, 3, 1)
  n_0 <- 2
  lx <- 20

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

test_that("mu = 0", {
  pars <- c(0.2, 0, 1, 0.1)
  brts <- c(3)
  cond <- 1
  n_0 <- 2
  lx <- 30
  test_p <- cond_prob_p(
    pars = pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    lx = lx
  )
  test_q <- cond_prob_q(
    pars = pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    lx = lx
  )
  testthat::expect_equal(
    test_p,
    1,
    tolerance = 1e-3
  )
  testthat::expect_equal(
    test_q,
    1,
    tolerance = 1e-3
  )
})

test_that("bd", {
  pars <- c(0.2, 0.15, 0, 0)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  lx <- 27
  mu_vec <- seq(from = 0.1, to = pars[1], length.out = 2)
  for (m in seq_along(mu_vec)) {
    pars[2] <- mu_vec[m]
    test0 <- exp(
      DDD::bd_loglik(
        pars1 = c(pars[1], pars[2], 0, 0),
        pars2 = c(0, 0, 0, 0, n_0),
        brts = brts,
        missnumspec = 0
      ) -
        DDD::bd_loglik(
          pars1 = c(pars[1], pars[2], 0, 0),
          pars2 = c(0, 1, 0, 0, n_0),
          brts = brts,
          missnumspec = 0
        )
    )
    test_p <- cond_prob_p(
      pars = pars,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      lx = lx
    )
    test_q <- cond_prob_p(
      pars = pars,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      lx = lx
    )
    testthat::expect_equal(
      test0,
      test_p,
      tolerance = 1e-3
    )
    testthat::expect_equal(
      test0,
      test_q,
      tolerance = 1e-3
    )
  }
})

test_that("accurate and fast - gentle parameters", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  pars <- c(0.2, 0.1, 1.0, 0.15)
  brts <- c(10)
  n_0 <- 2
  cond <- 1
  n_sims <- 1e3
  prob_cond_sim <- cond_prob_sim(
    pars = pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    n_sims = n_sims
  )
  lx <- 25
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

  # conditioning time must be smaller than the time for full likelihood
  expect_true(time_cond_p < time_likelihood)
  expect_true(time_cond_q < time_likelihood)

  # conditional likelihood must be "close" to simulations
  expect_equal(prob_cond_p, prob_cond_sim, tolerance = 1e-2)
  expect_equal(prob_cond_q, prob_cond_sim, tolerance = 1e-2)
})

test_that("accurate and fast - harder parameters", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  pars <- c(0.2, 0.1, 1.5, 0.15)
  brts <- c(10)
  n_0 <- 2
  cond <- 1
  n_sims <- 1e4
  prob_cond_sim <- cond_prob_sim(
    pars = pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    n_sims = n_sims
  )
  lx <- 50
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

  # conditional likelihood must be "close" to simulations
  expect_equal(prob_cond_p, prob_cond_sim, tolerance = 1e-3)
  expect_equal(prob_cond_q, prob_cond_sim, tolerance = 1e-3)
})
