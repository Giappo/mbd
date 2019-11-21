context("condprob")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

# diagonal in p_nu_matrix is (1 - q) ^ m (+ k) ----
test_that("diagonal in p_nu_matrix is (1 - q) ^ m (+ k)", {

  q <- 0.2
  lx <- 8

  eq <- "p_eq"
  log_nu_mat <- mbd::condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- mbd::condprob_log_q_mat(lx = lx, eq = eq, q = q)
  nu_q_mat <- exp(log_nu_mat + log_q_mat)
  testthat::expect_equal(
    diag(nu_q_mat),
    (1 - q) ^ (0:(lx - 1))
  )

  eq <- "q_eq"
  log_nu_mat <- mbd::condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- mbd::condprob_log_q_mat(lx = lx, eq = eq, q = q)
  nu_q_mat <- exp(log_nu_mat + log_q_mat)
  testthat::expect_equal(
    diag(nu_q_mat),
    (1 - q) ^ (0:(lx - 1) + 1) # k is equal to 1
  )

})

# right parmsvec in FORTRAN and R ----
test_that("right parmsvec in FORTRAN and R", {

  pars <- c(0.3, 0.1, 1.7, 0.13)
  lx <- 6
  lx2 <- lx ^ 2

  # FORTRAN
  fortran <- TRUE
  ## P equation
  eq <- "p_eq"
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = fortran)
  testthat::expect_equal(length(parmsvec), 3 + lx2)
  ## Q equation
  eq <- "p_eq"
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = fortran)
  testthat::expect_equal(length(parmsvec), 3 + lx2)

  # R
  fortran <- FALSE
  ## P equation
  eq <- "p_eq"
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = fortran)
  testthat::expect_true(parmsvec$lambda == pars[1])
  testthat::expect_true(parmsvec$mu == pars[2])
  testthat::expect_true(parmsvec$nu == pars[3])
  testthat::expect_equal(dim(parmsvec$nu_matrix), c(lx, lx))
  testthat::expect_equal(parmsvec$m1, t(parmsvec$m2))
  testthat::expect_equal(dim(parmsvec$empty_mat), c(lx + 2, lx + 2))
  ## Q equation
  eq <- "p_eq"
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = fortran)
  testthat::expect_true(parmsvec$lambda == pars[1])
  testthat::expect_true(parmsvec$mu == pars[2])
  testthat::expect_true(parmsvec$nu == pars[3])
  testthat::expect_equal(dim(parmsvec$nu_matrix), c(lx, lx))
  testthat::expect_equal(parmsvec$m1, t(parmsvec$m2))
  testthat::expect_equal(dim(parmsvec$empty_mat), c(lx + 2, lx + 2))

})

# right differentials in R ----
test_that("right differentials in R", {

  pars <- c(0.2, 0.1, 1.4, 0.12)
  lx <- 40

  # test for the P-equation
  eq <- "p_eq"
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = FALSE)
  pp <- matrix(0, lx, lx)
  pp[2, 2] <- 1
  pvec <- matrix(pp, lx ^ 2, 1)
  pp2 <- parmsvec$empty_mat
  mm <- 1 + 1:lx
  pp2[mm, mm] <- pp

  dp_la <- mbd::condprob_dp_lambda(
    pp = pp, pp2 = pp2, m1 = parmsvec$m1, m2 = parmsvec$m2, mm = mm
  )
  dp_mu <- mbd::condprob_dp_mu(
    pp = pp, pp2 = pp2, m1 = parmsvec$m1, m2 = parmsvec$m2, mm = mm
  )
  dp_nu <- mbd::condprob_dp_nu(pp = pp, nu_matrix = parmsvec$nu_matrix)
  testthat::expect_equal(sum(dp_la), 0, tolerance = 10 * .Machine$double.eps)
  testthat::expect_equal(sum(dp_mu), 0, tolerance = 10 * .Machine$double.eps)
  testthat::expect_equal(sum(dp_nu), 0, tolerance = 10 * .Machine$double.eps)
  dp <- mbd::condprob_dp(
    pvec = pvec,
    lambda = parmsvec$lambda,
    mu = parmsvec$mu,
    nu = parmsvec$nu,
    nu_matrix = parmsvec$nu_matrix,
    m1 = parmsvec$m1,
    m2 = parmsvec$m2,
    empty_mat = parmsvec$empty_mat
  )
  testthat::expect_equal(length(dp), lx ^ 2)
  testthat::expect_equal(sum(dp), 0, tolerance = 1e-5 * max(dp))

  # test for Q-equation
  eq <- "q_eq"
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = FALSE)
  qvec <- c(1, rep(0, lx ^ 2 - 1))
  dq <- mbd::condprob_dq(
    qvec = qvec,
    lambda = parmsvec$lambda,
    mu = parmsvec$mu,
    nu = parmsvec$nu,
    nu_matrix = parmsvec$nu_matrix,
    m1 = parmsvec$m1,
    m2 = parmsvec$m2,
    empty_mat = parmsvec$empty_mat
  )
  testthat::expect_equal(length(dq), lx ^ 2)

})

# full P_{n1, n2} and Q_{m1, m2} distributions ----
test_that("full P_{n1, n2} and Q_{m1, m2} distributions", {

  pars <- c(0.3, 0.15, 1.8, 0.11)
  lx <- 18
  brts <- c(8)

  # P equation
  eq <- "p_eq"
  ## R
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = FALSE)
  p_r <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = mbd::condprob_dp_rhs
  )
  testthat::expect_equal(p_r, t(p_r))
  ## FORTRAN
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = TRUE)
  p_fortran <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = "mbd_runmodpcp"
  )
  testthat::expect_equal(p_fortran, t(p_fortran))
  testthat::expect_equal(p_fortran, p_r)

  # P equation
  eq <- "q_eq"
  ## R
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = FALSE)
  q_r <- mbd::condprob_q_m1_m2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = mbd::condprob_dq_rhs
  )
  testthat::expect_equal(q_r, t(q_r))
  ## FORTRAN
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = TRUE)
  q_fortran <- mbd::condprob_q_m1_m2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = "mbd_runmodpcq"
  )
  testthat::expect_equal(q_fortran, t(q_fortran))

  # compare the two
  testthat::expect_equal(q_fortran, q_r)

})

# P_{n1, n2} sums up to one ----
test_that("P_{n1, n2} sums up to one", {

  pars <- c(0.2, 0.1, 2.5, 0.2)
  lx <- 22
  brts <- c(2)
  eq <- "p_eq"

  p_n1_n2 <- mbd::condprob_p_n1_n2(
    rhs_function = "mbd_runmodpcp",
    brts = brts,
    lx = lx,
    parmsvec = mbd::create_fast_parmsvec(
      pars = pars,
      lx = lx,
      eq = eq,
      fortran = TRUE
    )
  )

  testthat::expect_gt(
    sum(p_n1_n2),
    0.99
  )

})

# FORTRAN vs R: same result but FORTRAN is faster ----
test_that("FORTRAN vs R: same result but FORTRAN is faster", {

  brts <- c(8)
  pars <- c(0.2, 0.1, 1.2, 0.12)
  lx <- 20

  # test for the P-equation
  eq <- "p_eq"
  tp_fortran <- system.time(
    pc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE
    )
  )[[3]]
  tp_r <- system.time(
    pc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-5)
  testthat::expect_true(tp_fortran <= tp_r)

  # test for the Q-equation
  eq <- "q_eq"
  tq_fortran <- system.time(
    qc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE
    )
  )[[3]]
  tq_r <- system.time(
    qc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE
    )
  )[[3]]
  testthat::expect_equal(qc_r, qc_fortran, tolerance = 1e-5)
  testthat::expect_true(tq_fortran <= tq_r)

})

# FORTRAN vs R: hard test ----
test_that("FORTRAN vs R: hard test", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  brts <- c(10)
  pars <- c(0.2, 0.1, 1.5, 0.14)
  lx <- 50

  # test for the P-equation
  eq <- "p_eq"
  tp_fortran <- system.time(
    pc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE
    )
  )[[3]]
  tp_r <- system.time(
    pc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-5)
  testthat::expect_true(tp_fortran <= tp_r)

  # test for the Q-equation
  eq <- "q_eq"
  tq_fortran <- system.time(
    qc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE
    )
  )[[3]]
  tq_r <- system.time(
    qc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE
    )
  )[[3]]
  testthat::expect_equal(qc_r, qc_fortran, tolerance = 1e-5)
  testthat::expect_true(tq_fortran <= tq_r)

})

# condprob for mu = 0 ----
test_that("condprob for mu = 0", {

  pars <- c(0.2, 0, 1, 0.1)
  brts <- c(3)
  cond <- 1
  n_0 <- 2
  lx <- 25
  testthat::expect_equal(
    mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = "p_eq",
      fortran = TRUE
    ),
    1,
    tolerance = 1e-3
  )
  testthat::expect_equal(
    mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = "q_eq",
      fortran = TRUE
    ),
    1,
    tolerance = 1e-3
  )
})

# bd: nu = q = 0 ----
test_that("nu = q = 0", {

  pars <- c(0.2, 0.15, 0, 0)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  lx <- 15
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 2 + is_on_ci())
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
    for (eq in c("p_eq", "q_eq")) {
      for (fortran in c(TRUE, FALSE)) {
        test <- mbd::calculate_condprob(
          pars = pars,
          brts = brts,
          lx = lx,
          eq = eq,
          fortran = fortran
        )
        testthat::expect_equal(
          test0,
          test,
          tolerance = 1e-3
        )
      }
    }
  }

})

# bd: nu = 0 ----
test_that("nu = 0", {

  pars <- c(0.2, 0.15, 0, 0.5)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  lx <- 15
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 2 + is_on_ci())
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
    for (eq in c("p_eq", "q_eq")) {
      for (fortran in c(TRUE, FALSE)) {
        test <- mbd::calculate_condprob(
          pars = pars,
          brts = brts,
          lx = lx,
          eq = eq,
          fortran = fortran
        )
        testthat::expect_equal(
          test0,
          test,
          tolerance = 1e-3
        )
      }
    }
  }

})

# bd: q = 0 ----
test_that("q = 0", {

  pars <- c(0.2, 0.15, 3, 0)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  lx <- 15
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 2 + is_on_ci())
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
    for (eq in c("p_eq", "q_eq")) {
      for (fortran in c(TRUE, FALSE)) {
        test <- mbd::calculate_condprob(
          pars = pars,
          brts = brts,
          lx = lx,
          eq = eq,
          fortran = fortran
        )
        testthat::expect_equal(
          test0,
          test,
          tolerance = 1e-3
        )
      }
    }
  }

})

# condprob_select_eq ----
test_that("condprob_select_eq", {

  max_seed <- 2 + is_on_ci() * 23
  for (seed in 1:max_seed) {
    set.seed(seed)
    lambda <- runif(n = 1, min = 0.05, max = 0.3)
    mu <- runif(n = 1, min = 0, max = lambda)
    nu <- runif(n = 1, min = 0.3, max = 2.3)
    q <- runif(n = 1, min = 0.05, max = 0.35)
    pars <- c(lambda, mu, nu, q)
    brts <- c(8)

    pc_sim <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = 1e2,
      eq = "sim"
    )
    if (pc_sim > 0.45 && pc_sim < 0.55) {
      pc_sim <- mbd::calculate_condprob(
        pars = pars,
        brts = brts,
        lx = 1e3,
        eq = "sim"
      )
    }
    if (pc_sim >= 0.47 && pc_sim <= 0.53) {
      pc_sim <- mbd::calculate_condprob(
        pars = pars,
        brts = brts,
        lx = 1e4,
        eq = "sim"
      )
    }
    if (pc_sim >= 0.495 && pc_sim <= 0.51) {
      pc_sim <- mbd::calculate_condprob(
        pars = pars,
        brts = brts,
        lx = 1e5,
        eq = "sim"
      )
    }
    if (pc_sim > 0.5) {
      right_eq <- "p_eq"
    } else {
      right_eq <- "q_eq"
    }
    t_select <- system.time(
      eq <- mbd::condprob_select_eq(pars = pars, brts = brts, fortran = TRUE)
    )[[3]]
    testthat::expect_equal(eq, right_eq)
    testthat::expect_lt(t_select, 300) # select in less than 5 mins
  }

})

# probcond_select_eq: nasty case ----
test_that("probcond_select_eq: nasty case", {

  pars <- c(1.50, 0.15, 1.35, 0.09)
  brts <- c(5, 2.53, 0.79)
  fortran <- TRUE
  lx <- 23

  pc_sim <- mbd::condprob_sim(
    brts = brts,
    parmsvec = mbd::create_fast_parmsvec(
      pars = pars,
      lx = lx,
      eq = "sim",
      fortran = TRUE
    ),
    lx = 50,
    saveit = TRUE,
    starting_seed = 1
  )
  if (pc_sim > 0.5) {
    right_eq <- "p_eq"
  } else {
    right_eq <- "q_eq"
  }

  eq <- mbd::condprob_select_eq(
    pars = pars,
    brts = brts,
    fortran = fortran
  )
  testthat::expect_equal(eq, right_eq)

})

# probcond_select_eq: more nasty cases ----
test_that("probcond_select_eq: more nasty cases", {

  max_seed <- 2 + is_on_ci() * 23
  for (seed in 1:max_seed) {
    if (seed == 15 || seed == 27) next # too slow for simulations
    set.seed(seed)
    lambda <- runif(n = 1, min = 0.05, max = 1.5)
    mu <- runif(n = 1, min = 0, max = lambda)
    nu <- runif(n = 1, min = 0.3, max = 2.9)
    q <- runif(n = 1, min = 0.05, max = 0.4)
    pars <- c(lambda, mu, nu, q)
    brts <- c(5)

    pc_sim <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = 200,
      eq = "sim"
    )
    if (pc_sim > 0.48 && pc_sim < 0.52) {
      pc_sim <- mbd::calculate_condprob(
        pars = pars,
        brts = brts,
        lx = 1e3,
        eq = "sim"
      )
    }
    if (pc_sim > 0.49 && pc_sim < 0.51) {
      pc_sim <- mbd::calculate_condprob(
        pars = pars,
        brts = brts,
        lx = 1e4,
        eq = "sim"
      )
    }
    if (pc_sim > 0.5) {
      right_eq <- "p_eq"
    } else {
      right_eq <- "q_eq"
    }
    t_select <- system.time(
      eq <- mbd::condprob_select_eq(pars = pars, brts = brts, fortran = TRUE)
    )[[3]]
    testthat::expect_equal(eq, right_eq)
    testthat::expect_lt(t_select, 300) # select in less than 5 mins
  }

})

# probcond vs probcond_sim ----
test_that("probcond vs probcond_sim", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  pars <- c(0.2, 0.1, 1.0, 0.15)
  age <- 10
  brts <- c(age)

  n_sims <- 1e5
  pc_sim <- mbd::calculate_condprob(
    pars = pars,
    brts = brts,
    lx = n_sims,
    eq = "sim"
  )

  fortran <- TRUE
  time_pc <- system.time(
    pc <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = 70,
      fortran = fortran
    )
  )[[3]]

  # conditional likelihood must be "close" to simulations
  testthat::expect_equal(
    pc,
    pc_sim,
    tolerance = 1e-2
  )

  # conditioning time must be less than the time for full likelihood
  testthat::expect_true(time_pc < 300)

})
