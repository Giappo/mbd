context("condprob")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

print_from_global <- function(var = "seed") {
  if (var %in% ls(.GlobalEnv)) {
    cat(var, "is", get(var), "\n")
  }
}

# colSums respect constraints ----
test_that("colSums respect constraints", {

  q_vec <- 0.1 * 1:5
  for (q in q_vec) {
    absorb <- TRUE
    lambda <- 0.8
    mu <- 0.2
    nu <- 1

    pars <- c(lambda, mu, nu, q)
    lx <- 8

    eq <- "p_eq"
    nu_q_mat <- mbd::condprob_nu_matrix_p(pars = pars, lx = lx, absorb = absorb)
    testthat::expect_equal(
      unname(colSums(nu_q_mat)),
      rep(1, nrow(nu_q_mat))
    )

    eq <- "q_eq"
    nu_q_mat <- mbd::condprob_nu_matrix_q(pars = pars, lx = lx, absorb = absorb)
    testthat::expect_equal(
      # diag(nu_q_mat),
      # (1 - q) ^ (0:(lx - 1) + 1) # k is equal to 1
      unname(colSums(nu_q_mat)),
      rep(1 + q, nrow(nu_q_mat))
    )
  }

})

# right parmsvec in FORTRAN and R ----
testthat::test_that("right parmsvec in FORTRAN and R", {

  absorb <- TRUE
  pars <- c(0.3, 0.1, 1.7, 0.13)
  lx <- 6
  lx2 <- lx ^ 2

  # FORTRAN
  fortran <- TRUE
  ## P equation
  eq <- "p_eq"
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  testthat::expect_equal(length(parmsvec), 3 + lx2)
  ## Q equation
  eq <- "p_eq"
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  testthat::expect_equal(length(parmsvec), 3 + lx2)

  # R
  fortran <- FALSE
  ## P equation
  eq <- "p_eq"
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  testthat::expect_true(parmsvec$lambda == pars[1])
  testthat::expect_true(parmsvec$mu == pars[2])
  testthat::expect_true(parmsvec$nu == pars[3])
  testthat::expect_equal(dim(parmsvec$nu_matrix), c(lx, lx))
  testthat::expect_equal(parmsvec$m1, t(parmsvec$m2))
  testthat::expect_equal(dim(parmsvec$empty_mat), c(lx + 2, lx + 2))
  ## Q equation
  eq <- "p_eq"
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  testthat::expect_true(parmsvec$lambda == pars[1])
  testthat::expect_true(parmsvec$mu == pars[2])
  testthat::expect_true(parmsvec$nu == pars[3])
  testthat::expect_equal(dim(parmsvec$nu_matrix), c(lx, lx))
  testthat::expect_equal(parmsvec$m1, t(parmsvec$m2))
  testthat::expect_equal(dim(parmsvec$empty_mat), c(lx + 2, lx + 2))

})

# right differentials in R----
testthat::test_that("right differentials in R", {

  absorb <- TRUE
  fortran <- FALSE
  pars <- c(0.2, 0.1, 1.4, 0.12)
  lx <- 40

  # test for the P-equation
  eq <- "p_eq"
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  pp <- matrix(0, lx, lx)
  pp[2, 2] <- 1
  pvec <- matrix(pp, lx ^ 2, 1)
  pp2 <- parmsvec$empty_mat
  mm <- 1 + 1:lx
  pp2[mm, mm] <- pp

  dp_la <- mbd::condprob_dp_absorb_lambda(
    pp = pp, pp2 = pp2, m1 = parmsvec$m1, m2 = parmsvec$m2, mm = mm
  )
  dp_mu <- mbd::condprob_dp_absorb_mu(
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
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
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
testthat::test_that("full P_{n1, n2} and Q_{m1, m2} distributions", {

  absorb <- TRUE
  pars <- c(0.3, 0.15, 1.8, 0.11)
  lx <- 30
  brts <- c(8)

  # P equation
  eq <- "p_eq"
  ## R
  fortran <- FALSE
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  p_r <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = mbd::condprob_dp_absorb_rhs
  )
  testthat::expect_equal(p_r, t(p_r))
  testthat::expect_true(all(p_r <= 1))
  ## FORTRAN
  fortran <- TRUE
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  p_fortran <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = "mbd_runmodpcp_abs"
  )
  testthat::expect_equal(p_fortran, t(p_fortran))
  testthat::expect_equal(p_fortran, p_r, tolerance = 1e-3 * p_r)

  # Q equation
  eq <- "q_eq"
  ## R
  fortran <- FALSE
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  q_r <- mbd::condprob_q_m1_m2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = mbd::condprob_dq_absorb_rhs
  )
  testthat::expect_equal(q_r, t(q_r))
  testthat::expect_true(all(q_r <= 1))
  ## FORTRAN
  fortran <- TRUE
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    eq = eq,
    lx = lx,
    absorb = absorb,
    fortran = fortran
  )
  q_fortran <- mbd::condprob_q_m1_m2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = "mbd_runmodpcq_abs"
  )
  testthat::expect_equal(q_fortran, t(q_fortran))

  # compare the two
  testthat::expect_equal(q_fortran, q_r, tolerance = 1e-3 * q_r)

})

# P_{n1, n2} sums up to one ----
testthat::test_that("P_{n1, n2} sums up to one", {

  absorb <- TRUE
  pars <- c(0.2, 0.1, 2.5, 0.2)
  lx <- 30
  brts <- c(2)
  eq <- "p_eq"
  fortran <- TRUE

  p_n1_n2 <- mbd::condprob_p_n1_n2(
    rhs_function = "mbd_runmodpcp_abs",
    brts = brts,
    lx = lx,
    parmsvec = mbd::condprob_parmsvec(
      pars = pars,
      eq = eq,
      lx = lx,
      absorb = absorb,
      fortran = fortran
    )
  )

  testthat::expect_gt(
    sum(p_n1_n2),
    0.99
  )
  testthat::expect_lte(
    sum(p_n1_n2),
    1
  )

  # no absorb
  absorb <- TRUE
  p_n1_n2 <- mbd::condprob_p_n1_n2(
    rhs_function = "mbd_runmodpcp_abs",
    brts = brts,
    lx = lx,
    parmsvec = mbd::condprob_parmsvec(
      pars = pars,
      eq = eq,
      lx = lx,
      absorb = absorb,
      fortran = fortran
    )
  )

  testthat::expect_gt(
    sum(p_n1_n2),
    0.99
  )
  testthat::expect_lte(
    sum(p_n1_n2),
    1
  )

})

# FORTRAN vs R: same result but FORTRAN is faster ----
testthat::test_that("FORTRAN vs R: same result but FORTRAN is faster", {

  absorb <- TRUE
  brts <- c(5)
  pars <- c(0.2, 0.1, 1.2, 0.12)
  lx <- 22

  # test for the P-equation
  eq <- "p_eq"
  tp_fortran <- system.time(
    pc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE,
      absorb = absorb
    )
  )[[3]]
  tp_r <- system.time(
    pc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE,
      absorb = absorb
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-3)
  testthat::expect_true(tp_fortran <= tp_r)

  # test for the Q-equation
  eq <- "q_eq"
  tq_fortran <- system.time(
    qc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE,
      absorb = absorb
    )
  )[[3]]
  tq_r <- system.time(
    qc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE,
      absorb = absorb
    )
  )[[3]]
  testthat::expect_equal(qc_r, qc_fortran, tolerance = 1e-3)
  testthat::expect_true(tq_fortran <= tq_r)

})

# FORTRAN vs R: hard test ----
test_that("FORTRAN vs R: hard test", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  absorb <- TRUE
  brts <- c(7)
  pars <- c(0.2, 0.15, 1.5, 0.10)
  lx <- 38

  # test for the P-equation
  eq <- "p_eq"
  tp_fortran <- system.time(
    pc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE,
      absorb = absorb
    )
  )[[3]]
  tp_r <- system.time(
    pc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE,
      absorb = absorb
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-3)
  testthat::expect_true(tp_fortran <= tp_r)

  # test for the Q-equation
  eq <- "q_eq"
  tq_fortran <- system.time(
    qc_fortran <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = TRUE,
      absorb = absorb
    )
  )[[3]]
  tq_r <- system.time(
    qc_r <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = FALSE,
      absorb = absorb
    )
  )[[3]]
  testthat::expect_equal(qc_r, qc_fortran, tolerance = 1e-3)
  testthat::expect_true(tq_fortran <= tq_r)

})

# condprob for mu = 0 ----
test_that("condprob for mu = 0", {

  absorb <- TRUE
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
      fortran = TRUE,
      absorb = absorb
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
      fortran = TRUE,
      absorb = absorb
    ),
    1,
    tolerance = 1e-3
  )
})

# bd: nu = q = 0 ----
test_that("nu = q = 0", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  absorb <- TRUE
  pars <- c(0.2, 0.15, 0, 0)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  lx <- 15
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 2)
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
          fortran = fortran,
          absorb = absorb
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

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  absorb <- TRUE
  pars <- c(0.2, 0.15, 0, 0.5)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  lx <- 15
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 2)
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
          fortran = fortran,
          absorb = absorb
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

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  absorb <- TRUE
  pars <- c(0.2, 0.15, 3, 0)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  lx <- 15
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 2)
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
          fortran = fortran,
          absorb = absorb
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

# probcond vs probcond_sim ----
test_that("probcond vs probcond_sim", {

  skip("Use Nee instead")

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  absorb <- TRUE
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
      fortran = fortran,
      absorb = absorb
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
