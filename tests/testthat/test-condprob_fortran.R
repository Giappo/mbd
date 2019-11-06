context("condprob_fortran")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("the right parmsvecs and differentials are returned", {

  pars <- c(0.2, 0.1, 1.4, 0.12)
  lx <- 38

  # test for Q-equation
  eq <- "q_eq"
  q_matrices <- cond_prob_q_matrices(q = pars[4], lx = lx)
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  expect_equal(
    exp(log_nu_mat + log_q_mat),
    q_matrices$nu_matrix
  )
  parmsvec <- condprob_parmsvec(
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    pars = pars,
    lx = lx,
    eq = eq
  )
  expect_equal(
    length(parmsvec),
    length(pars) - 1 + lx ^ 2
  )

  # differential
  qq <- matrix(0, lx, lx)
  qq[2:(floor(lx / 2)), 2:(floor(lx / 2))] <- 1; qq <- qq / sum(qq)
  qvec <- matrix(qq, lx ^ 2, 1)
  qvec <- 0:(lx ^ 2 - 1); qvec <- qvec / sum(qvec)

  dq1 <- condprob_dq(qvec = qvec, parmsvec = parmsvec)
  dq2 <- cond_prob_q_rhs1(
    qvec = qvec,
    lambda = pars[1],
    mu = pars[2],
    nu = pars[3],
    nu_matrix = q_matrices$nu_matrix,
    m1 = q_matrices$m1,
    m2 = q_matrices$m2,
    empty_qq = q_matrices$empty_qq,
    t = 0,
    k = 1
  )
  expect_equal(dq1, dq2)

  # test for P-equation
  eq <- "p_eq"
  p_matrices <- cond_prob_p_matrices(q = pars[4], lx = lx)
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  expect_equal(
    exp(log_nu_mat + log_q_mat),
    p_matrices$nu_matrix
  )
  parmsvec <- condprob_parmsvec(
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    pars = pars,
    lx = lx,
    eq = eq
  )
  expect_equal(
    length(parmsvec),
    length(pars) - 1 + lx ^ 2
  )

  # differential
  pp <- matrix(0, lx, lx)
  pp[2:(floor(lx / 2)), 2:(floor(lx / 2))] <- 1; pp <- pp / sum(pp)
  pvec <- matrix(pp, lx ^ 2, 1)
  pvec <- 0:(lx ^ 2 - 1); pvec <- pvec / sum(pvec)

  dp1 <- condprob_dp(pvec = pvec, parmsvec = parmsvec)
  dp2 <- cond_prob_p_rhs1(
    pvec = pvec,
    lambda = pars[1],
    mu = pars[2],
    nu = pars[3],
    nu_matrix = p_matrices$nu_matrix,
    m1 = p_matrices$m1,
    m2 = p_matrices$m2,
    empty_pp = p_matrices$empty_pp,
    t = 0
  )
  expect_equal(dp1, dp2)
})

test_that("pc integrated", {

  pars <- c(0.2, 0.1, 1.4, 0.12)
  lx <- 60
  brts <- c(10); age <- max(brts)

  # P EQUATION
  ## R code
  t_p_r <- system.time(
    p_m1_m2 <- prob_cond_get_p_m1_m2(
      pars = pars,
      brts = brts,
      matrices = cond_prob_p_matrices(q = pars[4], lx = lx),
      rhs_function = cond_prob_p_rhs2
    )
  )[[3]]
  expect_equal(p_m1_m2, t(p_m1_m2))

  ## FORTRAN code
  eq <- "p_eq"
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  parmsvec <- condprob_parmsvec(
    pars = pars,
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    lx = lx,
    eq = eq
  )
  pp <- matrix(0, lx, lx); pp[2, 2] <- 1; pvec <- matrix(pp, lx ^ 2, 1)
  t_p_fortran <- system.time(
    p_fortran <- mbd_solve(
      vector = pvec,
      time_interval = age,
      func = "mbd_runmodpcp",
      parms = parmsvec
    )
  )[[3]]
  dim(p_fortran) <- c(lx, lx)
  expect_equal(p_fortran, t(p_fortran))
  expect_equal(p_m1_m2, p_fortran)

  # Q EQUATION
  ## R code
  t_q_r <- system.time(
    q_m1_m2 <- prob_cond_get_q_m1_m2(
      pars = pars,
      brts = brts,
      matrices = cond_prob_q_matrices(q = pars[4], lx = lx),
      rhs_function = cond_prob_q_rhs2
    )
  )[[3]]
  expect_equal(q_m1_m2, t(q_m1_m2))

  ## FORTRAN code
  eq <- "q_eq"
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  parmsvec <- condprob_parmsvec(
    pars = pars,
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    lx = lx,
    eq = eq
  )
  qq <- matrix(0, lx, lx); qq[1, 1] <- 1; qvec <- matrix(qq, lx ^ 2, 1)
  t_q_fortran <- system.time(
    q_fortran <- mbd_solve(
      vector = qvec,
      time_interval = age,
      func = "mbd_runmodpcq",
      parms = parmsvec
    )
  )[[3]]
  dim(q_fortran) <- c(lx, lx)
  expect_equal(q_fortran, t(q_fortran))
  expect_equal(q_m1_m2, q_fortran)

  skip("Fortran is slower than R!")
  expect_true(t_p_fortran <= t_p_r)
  expect_true(t_q_fortran <= t_q_r)
})

test_that("pc integrated - high 'mbness'", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  pars <- c(0.2, 0.1, 2.5, 0.4)
  lx <- 60
  brts <- c(10); age <- max(brts)

  # P EQUATION
  ## R code
  t_p_r <- system.time(
    p_m1_m2 <- prob_cond_get_p_m1_m2(
      pars = pars,
      brts = brts,
      matrices = cond_prob_p_matrices(q = pars[4], lx = lx),
      rhs_function = cond_prob_p_rhs2
    )
  )[[3]]
  expect_equal(p_m1_m2, t(p_m1_m2))

  ## FORTRAN code
  eq <- "p_eq"
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  parmsvec <- condprob_parmsvec(
    pars = pars,
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    lx = lx,
    eq = eq
  )
  pp <- matrix(0, lx, lx); pp[2, 2] <- 1; pvec <- matrix(pp, lx ^ 2, 1)
  t_p_fortran <- system.time(
    p_fortran <- mbd_solve(
      vector = pvec,
      time_interval = age,
      func = "mbd_runmodpcp",
      parms = parmsvec
    )
  )[[3]]
  dim(p_fortran) <- c(lx, lx)
  expect_equal(p_fortran, t(p_fortran))
  expect_equal(p_m1_m2, p_fortran)

  # Q EQUATION
  ## R code
  t_q_r <- system.time(
    q_m1_m2 <- prob_cond_get_q_m1_m2(
      pars = pars,
      brts = brts,
      matrices = cond_prob_q_matrices(q = pars[4], lx = lx),
      rhs_function = cond_prob_q_rhs2
    )
  )[[3]]
  expect_equal(q_m1_m2, t(q_m1_m2))

  ## FORTRAN code
  eq <- "q_eq"
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  parmsvec <- condprob_parmsvec(
    pars = pars,
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    lx = lx,
    eq = eq
  )
  qq <- matrix(0, lx, lx); qq[1, 1] <- 1; qvec <- matrix(qq, lx ^ 2, 1)
  t_q_fortran <- system.time(
    q_fortran <- mbd_solve(
      vector = qvec,
      time_interval = age,
      func = "mbd_runmodpcq",
      parms = parmsvec
    )
  )[[3]]
  dim(q_fortran) <- c(lx, lx)
  expect_equal(q_fortran, t(q_fortran))
  expect_equal(q_m1_m2, q_fortran)

  skip("Fortran is slower than R!")
  expect_true(t_p_fortran <= t_p_r)
  expect_true(t_q_fortran <= t_q_r)
})
