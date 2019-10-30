context("condprob_fortran")

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
