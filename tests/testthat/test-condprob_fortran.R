context("condprob_fortran")

test_that("the right parmsvec is returned", {

  pars <- c(0.2, 0.1, 1.4, 0.12)
  lx_top <- 400
  lx <- 8

  m1_mat <- condprob_m1_mat(lx_top)
  m2_mat <- condprob_m2_mat(m1_mat)
  empty_mat <- condprob_empty_mat(lx_top)

  # test for Q-equation
  eq <- "q_eq"
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  expect_equal(
    exp(log_nu_mat + log_q_mat),
    cond_prob_q_matrices(q = pars[4], lx = lx)$nu_matrix
  )
  parmsvec <- condprob_parmsvec(
    m1_mat = m1_mat,
    m2_mat = m2_mat,
    empty_mat = empty_mat,
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    pars = pars,
    lx = lx,
    eq = eq
  )
  expect_equal(
    length(parmsvec),
    length(pars) - 1 + lx ^ 2 + lx ^ 2 + lx ^ 2 + (lx + 2) ^ 2
  )

  # test for P-equation
  eq <- "p_eq"
  log_nu_mat <- condprob_log_nu_mat(lx = lx, eq = eq)
  log_q_mat <- condprob_log_q_mat(lx = lx, q = pars[4], eq = eq)
  expect_equal(
    exp(log_nu_mat + log_q_mat),
    cond_prob_p_matrices(q = pars[4], lx = lx)$nu_matrix
  )
  parmsvec <- condprob_parmsvec(
    m1_mat = m1_mat,
    m2_mat = m2_mat,
    empty_mat = empty_mat,
    log_nu_mat = log_nu_mat,
    log_q_mat = log_q_mat,
    pars = pars,
    lx = lx,
    eq = eq
  )
  expect_equal(
    length(parmsvec),
    length(pars) - 1 + lx ^ 2 + lx ^ 2 + lx ^ 2 + (lx + 2) ^ 2
  )

})
