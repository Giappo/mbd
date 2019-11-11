context("condprob")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

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

test_that("right differentials in R", {

  pars <- c(0.2, 0.1, 1.4, 0.12)
  lx <- 5

  # test for P-equation
  eq <- "p_eq"
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = FALSE)
  pvec <- c(1, rep(0, lx ^ 2 - 1))
  dp <- condprob_dp(
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
  dq <- condprob_dq(
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

test_that("full P_{n1, n2} and Q_{m1, m2} distributions", {

  pars <- c(0.3, 0.15, 1.8, 0.11)
  lx <- 10
  brts <- c(10)

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
  q_r <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = mbd::condprob_dq_rhs
  )
  testthat::expect_equal(q_r, t(q_r))
  ## FORTRAN
  parmsvec <-
    mbd::create_fast_parmsvec(pars = pars, lx = lx, eq = eq, fortran = TRUE)
  q_fortran <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    lx = lx,
    rhs_function = "mbd_runmodpcq"
  )
  testthat::expect_equal(q_fortran, t(q_fortran))
  testthat::expect_equal(q_fortran, q_r)
})
test_that("FORTRAN and R return same result but FORTRAN is faster", {

  brts <- c(10)
  pars <- c(0.2, 0.1, 1.2, 0.12)
  lx <- 25

  # P-equation
  eq <- "p_eq"

  fortran <- TRUE
  t_fortran <- system.time(
    pc_fortran <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  fortran <- FALSE
  t_r <- system.time(
    pc_r <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-5)
  testthat::expect_true(t_fortran <= t_r)

  # Q-equation
  eq <- "q_eq"

  fortran <- TRUE
  t_fortran <- system.time(
    pc_fortran <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  fortran <- FALSE
  t_r <- system.time(
    pc_r <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-5)
  testthat::expect_true(t_fortran <= t_r)
})

test_that("FORTRAN and R - hard test", {

  if (!is_on_ci()) {skip("To be performed on ci.")}

  brts <- c(10)
  pars <- c(0.2, 0.1, 1.5, 0.14)
  lx <- 60

  # P-equation
  eq <- "p_eq"

  fortran <- TRUE
  t_fortran <- system.time(
    pc_fortran <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  fortran <- FALSE
  t_r <- system.time(
    pc_r <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-5)
  testthat::expect_true(t_fortran <= t_r)

  # Q-equation
  eq <- "q_eq"

  fortran <- TRUE
  t_fortran <- system.time(
    pc_fortran <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  fortran <- FALSE
  t_r <- system.time(
    pc_r <- calculate_condprob(
      pars = pars,
      brts = brts,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )[[3]]
  testthat::expect_equal(pc_r, pc_fortran, tolerance = 1e-5)
  testthat::expect_true(t_fortran <= t_r)
})
