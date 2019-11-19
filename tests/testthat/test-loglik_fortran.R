context("loglik_fortran")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("fortran and r return the same result, but fortran is faster", {

  pars <- c(0.2, 0.1, 1.4, 0.12)
  n_0 <- 2
  cond <- 1
  lx <- 15
  brts <- c(9, 8, 7.6, 7.6, 6, 6, 6, 4, 3, 2, 2)

  t_r <- system.time(
    loglik_r <- mbd::mbd_loglik(
      fortran = FALSE,
      pars = pars,
      brts = brts,
      n_0 = n_0,
      cond = cond,
      q_threshold = 0,
      lx = lx
    )
  )[[3]]
  t_fortran <- system.time(
    loglik_fortran <- mbd::mbd_loglik(
      fortran = TRUE,
      pars = pars,
      brts = brts,
      n_0 = n_0,
      cond = cond,
      q_threshold = 0,
      lx = lx
    )
  )[[3]]
  testthat::expect_equal(loglik_r, loglik_fortran, tolerance = 1e-5)
  testthat::expect_true(t_fortran <= t_r)

})
