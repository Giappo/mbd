context("mbd_ml")

test_that("use", {
  set.seed(10)
  brts <- c(10, 9, 7, 6, 5, 3, 2, 1)
  test_pars <- c(0.3, 0.1, 1, 0.10)
  optim_pars <- c(FALSE, TRUE, FALSE, FALSE)
  n_0 <- 2
  cond <- 1

  mbd_out <- mbd::mbd_ml(
    start_pars = c(0.2, 0.15, 1, 0.1),
    true_pars = test_pars,
    optim_pars = optim_pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    verbose = FALSE
  )

  testthat::expect_equal(mbd_out$lambda, test_pars[1])
  testthat::expect_equal(mbd_out$nu, test_pars[3])
  testthat::expect_equal(mbd_out$q, test_pars[4])
  testthat::expect_equal(mbd_out$conv, 0)
  testthat::expect_equal(mbd_out$df, sum(optim_pars))
  testthat::expect_true(mbd_out$loglik <= 0)
})

test_that("compare results from bd and mbd in case of nu = q = 0", {
  set.seed(10)
  brts <- c(10, 9, 7, 6, 5, 3)
  test_pars <- c(0.3, 0.1, 0, 0.10)
  optim_pars <- c(FALSE, TRUE, FALSE, FALSE)
  start_pars <- c(0.2, 0.15, 1, 0.1)
  n_0 <- 2
  cond <- 1
  mbd_out <- mbd::mbd_ml(
    start_pars = start_pars,
    true_pars = test_pars,
    optim_pars = optim_pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    verbose = FALSE
  )

  if (rappdirs::app_dir()$os != "win") {
    sink("/dev/null")
  } else {
    sink(rappdirs::user_cache_dir())
  }
  bd_out <- DDD::bd_ML(
    brts = brts,
    idparsopt = (1:4)[optim_pars],
    idparsfix = (1:4)[!optim_pars],
    initparsopt = start_pars[(1:4)[optim_pars]],
    parsfix = test_pars[(1:4)[!optim_pars]],
    soc = n_0,
    cond = cond
  )
  sink()

  testthat::expect_true(abs(unname(mbd_out$mu) - unname(bd_out$mu0)) < 10 ^
                          -3)
})

test_that("mbd_ml can be silent", {
  set.seed(10)
  brts <- c(10, 9, 7, 6, 5, 3)
  test_pars <- c(0.3, 0.1, 0, 0.10)
  optim_pars <- c(FALSE, TRUE, FALSE, FALSE)
  n_0 <- 2
  cond <- 1
  testthat::expect_silent(
    mbd::mbd_ml(
      start_pars = c(0.2, 0.15, 1, 0.1),
      true_pars = test_pars,
      optim_pars = optim_pars,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    )
  )
})
