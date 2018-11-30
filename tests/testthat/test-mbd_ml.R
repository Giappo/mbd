context("mbd_ml")

test_that("use", {
  brts <- c(10, 9, 7, 6, 5)
  start_pars <- c(0.2, 0.15, 1, 0.1)
  n_0 <- 2
  cond <- 1

  test <- mbd::mbd_ml(
    start_pars = start_pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    verbose = FALSE
  )

  for (var_name in get_mbd_param_names()) {
    testthat::expect_true(test[var_name] >= 0)
  }
  testthat::expect_true(test$q <= 1)
  testthat::expect_true(is.numeric(test$loglik) == TRUE)
  testthat::expect_true(test$loglik <= 0)
  testthat::expect_true(test$df == length(start_pars))
  testthat::expect_true(test$conv == 0)
})

test_that("compare results from bd and mbd in case of nu = q = 0", {
  brts <- c(10, 9, 7, 6, 5, 4, 3, 2, 1)
  start_pars <- c(0.2, 0.15, 0, 0)
  n_0 <- 2
  cond <- 1
  optim_ids <- c(TRUE, TRUE, FALSE, FALSE)

  mbd_out <- mbd::mbd_ml(
    start_pars = start_pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    optim_ids = optim_ids,
    verbose = FALSE
  )

  if (rappdirs::app_dir()$os != "win") {
    sink("/dev/null")
  } else {
    sink(rappdirs::user_cache_dir())
  }
  bd_out <- DDD::bd_ML(
    brts = brts,
    idparsopt = 1:2,
    initparsopt = start_pars[1:2],
    soc = n_0,
    cond = cond
  )
  sink()

  testthat::expect_true(
    abs(
      unname(mbd_out$lambda) - unname(bd_out$lambda0)
    ) < 1e-3
  )
  testthat::expect_true(
    abs(
      unname(mbd_out$mu) - unname(bd_out$mu0)
    ) < 1e-3
  )
})

test_that("mbd_ml can be silent", {
  set.seed(10)
  brts <- c(10, 9, 7, 6, 5, 3)
  optim_ids <- c(FALSE, TRUE, FALSE, FALSE)
  n_0 <- 2
  cond <- 1
  testthat::expect_silent(
    mbd::mbd_ml(
      start_pars = c(0.2, 0.15, 1, 0.1),
      optim_ids = optim_ids,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    )
  )
})

test_that("abuse", {
  brts <- c(10, 9, 7, 6, 5)
  start_pars <- c(0.2, 0.15, 1, 0.1)
  n_0 <- 2
  cond <- 1

  testthat::expect_error(
    mbd::mbd_ml(
      start_pars = c(0.2, 0.15, -1, 0.1),
      brts = brts,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    ),
    "you cannot start from negative parameters"
  )

  testthat::expect_error(
    mbd::mbd_ml(
      start_pars = start_pars,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE,
      true_pars = c(start_pars[1:3], 0.2)
    ),
    "for fixed parameters start from the true values"
  )
})
