context("mbd_ml")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("use", {

  # if !is_on_ci() {
  #   skip "This is long"
  # }

  brts <- c(10, 9, 7, 6, 5)
  start_pars <- c(0.2, 0.15, 1, 0.1)
  n_0 <- 2
  cond <- 1

  test <- mbd::mbd_ml(
    start_pars = start_pars,
    brts = brts,
    cond = cond,
    n_0 = n_0,
    verbose = FALSE,
    lx = 100,
    methode = "lsodes"
  )

  for (var_name in get_param_names()) { # nolint internal function
    testthat::expect_true(test[var_name] >= 0)
  }
  testthat::expect_true(test$q <= 1)
  testthat::expect_true(is.numeric(test$loglik) == TRUE)
  testthat::expect_true(test$loglik <= 0)
  testthat::expect_true(test$df == length(start_pars))
  testthat::expect_true(test$conv == 0)
})

test_that("compare results from bd and mbd in case of nu = q = 0", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

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

  x <- capture.output(bd_out <- DDD::bd_ML(
    brts = brts,
    idparsopt = 1:2,
    initparsopt = start_pars[1:2],
    soc = n_0,
    cond = cond
  ))
  rm(x)

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

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

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

test_that("mbd_ml can produce output", {

  if (!is_on_ci()) {
    skip("To be performed on ci.")
  }

  brts <- c(10, 9, 7, 6, 5, 3)
  optim_ids <- c(FALSE, TRUE, FALSE, FALSE)
  n_0 <- 2
  cond <- 1
  output <- capture.output(
    mbd::mbd_ml(
      start_pars = c(0.2, 0.15, 1, 0.1),
      optim_ids = optim_ids,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      verbose = TRUE
    )
  )
  testthat::expect_true(
    length(output) > 0
  )
  testthat::expect_false(
    is.null(output)
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
    "You cannot start from negative parameters!"
  )
  testthat::expect_output(
    suppressWarnings(mbd::mbd_ml(
      start_pars = c(60, 50, 100, 0.5),
      brts = brts,
      cond = cond,
      n_0 = n_0,
      verbose = FALSE
    )),
    "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values." # nolint
  )
  testthat::expect_error(
    test <- mbd_ml(
      loglik_function = mbd_loglik,
      brts = brts,
      start_pars = c(-1, 0.1, 2, 0.1),
      cond = cond,
      n_0 = 2,
      verbose = FALSE
    ),
    "You cannot start from negative parameters!"
  )
})
