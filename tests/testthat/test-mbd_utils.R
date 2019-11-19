context("mbd_utils")

test_that("get_pkg_name", {
  testthat::expect_true(
    get_pkg_name() == "mbd"
  )
})

test_that("get_function_names", {
  # use
  testthat::expect_true(
    mbd::get_function_names(
      loglik_functions = "mbd_loglik"
    ) == "mbd_loglik"
  )
  testthat::expect_true(
    mbd::get_function_names(
      loglik_functions = c("mbd_loglik")
    ) == "mbd_loglik"
  )
  testthat::expect_true(
    all(
      mbd::get_function_names(
        loglik_functions = c("mbd_loglik", "pmb_loglik")
      ) == c("mbd_loglik", "pmb_loglik")
    )
  )
  testthat::expect_true(
    all(
      mbd::get_function_names(
        loglik_functions = mbd::mbd_logliks_function()
      ) == c("mbd_loglik", "pmb_loglik")
    )
  )
  testthat::expect_true(
    mbd::get_function_names(
      loglik_functions = mbd_loglik
    ) == "mbd_loglik"
  )
  testthat::expect_true(
    all(
      mbd::get_function_names(
        loglik_functions = c(mbd_loglik, pmb_loglik)
      ) == c("mbd_loglik", "pmb_loglik")
    )
  )
  # abuse
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    get_pkg_name(),
    "!"
  )
  testthat::expect_error(
    mbd::get_function_names(
      loglik_functions = "nonsense"
    ),
    error_message
  )
  testthat::expect_error(
    mbd::get_function_names(
      loglik_functions = c("nonsense1", "nonsense2")
    ),
    error_message
  )
  testthat::expect_error(
    mbd::get_function_names(
      loglik_functions = grepl
    ),
    error_message
  )
  testthat::expect_error(
    mbd::get_function_names(
      loglik_functions = c(exp, grepl)
    ),
    error_message
  )
})

test_that("get_model_names", {
  # use
  testthat::expect_silent(
    mbd::get_model_names(
      function_names = mbd_logliks_function(),
      verbose = FALSE
    )
  )
  testthat::expect_output(
    mbd::get_model_names(
      function_names = mbd::mbd_logliks_function(),
      verbose = TRUE
    ),
    "You are using the functions: mbd pmb"
  )
  # abuse
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    get_pkg_name(),
    "!"
  )
  testthat::expect_error(
    mbd::get_model_names(
      function_names = "nonsense"
    ),
    error_message
  )
})

test_that("cat2", {
  testthat::expect_output(
    mbd::cat2(
      message = "test",
      verbose = TRUE
    )
  )
})

test_that("print_info", {
  brts <- c(3, 2, 1)
  n_0 <- 2
  cond <- 1
  testthat::expect_output(
    mbd::print_info(
      brts = brts,
      n_0 = n_0,
      cond = cond,
      verbose = TRUE
    )
  )
  brts <- list(c(3, 2, 1), c(2.5, 1.5, 0.5))
  n_0s <- c(2, 1)
  cond <- 1
  testthat::expect_output(
    mbd::print_info(
      brts = brts,
      n_0 = n_0,
      cond = cond,
      verbose = TRUE
    )
  )
  testthat::expect_silent(
    mbd::print_info(
      brts = brts,
      n_0 = n_0,
      cond = cond,
      verbose = FALSE
    )
  )
})

test_that("mbd_logliks_experiment", {
  loglik_func <- get(
    mbd::mbd_logliks_experiment()
  )
  testthat::expect_true(
    is.numeric(
      loglik_func(
        pars = c(0.2, 0.1, 2, 0.1),
        brts = c(10, 8, 6, 6)
      )
    )
  )
})
