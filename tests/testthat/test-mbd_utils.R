context("mbd_utils")

# get_pkg_name -----
test_that("get_pkg_name", {
  testthat::expect_true(
    get_pkg_name() == "mbd"
  )
})

# get_function_names -----
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

# get_model_names -----
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

# cat2 -----
test_that("cat2", {
  testthat::expect_output(
    mbd::cat2(
      message = "test",
      verbose = TRUE
    )
  )
  testthat::expect_silent(
    mbd::cat2(
      message = "test",
      verbose = FALSE
    )
  )
})

# print_info -----
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

# mbd_logliks_experiment -----
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

# count_percentage_mb_events -----
test_that("count_percentage_mb_events", {

  brts <- c(10, 8, 7, 7, 7, 5, 3, 3, 3, 3)
  percentage <- mbd::count_percentage_mb_events(brts = brts)
  testthat::expect_true(
    percentage >= 0 && percentage <= 1
  )

})

# count_n_mb_events -----
test_that("count_n_mb_events", {

  brts <- c(10, 8, 7, 7, 7, 5, 3, 3, 3, 3)
  number_of_events <- mbd::count_n_mb_events(brts = brts)
  testthat::expect_true(number_of_events >= 0)

})

# file_path -----
test_that("file_path", {

  testthat::expect_equal(
    mbd::file_path("pippo", "/baudo"), # nolint indeed the test I need
    "pippo/baudo"
  )

})
