context("mbd_utils")

test_that("get_pkg_name", {
  testthat::expect_true(
    get_pkg_name() == "mbd"
  )
})

test_that("get_function_names & get_model_names", {
  #use
  testthat::expect_true(
    all(
      get_function_names(
        models = mbd_logliks_function()
      ) == c("mbd_loglik", "pmb_loglik")
    )
  )
  testthat::expect_true(
    get_function_names(
      models = mbd_loglik
    ) == "mbd_loglik"
  )
  testthat::expect_silent(
    get_model_names(
      function_names = mbd_logliks_function(),
      verbose = FALSE
    )
  )
  testthat::expect_output(
    get_model_names(
      function_names = mbd_logliks_function(),
      verbose = TRUE
    ),
    "You are using the functions: mbd pmb"
  )
  #abuse
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    get_pkg_name(),
    "!"
  )
  testthat::expect_error(
    get_function_names(
      models = "nonsense"
    ),
    error_message
  )
  testthat::expect_error(
    get_function_names(
      models = c("nonsense1", "nonsense2")
    ),
    error_message
  )
  testthat::expect_error(
    get_function_names(
      models = "grepl"
    ),
    error_message
  )
  testthat::expect_error(
    get_function_names(
      models = c(exp, grepl)
    ),
    error_message
  )
  testthat::expect_error(
    get_model_names(
      function_names = "nonsense"
    ),
    error_message
  )
})
