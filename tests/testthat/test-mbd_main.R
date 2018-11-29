context("mbd_main")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("use", {
  seed_interval <- 8:(8 + 3 * is_on_ci())
  for (seed in seed_interval) {
    # seed = 8 is critical
    sim_pars <- c(0.2, 0.15, 1.2, 0.1)
    cond <- 1
    test <- mbd_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      n_0 = 2,
      age = 10,
      tips_interval = c(0, Inf),
      start_pars = c(0.3, 0.1, 1, 0.1),
      loglik_function = mbd_loglik,
      verbose = FALSE
    )
    testthat::expect_true(
      is.data.frame(test)
    )
    testthat::expect_true(
      length(test$sim_lambda) > 0
    )
    testthat::expect_true(
      length(test$sim_mu) > 0
    )
    testthat::expect_true(
      length(test$sim_nu) > 0
    )
    testthat::expect_true(
      length(test$sim_q) > 0
    )
    testthat::expect_true(
      length(test$lambda) > 0
    )
    testthat::expect_true(
      length(test$mu) > 0
    )
    testthat::expect_true(
      length(test$nu) > 0
    )
    testthat::expect_true(
      length(test$q) > 0
    )
    testthat::expect_true(
      length(test$loglik) > 0
    )
    testthat::expect_true(
      length(test$df) > 0
    )
    testthat::expect_true(
      length(test$conv) > 0
    )
    testthat::expect_true(
      length(test$tips) > 0
    )
    testthat::expect_true(
      length(test$seed) > 0
    )
    testthat::expect_true(
      length(test$model) > 0
    )

    # test file saving
    if (.Platform$OS.type == "windows") {
      sim_path <- system.file("extdata", package = mbd_pkg_name())
    } else {
      sim_path <- getwd()
    }
    # check data_path folder existence
    data_path <- file.path(sim_path, "data")
    testthat::expect_true(
      file.exists(data_path)
    )
    # check data file existence
    data_file_name <- file.path(
      data_path,
      paste0(mbd_pkg_name(), "_sim_", seed, ".RData")
    )
    testthat::expect_true(
      file.exists(data_file_name)
    )
    # check results file existence
    results_file_name <- file.path(
      sim_path,
      paste0(mbd_pkg_name(), "_mle_", seed, ".txt")
    )
    testthat::expect_true(
      file.exists(results_file_name)
    )
    # check if saved results are the right ones
    testthat::expect_equal(
      utils::read.csv(results_file_name)[, -1],
      test
    )
  }
})
