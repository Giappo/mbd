context("mbd_main")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

test_that("use", {

  if (!is_on_ci()) {
    skip("This is long")
  }

  seed_interval <- 8:10
  for (seed in seed_interval) {
    # seed = 8 is critical
    sim_pars <- c(0.2, 0.15, 1.2, 0.1)
    cond <- 1
    n_0 <- 2
    t_0s <- age <- 10
    loglik_functions <- mbd_loglik
    optim_ids <- c(TRUE, TRUE, TRUE, TRUE)

    test <- mbd_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      n_0 = n_0,
      age = age,
      tips_interval = c(0, Inf),
      start_pars = c(0.3, 0.1, 1, 0.1),
      loglik_functions = loglik_functions,
      optim_ids = optim_ids,
      verbose = FALSE
    )

    testthat::expect_true(
      is.data.frame(test)
    )
    testthat::expect_true(
      length(test$sim_lambda) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$sim_mu) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$sim_nu) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$sim_q) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$lambda) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$mu) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$nu) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$q) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$loglik) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$df) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$conv) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$tips) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$cond) == length(loglik_functions)
    )
    testthat::expect_true(
      all(
        is.numeric(test$cond)
      )
    )
    testthat::expect_true(
      length(test$seed) == length(loglik_functions)
    )
    testthat::expect_true(
      length(unique(test$seed)) == 1
    )
    testthat::expect_true(
      all(
        test$n_0 == n_0
      )
    )
    testthat::expect_true(
      all(
        test$t_0 == age
      )
    )
    testthat::expect_true(
      length(test$optim_lambda) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$optim_mu) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$optim_nu) == length(loglik_functions)
    )
    testthat::expect_true(
      length(test$optim_q) == length(loglik_functions)
    )
    testthat::expect_true(
      all(
        c(
          test$optim_lambda,
          test$optim_mu,
          test$optim_nu,
          test$optim_q
        ) %in% c("TRUE", "FALSE")
      )
    )
    testthat::expect_true(
      length(test$model) == length(loglik_functions)
    )

    # test file saving
    pkg_name <- get_pkg_name() # nolint internal function
    if (.Platform$OS.type == "windows") {
      project_folder <- system.file("extdata", package = get_pkg_name())
    } else {
      project_folder <- getwd()
    }
    # check data folder existence
    data_folder <- file.path(project_folder, "data")
    testthat::expect_true(
      file.exists(data_folder)
    )
    # check results folder existence
    results_folder <- file.path(project_folder, "results")
    testthat::expect_true(
      file.exists(results_folder)
    )
    # check data file existence
    data_file_name <- create_data_file_name( # nolint internal function
      data_folder = data_folder,
      sim_pars = sim_pars,
      optim_ids = optim_ids,
      cond = cond,
      n_0 = n_0,
      t_0s = t_0s,
      seed = seed
    )
    testthat::expect_true(
      file.exists(data_file_name)
    )
    # check results file existence
    results_file_name <- create_results_file_name( # nolint internal function
      results_folder = results_folder,
      sim_pars = sim_pars,
      optim_ids = optim_ids,
      cond = cond,
      n_0 = n_0,
      t_0s = t_0s,
      seed = seed
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

test_that("it works also for a subset of parameters", {
  seed <- 10
  sim_pars <- c(0.3, 0.2, 1.5, 0.1)
  cond <- 1
  n_0 <- 2
  t_0s <- age <- 6
  loglik_functions <- mbd_loglik
  optim_ids <- c(TRUE, FALSE, FALSE, FALSE)

  test <- mbd_main(
    seed = seed,
    sim_pars = sim_pars,
    cond = cond,
    n_0 = n_0,
    age = age,
    tips_interval = c(0, Inf),
    start_pars = c(0.3, 0.1, 1, 0.1),
    loglik_functions = loglik_functions,
    optim_ids = optim_ids,
    verbose = FALSE
  )

  testthat::expect_true(
    is.data.frame(test)
  )
  testthat::expect_true(
    all(test$sim_mu_m == test$mu_m)
  )
  testthat::expect_true(
    all(test$sim_lambda_s == test$lambda_s)
  )
  testthat::expect_true(
    all(test$sim_mu_s == test$mu_s)
  )
  pkg_name <- get_pkg_name() # nolint internal function
  # test file saving
  if (.Platform$OS.type == "windows") {
    project_folder <- system.file("extdata", package = pkg_name)
  } else {
    project_folder <- getwd()
  }
  # check data folder existence
  data_folder <- file.path(project_folder, "data")
  testthat::expect_true(
    file.exists(data_folder)
  )
  # check results folder existence
  results_folder <- file.path(project_folder, "results")
  testthat::expect_true(
    file.exists(results_folder)
  )
  # check data file existence
  data_file_name <- create_data_file_name( # nolint internal function
    data_folder = data_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  testthat::expect_true(
    file.exists(data_file_name)
  )
  suppressWarnings(file.remove(data_file_name))
  # check results file existence
  results_file_name <- create_results_file_name( # nolint internal function
    results_folder = results_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  testthat::expect_true(
    file.exists(results_file_name)
  )
  suppressWarnings(file.remove(results_file_name))
})

test_that("abuse", {
  seed <- 1
  sim_pars <- c(0.3, 0.2, 2, 0.1)
  cond <- 1
  optim_ids <- c(TRUE, FALSE, FALSE, FALSE)
  n_0 <- 2
  age <- 10
  testthat::expect_error(
    test <- mbd_main(
      seed = seed,
      sim_pars = sim_pars,
      cond = cond,
      n_0 = n_0,
      age = age,
      tips_interval = c(0, Inf),
      start_pars = c(0.3, 0.1, 1, 0.1),
      loglik_functions = "nonsense",
      optim_ids = optim_ids,
      verbose = FALSE
    ),
    paste0(
      "This is not a likelihood function provided by ",
      get_pkg_name(),
      "!"
    )
  )
})
