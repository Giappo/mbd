context("mbd_checks")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

# check_n_0 ----
test_that("check_n_0", {
  testthat::expect_error(
    mbd::check_n_0(n_0 = "pippo"),
    "'n_0' must be numeric."
  )
  testthat::expect_error(
    mbd::check_n_0(n_0 = 3),
    "'n_0' must be either 1 or 2."
  )
})

# check_seed ----
test_that("check_seed", {
  testthat::expect_error(
    mbd::check_seed(seed = "pippo"),
    "'seed' must be integer or NA."
  )
  testthat::expect_error(
    mbd::check_seed(seed = 0.75),
    "'seed' must be integer or NA."
  )
})

# check_pars ----
test_that("check_pars", {
  testthat::expect_error(
    mbd::check_pars(pars = c(2, 1)),
    "'pars' must have a length of four."
  )
  testthat::expect_error(
    mbd::check_pars(pars = c(NA, NaN, NA, NA)),
    "'pars' cannot contain NaNs"
  )
  testthat::expect_equal(
    mbd::check_pars(pars = c(-1, 2, 3, 5)),
    "wrong"
  )
})

# check_brts ----
test_that("check_brts", {
  testthat::expect_error(
    mbd::check_brts(brts = c(10, 10, 8), n_0 = 2),
    "Crown/stem age has to be reported only once in the branching times."
  )
  testthat::expect_error(
    mbd::check_brts(brts = c(10, 8, 8, 8, 8, 4), n_0 = 2),
    "At any time you cannot have more speciations than number of species."
  )
  testthat::expect_error(
    mbd::check_brts(brts = c(10, -8, 4), n_0 = 2),
    "'brts' values must be positive. Present time is taken at the 0."
  )
})

# check_cond ----
test_that("check_cond", {
  testthat::expect_error(
    mbd::check_cond(cond = 40, tips_interval = c(0, Inf), n_0 = 2),
    "This conditioning is not implemented."
  )
  testthat::expect_error(
    mbd::check_cond(cond = 1, tips_interval = c(5, 3), n_0 = 2),
    "'tips_interval' must contain two values, of which the second is larger."
  )
  testthat::expect_error(
    mbd::check_cond(cond = 1, tips_interval = c(2, -3), n_0 = 2)
  )
  testthat::expect_error(
    mbd::check_cond(cond = 0, tips_interval = c(2, 30), n_0 = 2),
    "If 'cond' == 0, you cannot put restrictions on tree tips."
  )
  testthat::expect_error(
    mbd::check_cond(cond = 1, tips_interval = c(1, 30), n_0 = 2),
    "If 'cond' == 1, you cannot require less than 'n_0' tips."
  )
})
