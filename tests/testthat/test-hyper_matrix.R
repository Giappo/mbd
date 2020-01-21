context("hyper_matrix")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

# big but not too big ----
test_that("big but not too big", {

  m <- mbd::hyper_matrix2(n_species = 2, k = 2, q = 0.1)
  testthat::expect_equal(m[1, 1], 0.81)
  testthat::expect_equal(m[1, 2], 0.00)
  testthat::expect_equal(m[1, 3], 0.00)
  testthat::expect_equal(m[2, 1], 0.36)
  testthat::expect_equal(m[2, 2], 0.729)
  testthat::expect_equal(m[2, 3], 0.00)
  testthat::expect_equal(m[3, 1], 0.04)
  testthat::expect_equal(m[3, 2], 0.405)
})

# is silent ----
test_that("is silent", {

  if (!is_on_ci()) {
    testthat::expect_true(TRUE)
  } else {
    testthat::expect_silent(
      mbd::hyper_matrix2(n_species = 2000, k = 2, q = 0.1)
    )
    gc()
  }
})

# abuse ----
test_that("abuse", {

  if (!is_on_ci()) {
    testthat::expect_true(TRUE)
  } else {
    testthat::expect_error(
      mbd::hyper_matrix2(n_species = 100000, k = 2, q = 0.1),
      "'n_species' must be below 46340"
    )
    gc()
  }
})
