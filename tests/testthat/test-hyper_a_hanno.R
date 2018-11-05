context("hyper_a_hanno")

test_that("big but not too big", {
  
  m <- mbd:::hyper_a_hanno(n_species = 2, k = 2, q = 0.1)
  testthat::expect_equal(m[1, 1], 0.81)
  testthat::expect_equal(m[1, 2], 0.00)
  testthat::expect_equal(m[1, 3], 0.00)
  testthat::expect_equal(m[2, 1], 0.36)
  testthat::expect_equal(m[2, 2], 0.729)
  testthat::expect_equal(m[2, 3], 0.00)
  testthat::expect_equal(m[3, 1], 0.04)
  testthat::expect_equal(m[3, 2], 0.405)
  testthat::expect_equal(m[3, 3], 0.6561)
})

test_that("big but not too big", {
  testthat::expect_silent(
    mbd:::hyper_a_hanno(n_species = 10000, k = 2, q = 0.1)
  )
  gc()
})

test_that("abuse", {
  testthat::expect_error(
    mbd:::hyper_a_hanno(n_species = 100000, k = 2, q = 0.1),
    "'n_species' must be below 46340"
  )
  gc()
})
