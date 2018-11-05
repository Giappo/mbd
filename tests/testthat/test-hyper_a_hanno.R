context("hyper_a_hanno")

# big but not too big ----
test_that("big but not too big", {

  skip("Just testing style now")

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

# is silent ----
test_that("is silent", {

  skip("Just testing style now")

  testthat::expect_silent(
    mbd:::hyper_a_hanno(n_species = 10000, k = 2, q = 0.1)
  )
  gc()
})

# abuse ----
test_that("abuse", {

  skip("Just testing style now")

  testthat::expect_error(
    mbd:::hyper_a_hanno(n_species = 100000, k = 2, q = 0.1),
    "'n_species' must be below 46340"
  )
  gc()
})
