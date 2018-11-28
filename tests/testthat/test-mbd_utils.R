context("mbd_utils")

test_that("mbd_pkg_name", {
  testthat::expect_true(
    mbd_pkg_name() == "mbd"
  )
})
