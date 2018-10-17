context("get_mbd_param_names")

test_that("use", {
  created <- get_mbd_param_names()
  expected <- c("lambda", "mu", "nu", "q")
  expect_equal(created, expected)
})
