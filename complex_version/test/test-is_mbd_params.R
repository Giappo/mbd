context("is_mbd_params")

test_that("use", {
  expect_true(is_mbd_params(x = create_mbd_params(0.1, 0.2, 0.3, 0.4)))
  expect_false(is_mbd_params("nonsense"))
  expect_false(is_mbd_params(NULL))
  expect_false(is_mbd_params(NA))
})
