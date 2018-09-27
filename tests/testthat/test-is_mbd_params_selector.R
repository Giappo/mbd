context("is_mbd_params_selector")

test_that("use", {
  expect_true(is_mbd_params_selector(create_mbd_params_selector()))
  expect_false(is_mbd_params_selector("nonsense"))
  expect_false(is_mbd_params_selector(NULL))
  expect_false(is_mbd_params_selector(NA))
})
