context("is_mbd_params_selector")

test_that("use", {
  expect_true(is_mbd_params_selector(x = create_mbd_params_selector()))
  expect_false(is_mbd_params_selector(x = "nonsense"))
  expect_false(is_mbd_params_selector(NULL))
  expect_false(is_mbd_params_selector(NA))
  expect_false(
    mbd:::is_mbd_params_selector(
      list(
        lambda = "nonsense",
        mu = TRUE,
        nu = TRUE,
        q = TRUE
      )
    )
  )
  expect_false(
    mbd:::is_mbd_params_selector(
      list(
        lambda = TRUE,
        mu = "nonsense",
        nu = TRUE,
        q = TRUE
      )
    )
  )
  expect_false(
    mbd:::is_mbd_params_selector(
      list(
        lambda = TRUE,
        mu = TRUE,
        nu = "nonsense",
        q = TRUE
      )
    )
  )
  expect_false(
    mbd:::is_mbd_params_selector(
      list(
        lambda = TRUE,
        mu = TRUE,
        nu = TRUE,
        q = "nonsense"
      )
    )
  )
})
