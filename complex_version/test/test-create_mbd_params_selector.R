context("create_mbd_params_selector")

test_that("use, by default all are false", {
  selector <- create_mbd_params_selector()
  expect_false(selector$lambda)
  expect_false(selector$mu)
  expect_false(selector$nu)
  expect_false(selector$q)
})

test_that("use, all true", {
  selector <- create_mbd_params_selector(
    lambda = TRUE,
    mu = TRUE,
    nu = TRUE,
    q = TRUE
  )
  expect_true(selector$lambda)
  expect_true(selector$mu)
  expect_true(selector$nu)
  expect_true(selector$q)
})

test_that("abuse", {
  expect_error(
    create_mbd_params_selector(
      lambda = "nonsense",
      mu = TRUE,
      nu = TRUE,
      q = TRUE
    ),
    "'lambda' must be either TRUE or FALSE"
  )
  expect_error(
    create_mbd_params_selector(
      lambda = TRUE,
      mu = "nonsense",
      nu = TRUE,
      q = TRUE
    ),
    "'mu' must be either TRUE or FALSE"
  )
  expect_error(
    create_mbd_params_selector(
      lambda = TRUE,
      mu = TRUE,
      nu = "nonsense",
      q = TRUE
    ),
    "'nu' must be either TRUE or FALSE"
  )
  expect_error(
    create_mbd_params_selector(
      lambda = TRUE,
      mu = TRUE,
      nu = TRUE,
      q = "nonsense"
    ),
    "'q' must be either TRUE or FALSE"
  )
})
