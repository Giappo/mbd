context("mbd_transform")

test_that("inputs greater than 1 make no sense", {
  transformed_pars <- c(1.5, 0.15, 0.5, 0.15)
  testthat::expect_error(
    mbd:::mbd_transform_back(transformed_pars),
    "inputs greater than 1 make no sense"
  )
})

test_that("double transform brings back the arguments", {
  test_pars <- c(0.2, 0.15, 0.5, 0.15)

  testthat::expect_true(all.equal(
    mbd:::mbd_transform_back(mbd:::mbd_transform_forward(test_pars)),
    test_pars
  ))
  testthat::expect_true(all.equal(
    mbd:::mbd_transform_forward(mbd:::mbd_transform_back(test_pars)),
    test_pars
  ))

  test_pars <- c(0.2, 0.15, 0, 0.15)
  testthat::expect_true(all.equal(
    mbd:::mbd_transform_back(mbd:::mbd_transform_forward(test_pars)),
    test_pars
  ))
  testthat::expect_true(all.equal(
    mbd:::mbd_transform_forward(mbd:::mbd_transform_back(test_pars)),
    test_pars
  ))

  test_pars <- c(0.2, 0.15, Inf, 0.15)
  testthat::expect_true(all.equal(
    mbd:::mbd_transform_back(mbd:::mbd_transform_forward(test_pars)),
    test_pars
  ))
})
