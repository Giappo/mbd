context("pars_transform")

test_that("inputs greater than 1 make no sense", {
  transformed_pars <- c(1.5, 0.15, 0.5, 0.15)
  pars <- pars_transform_back(transformed_pars) # nolint internal function
  testthat::expect_true(
    any(pars < 0)
  )
})

test_that("double transform brings back the arguments", {
  test_pars <- c(0.2, 0.15, 0.5, 0.15)

  testthat::expect_true(all.equal(
    pars_transform_back(pars_transform_forward(test_pars)), # nolint internal function
    test_pars
  ))
  testthat::expect_true(all.equal(
    pars_transform_forward(pars_transform_back(test_pars)), # nolint internal function
    test_pars
  ))

  test_pars <- c(0.2, 0.15, 0, 0.15)
  testthat::expect_true(all.equal(
    pars_transform_back(pars_transform_forward(test_pars)), # nolint internal function
    test_pars
  ))
  testthat::expect_true(all.equal(
    pars_transform_forward(pars_transform_back(test_pars)), # nolint internal function
    test_pars
  ))

  test_pars <- c(0.2, 0.15, Inf, 0.15)
  testthat::expect_true(all.equal(
    pars_transform_back(pars_transform_forward(test_pars)), # nolint internal function
    test_pars
  ))
})
