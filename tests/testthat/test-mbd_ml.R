context("mbd_ml")

test_that("use", {

  set.seed(10)
  test_pars <- c(0.3, 0.1, 0.1, 0.15)
  idparsfix <- c(1, 2, 3)

  out <- mbd_ml(
    brts = c(1, 2, 3),
    initparsopt = 0.11,
    idparsopt = 4,
    idparsfix = idparsfix,
    parsfix = test_pars[idparsfix],
    missnumspec = 0,
    cond = 1,
    soc = 2,
    verbose = FALSE
  )
  expect_equal(0.3, out$lambda)
  expect_equal(0.1, out$mu)
  expect_equal(0.1, out$nu)
  expect_equal(0.0010006043479098515256, out$q)
})

test_that("abuse", {

  testthat::expect_error(
    mbd:::mbd_ml(
      brts = "nonsense",
      initparsopt = 0.11,
      idparsopt = 4,
      idparsfix = c(1, 2, 3),
      parsfix = test_pars[idparsfix],
      missnumspec = 0, cond = 1, soc = 2
    ),
    "'brts' must be numeric"
  )
})
