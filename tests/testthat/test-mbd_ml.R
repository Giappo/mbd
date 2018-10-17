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
  expect_equal(1.2312452488755977663e-13, out$q)
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
  testthat::expect_error(
    mbd:::mbd_ml(
      brts = c(1, 2, 3),
      initparsopt = c(1, 2, 3, 4, 5, 6, 7, 8), # Too long
      idparsopt = 4,
      idparsfix = c(1, 2, 3),
      parsfix = test_pars[idparsfix],
      missnumspec = 0, cond = 1, soc = 2
    ),
    "lengths of 'idparsopt' and'initparsopt' must match"
  )
  testthat::expect_error(
    mbd:::mbd_ml(
      brts = c(1, 2, 3),
      initparsopt = 0.11,
      idparsopt = 4,
      idparsfix = c(1, 2, 3),
      parsfix = c(1, 2, 3, 4, 5, 6, 7, 8), # Too long
      missnumspec = 0, cond = 1, soc = 2
    ),
    "lengths of 'idparsfix' and'parsfix' must match"
  )
  testthat::expect_error(
    mbd:::mbd_ml(
      brts = c(1, 2, 3),
      initparsopt = 0.11,
      idparsopt = 4,
      idparsfix = c(1, 2, 3),
      parsfix = c(1, 2, 3, 4, 5, 6, 7, 8), # Too long
      missnumspec = 0, cond = 1, soc = 2
    ),
    "lengths of 'idparsfix' and'parsfix' must match"
  )
  testthat::expect_error(
    mbd:::mbd_ml(
      brts = c(1, 2, 3),
      initparsopt = 0.11,
      idparsopt = 4, # 4 present first time
      idparsfix = c(1, 2, 3, 4), # 4 present second time
      parsfix = c(1, 2, 3, 0.1),
      missnumspec = 0, cond = 1, soc = 2
    ),
    paste0(
      "IDs 1 to 4 must be present exactly once ",
      "in either 'idparsfix' or 'idparsopt'"
    )
  )

})
