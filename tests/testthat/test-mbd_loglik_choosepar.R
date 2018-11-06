context("mbd_loglik_choosepar")

test_that("mbd_loglik_choosepar yields a negative loglik", {

  set.seed(10)
  brts <- c(10, 9, 7, 6, 5, 3, 2, 1)
  test_pars <- c(0.3, 0.1, 1, 0.10)
  n_0 <- 2
  cond <- 1

  trparsopt <- c(0.13)
  trparsfix <- c(0.16, 0.66, 0.09)
  idparsopt <- 2
  idparsfix <- (1:4)[-idparsopt]

  out <- mbd_loglik_choosepar(
    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    brts = brts,
    n_0 = n_0,
    cond = cond
  )

  testthat::expect_true(out <= 0)
})

test_that("mbd_loglik_choosepar yields the same result as mbd_loglik", {

  set.seed(10)
  brts <- c(10, 9, 7, 6, 5, 3, 2, 1)
  pars <- c(0.3, 0.1, 1, 0.10)
  n_0 <- 2
  cond <- 1

  trpars <- mbd_transform_forward(pars)
  idparsopt <- 2
  idparsfix <- (1:4)[-idparsopt]
  trparsopt <- trpars[idparsopt]
  trparsfix <- trpars[idparsfix]

  cs_out <- mbd_loglik_choosepar(
    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    brts = brts,
    n_0 = n_0,
    cond = cond
  )
  ll_out <- mbd_loglik(
    pars = pars,
    brts = brts,
    n_0 = n_0,
    cond = cond
  )
  testthat::expect_equal(cs_out, ll_out)
})
