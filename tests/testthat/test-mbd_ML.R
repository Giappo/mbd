context("mbd_ML")

test_that("mbd_ML can be silent", {

  skip("Fix @Giappo")
  # https://github.com/richelbilderbeek/razzo_article/issues/10
  set.seed(10)
  test_pars <- c(0.3, 0.1, 0.1, 0.15)
  idparsfix <- c(1,2,3)

  expect_silent(
    mbd_ML(
      brts = c(1,2,3),
      initparsopt = 0.11,
      idparsopt = 4,
      idparsfix = idparsfix,
      parsfix = test_pars[idparsfix],
      missnumspec = 0,
      cond = 1,
      soc = 2,
      verbose = FALSE
    )
  )
})

test_that("use", {

  phylogeny <- ape::read.tree(text = "((A:1, B:1):2, C:3);")  
  # maximize the likelihood only for the parameter q
  brts <- ape::branching.times(phylogeny)
  lambda <- 0.1 # speciation rate
  mu <- 0.2 # extinction rate
  nu <- 0.3 # trigger rate
  q <- 0.4 # ?speciation rate during event? 
  init_pars <- c(lambda, mu, nu, q)
  idparsopt <- 4 # Only optimize the fourth parameter, q
  ids <- 1:4
  idparsfix <- ids[-idparsopt] # Fix all parameters except q
  parsfix <- init_pars[idparsfix] # Use the known values for the fixed parameters
  initparsopt <- init_pars[idparsopt] # Set an initial guess for q
  ml_classic <- mbd_ML(
    brts = brts, 
    initparsopt = initparsopt, 
    idparsopt = idparsopt, 
    parsfix = parsfix, 
    idparsfix = idparsfix, 
    soc = 2, # Condition on crown age 
    cond = 1, # Condition on stem/crown age and non-extinction
    verbose = FALSE
  )
  expect_equal(lambda, ml_classic$lambda)
  expect_equal(mu, ml_classic$mu)
  expect_equal(nu, ml_classic$nu)
  expect_equal(0.2201847468814873976, ml_classic$q)
})

test_that("abuse", {

  testthat::expect_error(
    mbd:::mbd_ML(
      brts = "nonsense",
      initparsopt = 0.11,
      idparsopt = 4,
      idparsfix = c(1,2,3), parsfix = test_pars[idparsfix], missnumspec = 0, cond = 1, soc = 2
    ),
    "'brts' must be numeric"
  )
})
