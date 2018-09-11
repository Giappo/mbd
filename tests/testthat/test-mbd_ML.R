context("mbd_ML")

test_that("mbd_ML can be silent", {

  skip("Fix @Giappo")
  set.seed(10)
  test_pars <- c(0.3, 0.1, 0.1, 0.15)
  idparsfix <- c(1,2,3)
  
  testthat::expect_silent(
    MBD:::mbd_ML(
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

test_that("abuse", {
  
  testthat::expect_error(
    MBD:::mbd_ML(
      brts = "nonsense", 
      initparsopt = 0.11,
      idparsopt = 4,
      idparsfix = c(1,2,3), parsfix = test_pars[idparsfix], missnumspec = 0, cond = 1, soc = 2
    ),
    "'brts' must be numeric"
  )
})
