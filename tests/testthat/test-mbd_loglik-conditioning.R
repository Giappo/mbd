context("mbd_loglik - conditioning")

# conditioned likelihood is lesser than unconditioned likelihood ----
test_that("conditioned likelihood is lesser than unconditioned likelihood", {

  pars <- c(0.2, 0.15, 2.0, 0.1)
  brts <- c(5, 4, 3, 3, 1)
  n_0   <- 2
  
  testthat::expect_lt(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0  = n_0, 
      cond = 0 
    ),
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      n_0  = n_0, 
      cond = 1 
    )
  )
})

# right conditioning ----
test_that("right conditioning", {
  
  pars <- c(0.2, 0.15, 0, 0.1)
  brts <- c(5, 4, 1)
  n_0  <- 2
  cond <- 1
  
  pars2 <- c(0, cond, 0, 0, n_0)
  testthat::expect_true(
    abs(
      mbd::mbd_loglik(
        pars = pars,
        brts = brts,
        n_0 = n_0, # Crown age
        cond = cond,
        lx = 300
      ) -
      DDD::bd_loglik(
        pars1 = pars[1:3], 
        pars2 = pars2, 
        brts = brts, 
        missnumspec = 0
      )
    ) <= 1e-5
  )
  
  pars2 <- c(0, cond, 0, 0, n_0)
  testthat::expect_true(
    abs(
      mbd::mbd_loglik(
        pars = pars,
        brts = brts,
        n_0 = n_0, # Crown age
        cond = cond,
        lx = 300
      ) -
      DDD::bd_loglik(
        pars1 = pars[1:3], 
        pars2 = pars2, 
        brts = brts, 
        missnumspec = 0
      )
    ) <= 1e-5
  )
})

# abuse ----
test_that("abuse", {

  expect_error(
    mbd::mbd_loglik(
      pars = c(0.1), # Too few
      brts = c(1, 2, 3),
      n_0 = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'pars' must have a length of four"
  )

  expect_error(
    mbd::mbd_loglik(
      pars = c(NaN, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'pars' cannot contain NaNs"
  )
  expect_equal(
    mbd::mbd_loglik(
      pars = c(Inf, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2, 
      cond = 1  
    ),
    -Inf
  )
  expect_equal(
    mbd::mbd_loglik(
      pars = c(-12.34, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, -12.34, 2.0, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, -12.34, 0.1),
      brts = c(1, 2, 3),
      n_0 = 2,
      cond = 1
    ),
    -Inf
  )
  expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, -12.34),
      brts = c(1, 2, 3),
      n_0 = 2, 
      cond = 1
    ),
    -Inf
  )
  expect_equal(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 12.34),
      brts = c(1, 2, 3),
      n_0 = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    -Inf
  )
  
})
