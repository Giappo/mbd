context("mbd_loglik")

# silent if verbose is FALSE ----
test_that("silent if verbose is FALSE", {

  testthat::skip("fix")
  testthat::expect_silent(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, 
      cond = 1 
    )
  )
})

# mbd_loglik yields the same results as mbd_loglik_basic for non pathologic inputs ----
test_that("mbd_loglik yields the same results as mbd_loglik_basic for non pathologic inputs", {

  pars <- c(0.2, 0.15, 2.0, 0.1)
  brts <- c(5, 4, 3, 3, 1)
  soc  <- 2
  cond <- 0
  precision <- 800
  
  testthat::expect_true(
    abs(
      mbd::mbd_loglik(
        pars = pars,
        brts = brts,
        soc  = soc,
        cond = cond, 
        lx0 = precision 
      ) -
      mbd::mbd_loglik_basic(
        pars = pars,
        brts = brts,
        soc  = soc, 
        cond = cond, 
        lx = precision
      )
    ) <= precision^-1
  )
  
  cond <- 1
  testthat::expect_true(
    abs(
      mbd::mbd_loglik(
        pars = pars,
        brts = brts,
        soc  = soc,
        cond = cond, 
        lx0 = precision 
      ) - 
        mbd::mbd_loglik_basic(
          pars = pars,
          brts = brts,
          soc  = soc, 
          cond = cond, 
          lx = precision
        ) 
    ) <= precision^-1
  )
})

# likelihoods using expo and lsoda are identical ----
test_that("likelihoods using expo and lsoda are identical", {
  
  pars <- c(0.23, 0.12, 0.5, 0.24)
  brts <- c(1, 2, 2)
  soc  <- 2
  cond <- 1
  loglik_expo  <- mbd::mbd_loglik(
    pars = pars, brts = brts, soc = soc, cond = cond, methode = "expo"
  )
  loglik_lsoda <- mbd::mbd_loglik(
    pars = pars, brts = brts, soc = soc, cond = cond, methode = "lsoda"
  )
  testthat::expect_equal(loglik_expo, loglik_lsoda)
})

# nu is zero ----
test_that("nu is zero", {

  # pars[3] is nu
  pars <- c(0.2, 0.1, 0, 0.1)
  brts <- c(1, 2, 3)
  soc  <- 2
  cond <- 0
  
  pars2 <- c(0, cond, 0, 0, soc)
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      soc = soc, # Crown age
      cond = cond
    ),
    DDD::bd_loglik(
      pars1 = pars[1:3], 
      pars2 = pars2, 
      brts = brts, 
      missnumspec = 0
    )
  )
})

# q is zero ----
test_that("q is zero", {
  
  # pars[4] is q
  pars <- c(0.2, 0.1, 2.0, 0)
  brts <- c(1, 2, 3)
  soc  <- 2
  cond <- 0
  
  pars2 <- c(0, cond, 0, 0, soc)
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      soc = soc, # Crown age
      cond = cond
    ),
    DDD::bd_loglik(
      pars1 = pars[1:3], 
      pars2 = pars2, 
      brts = brts, 
      missnumspec = 0
    )
  )
})

# nu and q are zero ----
test_that("nu and q are zero", {

  # pars[3] is nu
  # pars[4] is q
  pars <- c(0.2, 0.1, 0.0, 0.0)
  brts <- c(1, 2, 3)
  soc  <- 2
  cond <- 0
  
  pars2 <- c(0, cond, 0, 0, soc)
  testthat::expect_equal(
    mbd::mbd_loglik(
      pars = pars,
      brts = brts,
      soc = soc, # Crown age
      cond = cond
    ),
    DDD::bd_loglik(
      pars1 = pars[1:3], 
      pars2 = pars2, 
      brts = brts, 
      missnumspec = 0
    )
  )
})

# abuse ----
test_that("abuse", {

  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.1), # Too few
      brts = c(1, 2, 3),
      soc = 2, 
      cond = 1 
    ),
    "'pars' must have a length of four"
  )

  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(NaN, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'pars' cannot contain NaNs"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(Inf, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'pars' cannot contain Infs"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(-12.34, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'lambda' must be positive"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, -12.34, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'mu' must be positive"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, -12.34, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'nu' must be positive"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, -12.34),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'q' must be positive"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 12.34),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "'q' must be less or equal to one"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1),
      brts = c(1, 2, 3),
      soc = 2, # Crown age
      cond = 1,  # Non-extinction
      minimum_multiple_births = -1234
    ),
    "'minimum_multiple_births' must be positive"
  )
  testthat::expect_error(
    mbd::mbd_loglik(
      pars = c(0.2, 0.15, 2.0, 0.1),
      brts = c(1, 2, 2, 2, 3),
      soc = 2, # Crown age
      cond = 1  # Non-extinction
    ),
    "these branching times cannot be generated by a mbd process"
  )
})
