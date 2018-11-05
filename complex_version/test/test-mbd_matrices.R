context("mbd_matrices")

test_that("", {
  
  #test 1
  lambda <- 0.2; mu <- 0.1; nu <- 2; q <- 0.1
  pars <- c(0.2, 0.1, 0, 0.1)
  lx <- 8
  m1 <- mbd::create_a(pars,
                      k = 0,
                      max_number_of_species = lx - 1)
  
  m2 <- matrix(DDD::dd_loglik_M_aux(pars = c(lambda, mu , Inf), lx = lx, k = 0, ddep = 1), lx, lx)
    
  testthat::expect_true(all.equal(m1, m2))
  
  #test 2
  k  <- 6
  lx <- 100
  m1 <- mbd::create_a(pars,
                      k = k,
                      max_number_of_species = lx - 1)
  
  m2 <- matrix(DDD::dd_loglik_M_aux(pars = c(lambda, mu , Inf), lx = lx, k = k, ddep = 1), lx, lx)
  
  testthat::expect_true(all.equal(m1, m2))
  
})

test_that("", {
  
  expect_equal(, )
})

