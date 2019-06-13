context("mbd_matrices")

# a_matrix ----
test_that("a_matrix", {

  # Test 1
  lambda <- 0.2; mu <- 0.1; nu <- 0; q <- 0.1
  pars <- c(lambda, mu, nu, q)
  lx <- 8
  m1 <- mbd::create_a(
    pars,
    k = 0,
    lx = lx
  )
  m2 <- matrix(
    DDD:::dd_loglik_M_aux(
      pars = c(lambda, mu, Inf),
      lx = lx + 1,
      k = 0,
      ddep = 1),
    nrow = lx + 1,
    ncol = lx + 1
  )
  testthat::expect_true(all.equal(m1, m2))
  testthat::expect_true(all(dim(m1) == lx + 1))
  testthat::expect_true(all(diag(m1) <= 0))
  testthat::expect_true(
    all(m1[row(m1) < col(m1) - 1] == 0)
  )

  # Test 2
  k  <- 6
  lx <- 100
  m1 <- mbd::create_a(
    pars,
    k = 0,
    lx = lx
  )
  m2 <- matrix(
    DDD:::dd_loglik_M_aux(
      pars = c(lambda, mu, Inf),
      lx = lx + 1,
      k = 0,
      ddep = 1
    ),
    nrow = lx + 1,
    ncol = lx + 1
  )
  testthat::expect_true(all.equal(m1, m2))
  testthat::expect_true(all(dim(m1) == lx + 1))
  testthat::expect_true(all(diag(m1) <= 0))
  testthat::expect_true(
    all(m1[row(m1) < col(m1) - 1] == 0)
  )
})

# b_matrix ----
test_that("b_matrix", {

  # Test 1
  lambda <- 0.2; mu <- 0.1; nu <- 2; q <- 0.1
  pars <- c(lambda, mu, nu, q)
  k <- 0
  b <- 0
  lx <- 20
  m1 <- mbd::create_b(
    pars,
    k = k,
    b = b,
    lx = lx
  ); m1

  testthat::expect_true(
    all(m1[col(m1) > row(m1)] == 0)
  )
  testthat::expect_true(all(dim(m1) == lx + 1))
  testthat::expect_true(all(m1 >= 0))
  if (k == 0) {
    testthat::expect_true(
      all(
        m1[row(m1) > 2 * col(m1)] == 0
      )
    )
  }

  # Test 2
  lambda <- 0.2; mu <- 0.1; nu <- 2; q <- 0.1
  pars <- c(lambda, mu, nu, q)
  k <- 0
  b <- 1
  lx <- 8
  testthat::expect_error(
    m1 <- mbd::create_b(
      pars,
      k = k,
      b = b,
      lx = lx
    ),
    "you can't have more births than species present in the phylogeny"
  )

  # Test 3
  lambda <- 0.2; mu <- 0.1; nu <- 2; q <- 0.1
  pars <- c(lambda, mu, nu, q)
  k <- 3
  b <- 1
  lx <- 30
  m1 <- mbd::create_b(
    pars,
    k = k,
    b = b,
    lx = lx
  ); m1
  testthat::expect_true(
    # b is lower triangular
    all(m1[col(m1) > row(m1)] == 0)
  )
  testthat::expect_true(all(dim(m1) == lx + 1))
  testthat::expect_true(all(m1 >= 0))
  testthat::expect_true(
    # you cannot have more speciations than k - b
    all(
      m1[row(m1) > 2 * col(m1) + (k - b)] == 0
    )
  )
})
