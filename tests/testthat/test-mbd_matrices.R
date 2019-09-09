context("mbd_matrices")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

# a_matrix ----
test_that("a_matrix: mbd vs ddd", {

  # mbd and ddd yield same matrices for nu = q = 0
  # here k = 0
  lambda <- 0.2; mu <- 0.1; nu <- 0; q <- 0.1
  pars <- c(lambda, mu, nu, q)
  lx <- 8
  m1 <- mbd::create_a(
    pars,
    k = 0,
    lx = lx,
    no_species_out_of_the_matrix = TRUE
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

  # mbd and ddd yield same matrices for nu = q = 0
  # here k != 0
  k  <- 6
  lx <- 100
  m1 <- mbd::create_a(
    pars,
    k = 0,
    lx = lx,
    no_species_out_of_the_matrix = TRUE
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

test_that("a_matrix: # entries are in accordance with the theory", {
  lambda <- 3; mu <- 1.5; nu <- 7; q <- 0.4
  pars <- c(lambda, mu, nu, q)
  k <- 2
  lx <- 50 + 150 * is_on_ci()
  no_species_out_of_the_matrix <- TRUE
  a_matrix <- mbd::create_a(
    pars,
    k = k,
    lx = lx,
    no_species_out_of_the_matrix = no_species_out_of_the_matrix
  )
  # lower triangular matrix: m > n
  for (m in 1:lx) {
    for (n in 0:(m - 1)) {
      j <- 0:min(m - n, k)
      entry_m_n <-
        (m == n + 1) * lambda * (m - 1 + 2 * k) +
        nu *
        (1 - q) ^ k *
        q ^ (m - n) *
        (1 - q) ^ (2 * n - m) *
        sum(
          (2 ^ j) * choose(k, j) * choose(n, m - n - j)
        )
      testthat::expect_equal(
        entry_m_n,
        a_matrix[m + 1, n + 1]
      )
    }
  }
  # main diagonal
  m <- 0:lx
  nu_terms <- rep(0, lx + 1)
  if (no_species_out_of_the_matrix == TRUE) {
    for (n in 0:(lx - 1)) {
      limit <- min(n + k, lx - n)
      avec <- 1:limit
      nu_terms[n + 1] <-
        sum(
          choose(n + k, avec) * (q ^ avec) * (1 - q) ^ (n + k - avec)
        )
    }
  } else {
    nu_terms <- (1 - (1 - q) ^ (m + k))
  }
  entry_m_m <-
    -nu * nu_terms +
    -mu * (m + k) +
    -lambda * (m + k) +
    no_species_out_of_the_matrix * lambda * c(rep(0, lx), 1) * (m + k)

  testthat::expect_equal(
    diag(a_matrix),
    entry_m_m
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
  )

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
  )
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

test_that("b_matrix: # entries are in accordance with the theory", {
  lambda <- 4; mu <- 2; nu <- 7.2; q <- 0.5
  pars <- c(lambda, mu, nu, q)
  k <- 8; b <- 3
  lx <- 50 + 150 * is_on_ci()
  b_matrix <- mbd::create_b(
    pars,
    k = k,
    lx = lx,
    b = b
  )
  # lower triangular matrix with diagonal: m >= n
  for (m in 0:lx) {
    for (n in 0:m) {
      a <- m - n
      j <- 0:min(m - n, k)
      entry_m_n <-
        (m == n) * (b == 1) * lambda * k +
        nu *
        choose(k, b) * q ^ b *
        (1 - q) ^ (k - b + m) *
        q ^ a *
        (1 - q) ^ (-2 * a) *
        sum(
          2 ^ j * choose(k - b, j) * choose(m - a, a - j)
        )
      testthat::expect_equal(
        entry_m_n,
        b_matrix[m + 1, n + 1]
      )
    }
  }
})
