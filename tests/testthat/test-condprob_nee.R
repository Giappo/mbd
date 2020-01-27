context("condprob_nee")

is_on_ci <- function() {
  is_it_on_appveyor <- Sys.getenv("APPVEYOR") != ""
  is_it_on_travis <- Sys.getenv("TRAVIS") != ""
  is_it_on_appveyor || is_it_on_travis # nolint internal function
}

print_from_global <- function(var = "seed") {
  if (var %in% ls(.GlobalEnv)) {
    cat(var, "is", get(var), "\n")
  }
}

# bd: nu = q = 0 ----
test_that("nu = q = 0", {

  pars <- c(0.2, 0.15, 0, 0)
  brts <- c(5)
  cond <- 1
  n_0 <- 2
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 4)
  for (m in seq_along(mu_vec)) {
    pars[2] <- mu_vec[m]
    test0 <- exp(
      DDD::bd_loglik(
        pars1 = c(pars[1], pars[2], 0, 0),
        pars2 = c(0, 0, 0, 0, n_0),
        brts = brts,
        missnumspec = 0
      ) -
        DDD::bd_loglik(
          pars1 = c(pars[1], pars[2], 0, 0),
          pars2 = c(0, 1, 0, 0, n_0),
          brts = brts,
          missnumspec = 0
        )
    )
    test <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      eq = "nee"
    )
    testthat::expect_equal(
      test0,
      test,
      tolerance = 1e-3
    )
  }

})

# bd: nu = 0 ----
test_that("nu = 0", {

  pars <- c(0.2, 0.15, 0, 0.5)
  brts <- c(7)
  cond <- 1
  n_0 <- 2
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 4)
  for (m in seq_along(mu_vec)) {
    pars[2] <- mu_vec[m]
    test0 <- exp(
      DDD::bd_loglik(
        pars1 = c(pars[1], pars[2], 0, 0),
        pars2 = c(0, 0, 0, 0, n_0),
        brts = brts,
        missnumspec = 0
      ) -
        DDD::bd_loglik(
          pars1 = c(pars[1], pars[2], 0, 0),
          pars2 = c(0, 1, 0, 0, n_0),
          brts = brts,
          missnumspec = 0
        )
    )
    test <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      eq = "nee"
    )
    testthat::expect_equal(
      test0,
      test,
      tolerance = 1e-3
    )
  }

})

# bd: q = 0 ----
test_that("q = 0", {

  pars <- c(0.2, 0.15, 3, 0)
  brts <- c(1)
  cond <- 1
  n_0 <- 2
  mu_vec <- seq(from = 0.05, to = pars[1], length.out = 4)
  for (m in seq_along(mu_vec)) {
    pars[2] <- mu_vec[m]
    test0 <- exp(
      DDD::bd_loglik(
        pars1 = c(pars[1], pars[2], 0, 0),
        pars2 = c(0, 0, 0, 0, n_0),
        brts = brts,
        missnumspec = 0
      ) -
        DDD::bd_loglik(
          pars1 = c(pars[1], pars[2], 0, 0),
          pars2 = c(0, 1, 0, 0, n_0),
          brts = brts,
          missnumspec = 0
        )
    )
    test <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      eq = "nee"
    )
    testthat::expect_equal(
      test0,
      test,
      tolerance = 1e-3
    )
  }

})

# nee vs sims ----
test_that("nee vs sims", {

  skip("Need to simulate first on local")

  max_seed <- 50
  min_n_sims <- 2e3
  max_n_sims <- 2e4
  for (seed in 1:max_seed) {
    set.seed(seed)
    print_from_global("seed")
    lambda <- runif(n = 1, min = 0.05, max = 0.3)
    mu <- runif(n = 1, min = 0, max = lambda)
    nu <- runif(n = 1, min = 0.3, max = 2.3)
    q <- runif(n = 1, min = 0.05, max = 0.35)
    pars <- c(lambda, mu, nu, q)
    brts <- c(8)

    pc <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      eq = "nee"
    )
    n_sims <- ceiling(
      max(
        min(1 / abs(pc - 0.5) ^ 2, max_n_sims),
        min_n_sims
      )
    )
    print_from_global("n_sims")
    pc_sim <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = n_sims,
      eq = "sim"
    )
    print_from_global("pc_sim")
    print_from_global("pc")
    testthat::expect_equal(pc, pc_sim, tolerance = 1e-2)
  }

  for (seed in 1:max_seed) {
    if (seed == 15 || seed == 27 || seed == 35) next # too slow for simulations
    set.seed(seed)
    print_from_global("seed")
    lambda <- runif(n = 1, min = 0.05, max = 1.5)
    mu <- runif(n = 1, min = 0, max = lambda)
    nu <- runif(n = 1, min = 0.3, max = 2.9)
    q <- runif(n = 1, min = 0.05, max = 0.4)
    pars <- c(lambda, mu, nu, q)
    brts <- c(5)

    pc <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      eq = "nee"
    )
    n_sims <- ceiling(
      max(
        min(1 / abs(pc - 0.5) ^ 2, max_n_sims),
        min_n_sims
      )
    )
    print_from_global("n_sims")
    pc_sim <- mbd::calculate_condprob(
      pars = pars,
      brts = brts,
      lx = n_sims,
      eq = "sim"
    )
    print_from_global("pc_sim")
    print_from_global("pc")
    testthat::expect_equal(pc, pc_sim, tolerance = 1e-2)
  }

})
