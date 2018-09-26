context("ExpovsLsoda")

test_that("PureBirth theoretical check", {
  skip("PureBirth theoretical check takes too long")
  test_size <- 10;

  soc <- 2; cond <- 1
  lambda <- sample(size = test_size * 2, x = (1:10) * 0.1, replace = TRUE)
  mu <- sample(size = test_size * 2, x = (1:10) * 0.01, replace = TRUE)
  nu <- sample(size = test_size * 2, x = (1:10) * 0.5, replace = TRUE)
  q <- sample(size = test_size * 2, x = (2:14) * 0.01, replace = TRUE)
  test_expo <- rep(0, test_size)
  test_lsoda <- rep(0, test_size)

  for (j in 1:test_size)
  {
    simpars <- c(lambda[j], mu[j], nu[j], q[j]);
    brts <- mbd:::mbd_sim(pars = simpars, soc = soc, age = 10, cond = cond)$brts
    testpars <- c(lambda[test_size + j], mu[test_size + j], nu[test_size + j], q[test_size + j])
    test_expo[j]  <- mbd::mbd_loglik(pars = testpars, brts = brts, soc = soc, cond = cond, methode = "expo")
    test_lsoda[j] <- mbd::mbd_loglik(pars = testpars, brts = brts, soc = soc, cond = cond, methode = "lsoda")
  }
  # print(prod(test_lsoda))
  # print(prod(test_expo) )

  for (i in 1:test_size){
  expect_equal(test_lsoda[i], test_expo[i])
  }
})

# warnings()


