context("TestSims")

test_that("Run sim for lambda = 0.0", {

  test_size <- 1
  test_sim <- vector("list")
  soc <- 2
  nu <- 0.2
  mu <- 0.1
  q <- 0.3

  test_sim <- mbd::mbd_sim0(
    pars = c(nu, mu, q), soc = soc, age = 10, cond = 1
  )$l_matrix

  #check if all the species are born before their parents die
  back_to_the_future_test  <-  1
  l_matrix <- test_sim
  parents <- abs(l_matrix[, 2])
  back_to_the_future_test <- back_to_the_future_test *
    prod(l_matrix[-c(1:soc), 1] >= l_matrix[parents[-c(1:soc)], 4])

  #check if the conditioning on survival works
  conditioning_test <- 1
  l_matrix <- test_sim
  check_signs <- sort(unique(sign(l_matrix[, 3])))
  conditioning_test <- conditioning_test * prod(check_signs == c(-1, 1))

  expect_equal(back_to_the_future_test, 1)
  expect_equal(conditioning_test, 1)
})
