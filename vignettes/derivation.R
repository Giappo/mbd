## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
lambda <- 0.2; mu <- 0.15; nu <- 1; q <- 0.1;
pars <- c(lambda, mu, nu, q)
lx <- 9
k <- 2
transition_matrix <- mbd:::create_a(pars = pars, lx = lx, k = k)
transition_matrix

## ------------------------------------------------------------------------
for (m in 1:lx) {
  for (n in 0:(m - 1)) {
    a <- m - n
    j <- 0:min(m - n, k)
    entry_m_n <- 
      (m - n == 1) * lambda * (m + 2 * k - 1) +
      nu *
      (1 - q) ^ (k + m) *
      q ^ a * 
      (1 - q) ^ (-2 * a) *
      sum(
        2 ^ j * choose(k, j) * choose(m - a, a - j)
      )
    testthat::expect_equal(
      entry_m_n,
      transition_matrix[m + 1, n + 1]
    )
  }
}

## ------------------------------------------------------------------------
t_1 <- 0
t_2 <- 2
q_1 <- c(1, rep(0, lx))
q_2 <- mbd:::a_operator(
  q_vector = q_1,
  transition_matrix = transition_matrix,
  time_interval = t_2 - t_1
)
print(paste(c("Vector Q at time t_1 is:", q_1), collapse = " "))
print(paste(c("Vector Q at time t_2 is:", signif(q_2, digits = 2)), collapse = " "))
plot(q_2, ylab = "Q(t_2)", xlab = "m")

## ------------------------------------------------------------------------
b <- 1
b_matrix <- mbd:::create_b(pars = pars, lx = lx, k = k, b = b)
b_matrix

## ------------------------------------------------------------------------
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

