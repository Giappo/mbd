## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
while (!require(devtools)) {install.packages("devtools")};
devtools::install_github("Giappo/mbd@debug");
library(mbd);

## ------------------------------------------------------------------------
methode <- "ode45" # ode method
loglik <- mbd_loglik(
  pars = pars,
  brts = brts,
  n_0 = n_0,
  cond = cond,
  tips_interval = tips_interval,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
lambda <- 0.2; mu <- 0.15; nu <- 1; q <- 0.1; # mbd parameters
pars <- c(lambda, mu, nu, q)
lx <- 7 # matrix dimension
k <- 2 # k-species
b <- 1 # births
n_matrix <- mbd:::create_n(pars = pars, k = k, b = b, lx = lx)
n_matrix

## ------------------------------------------------------------------------
lambda <- 0.2; mu <- 0.15; nu <- 1; q <- 0.1;
pars <- c(lambda, mu, nu, q)
lx <- 9
k <- 2
transition_matrix <- mbd:::create_a(pars = pars, lx = lx, k = k)
transition_matrix

## ------------------------------------------------------------------------
# m < n
for (n in 1:lx) {
  for (m in 0:(n - 1)) {
    testthat::expect_equal(
      (m == n - 1) * mu * (m + 1),
      transition_matrix[m + 1, n + 1]
    )
  }
}
# m = n
m <- 0:lx
testthat::expect_equal(
  diag(transition_matrix),
  -(lambda + mu) * (m + k) - nu * (1 - (1 - q) ^ (m + k))
)
# m > n
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
# m < n
testthat::expect_true(
  all(b_matrix[upper.tri(b_matrix, diag = FALSE)] == 0)
)
# m >= n
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

## ------------------------------------------------------------------------
brts <- c(10, 9, 7, 6, 5) # branching times
pars <- c(2, 0.05, 1.00, 0.1) # mbd parameters = c(lambda, mu, nu, q)
n_0 <- 2 # starting species
cond <- 1 # conditioning? 1 is yes, 0 is no.
tips_interval <- c(n_0 * (cond > 0), Inf) # how many tips are allowed in tree?
lx <- 200 # dimension of the transition matrix
abstol <- 1e-16 # absolute tolerance for ode
reltol <- 1e-10 # relative tolerance for ode

## ------------------------------------------------------------------------
methode <- "lsodes" # ode method
loglik <- mbd_loglik(
  pars = pars,
  brts = brts,
  n_0 = n_0,
  cond = cond,
  tips_interval = tips_interval,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
methode <- "lsoda" # ode method
loglik <- mbd_loglik(
  pars = pars,
  brts = brts,
  n_0 = n_0,
  cond = cond,
  tips_interval = tips_interval,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
brts <- c(10, 9, 7, 6, 5)
pars <- c(0.0204942104, 0.0001333249, 1.5728643216, 0.0787076385)
n_0 <- 2
cond <- 1
tips_interval <- c(n_0 * (cond > 0), Inf)
lx <- 200
methode <- "lsodes"
abstol <- 1e-16
reltol <- 1e-10

## ------------------------------------------------------------------------
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
methode <- "ode45" # ode method
loglik <- mbd_loglik(
  pars = pars,
  brts = brts,
  n_0 = n_0,
  cond = cond,
  tips_interval = tips_interval,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
methode <- "lsodes" # ode method
loglik <- mbd_loglik(
  pars = pars,
  brts = brts,
  n_0 = n_0,
  cond = cond,
  tips_interval = tips_interval,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
methode <- "lsoda" # ode method
loglik <- mbd_loglik(
  pars = pars,
  brts = brts,
  n_0 = n_0,
  cond = cond,
  tips_interval = tips_interval,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
brts <- c(10, 9, 7, 6, 5)
pars <- c(0.0204942104, 0.0001333249, 1.5728643216, 0.0787076385)
n_0 <- 2
cond <- 1
tips_interval <- c(n_0 * (cond > 0), Inf)
lx <- 200
methode <- "lsodes"
abstol <- 1e-16
reltol <- 1e-10

## ------------------------------------------------------------------------
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = 1000,
  tips_interval = tips_interval,
  methode = methode,
  abstol = abstol * 10 ^ -3,
  reltol = reltol * 10 ^ -3,
  debug_mode = TRUE
)
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = "ode45",
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = "lsoda",
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
out <- mbd:::mbd_calculate_q_vector(
  pars = pars,
  brts = brts[1],
  n_0 = n_0,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol
)
sum(out$q_f)
out$q_t

## ------------------------------------------------------------------------
methode <- "lsoda" # ode method
loglik <- mbd_loglik(
  pars = pars,
  brts = brts,
  n_0 = n_0,
  cond = cond,
  tips_interval = tips_interval,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
brts <- c(10, 9, 7, 6, 5)
pars <- c(0.0204942104, 0.0001333249, 1.5728643216, 0.0787076385)
n_0 <- 2
cond <- 1
tips_interval <- c(n_0 * (cond > 0), Inf)
lx <- 200
methode <- "lsodes"
abstol <- 1e-16
reltol <- 1e-10

## ------------------------------------------------------------------------
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = methode,
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = 1000,
  tips_interval = tips_interval,
  methode = methode,
  abstol = abstol * 10 ^ -3,
  reltol = reltol * 10 ^ -3,
  debug_mode = TRUE
)
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = "ode45",
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = "lsoda",
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
out <- mbd:::mbd_calculate_q_vector(
  pars = pars,
  brts = brts[1],
  n_0 = n_0,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol
)
sum(out$q_f)
out$q_t

## ------------------------------------------------------------------------
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = 1000,
  tips_interval = tips_interval,
  methode = methode,
  abstol = abstol * 10 ^ -3,
  reltol = reltol * 10 ^ -3,
  debug_mode = TRUE
)
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = "ode45",
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)
pc <- mbd:::calculate_conditional_prob(
  brts = brts,
  pars = pars,
  cond = cond,
  n_0 = n_0,
  lx = lx,
  tips_interval = tips_interval,
  methode = "lsoda",
  abstol = abstol,
  reltol = reltol,
  debug_mode = TRUE
)

## ------------------------------------------------------------------------
out <- mbd:::mbd_calculate_q_vector(
  pars = pars,
  brts = brts[1],
  n_0 = n_0,
  lx = lx,
  methode = methode,
  abstol = abstol,
  reltol = reltol
)
sum(out$q_f)
out$q_t

