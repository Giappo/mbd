probcond_nu_matrix <- function(
  q,
  lx
) {
  lq <- lx #bart's translation
  pq <- q #bart's translation
  
  cc <- matrix(0, nrow = lq, ncol = lq)
  for (m1 in 0:(lq - 1)) {
    for (a1 in 0:floor((m1 + 1) / 2)) {
      aux <- log(m1 + 1) +
        lgamma(m1 - a1 + 1) -
        lgamma(m1 - 2 * a1 + 2) -
        lgamma(a1 + 1)
      aux <- exp(aux)
      aux <- aux * pq ^ a1 * (1 - pq) ^ (m1 + 1 - 2 * a1)
      nu_matrix[m1 + 1, m1 - a1 + 1] <- aux
    }
  }
  nu_matrix
}

probcond_rhs1 <- function(
  qvec,
  lambda,
  mu,
  nu,
  nu_matrix,
  k
) {
  lq2 <- length(qvec)
  lq <- sqrt(lq2)
  qq <- matrix(qvec, nrow = lq, ncol = lq)
  m1 <- col(qq) - 1
  m2 <- row(qq) - 1
  
  dq1 <- (2 * k + m1 - 1) * cbind(rep(0, lq), qq[, 1:(lq - 1)]) +
    (2 * k + m2 - 1) * rbind(rep(0, lq), qq[1:(lq - 1), ]) -
    (2 * k + m1 + m2) * qq # ok
  
  dq2 <- (2 * k + m1 - 1) * cbind(qq[, 2:lq], rep(0, lq)) +
    (2 * k + m2 - 1) * rbind(qq[2:lq, ], rep(0, lq)) -
    (2 * k + m1 + m2) * qq # ok
  
  dq3 <- nu_matrix %*% qq %*% t(nu_matrix) - qq # first cc is m1, t(cc) is m2
  
  dq <- lambda * dq1 + mu * dq2 + nu * dq3
  
  dq <- matrix(dq, nrow = lq2, ncol = 1)
  dq
}
probcond_rhs2 <- function(t, x, parms) {
  list(probcond_rhs1(
    qvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    k = parms$k
  ))
}
probcond <- function(
  pars,
  brts
) {

  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  times <- c(0, max(abs(brts))
  la=0.0204942104;
  mu=0.0001333249;
  nu=1.5728643216;
  pq=0.0787076385;
  tt=10; # time between crown age and present
  lq=4; # maximal number of missing species
  
  # construct auxiliary matrix
  cc <- probcond2_matr(pq = pq, lq = lq); cc
  
  # integrate equations
  parms <- list()
  parms$la <- la
  parms$mu <- mu
  parms$nu <- nu
  parms$cc <- cc
  parms$kk <- 1
  x0 <- c(y = c(1, rep(0, lq ^ 2 - 1)))
  # x0 <- matrix(1:(lq ^ 2), nrow = lq^2, ncol = 1)
  times <- seq(0, tt, length.out = 2)
  X4 <- deSolve::ode(
    y = x0,
    times = times,
    func = probcond_rhs,
    parms = parms,
    method = "ode45",
    atol = 1e-100,
    rtol=1e-10,
    tcrit = tt
  )[2, -1]; X4
  Q4 <- matrix(X4, nrow = lq, ncol = lq)
  % compute conditioning probability
  m1 <- col(Q4) - 1
  m2 <- row(Q4) - 1
  Pc4 <- Q4 / ((m1 + 1) * (m2 + 1))
  sum(Pc4)
}
  