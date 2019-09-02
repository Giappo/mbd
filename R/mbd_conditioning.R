#' Auxilary function for cond_prob
#' @noRd
cond_prob_nu_matrix <- function(
  q,
  lq
) {
  
  pq <- q #bart's translation
  
  nu_matrix <- matrix(0, nrow = lq, ncol = lq)
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

#' Auxilary function for cond_prob
#' @noRd
cond_prob_rhs1 <- function(
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

#' Auxilary function for cond_prob
#' @noRd
cond_prob_rhs2 <- function(t, x, parms) {
  list(cond_prob_rhs1(
    qvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    k = parms$k
  ))
}

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the conditional probability
#' @author Giovanni Laudanno
#' @export
cond_prob <- function(
  pars,
  brts,
  cond,
  n_0,
  tips_interval = c(n_0 * (cond > 0), Inf),
  lx = 30,
  debug_mode = FALSE
) {
  if (n_0 != 2) {
    stop("This works only for n_0 == 2.")
  }
  if (cond == 0) {
    return(1)
  }
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  tt <- max(abs(brts)) # time between crown age and present
  times <- c(0, tt)
  
  lq <- lx # maximal number of missing species for both m1 and m2
  
  # construct auxiliary matrix
  nu_matrix <- cond_prob_nu_matrix(q = q, lq = lq)
  
  # integrate equations
  parms <- list()
  parms$lambda <- lambda
  parms$mu <- mu
  parms$nu <- nu
  parms$nu_matrix <- nu_matrix
  parms$kk <- 1
  q_0 <- c(y = c(1, rep(0, lq ^ 2 - 1)))
  
  ode_out <- deSolve::ode(
    y = q_0,
    times = times,
    func = cond_prob_rhs2,
    parms = parms,
    method = "lsoda",
    atol = 1e-100,
    rtol = 1e-10,
    tcrit = tt
  )[2, -1]
  q_m1_m2 <- matrix(ode_out, nrow = lq, ncol = lq)
  
  # compute conditioning probability
  m1 <- col(q_m1_m2) - 1
  m2 <- row(q_m1_m2) - 1
  p_m1_m2 <- q_m1_m2 / ((m1 + 1) * (m2 + 1))
  pc <- sum(p_m1_m2)
  
  if (!(pc >= 0 && pc <= 1)) {
    if (debug_mode == TRUE) {
      plot(
        q_t,
        xlab = "m",
        ylab = "Q_m^k(t_p - t_c)",
        main = paste0(
          "pars = ", pars[1], ", ", pars[2], ", ", pars[3], ", ", pars[4]
        ),
        sub = paste0(
          "Conditional probability = ", pc
        )
      )
      cat("The value of pc is: ", pc, "\n")
    } else {
      stop("problems: pc is wrong!")
    }
  } # debug
  
  pc
}
