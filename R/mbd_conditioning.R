#' Auxilary function for cond_prob
#' @author Giovanni Laudanno, Bart Haegeman
#' @noRd
cond_prob_matrices <- function(
  q,
  lq
) {

  pq <- q #bart's translation
  nu_matrix <- matrix(0, nrow = lq, ncol = lq)
  for (m1 in 0:(lq - 1)) {
    for (a1 in 0:floor((m1 + 1) / 2)) { # nolint lintrbot is math's enemy
      aux <- log(m1 + 1) +
        lgamma(m1 - a1 + 1) -
        lgamma(m1 - 2 * a1 + 2) -
        lgamma(a1 + 1)
      aux <- exp(aux)
      aux <- aux * pq ^ a1 * (1 - pq) ^ (m1 + 1 - 2 * a1)
      nu_matrix[m1 + 1, m1 - a1 + 1] <- aux
    }
  }

  empty_qq <- matrix(0, nrow = (lq + 2), ncol = (lq + 2))
  m1 <- col(nu_matrix) - 1
  m2 <- row(nu_matrix) - 1
  list(
    nu_matrix = nu_matrix,
    empty_qq = empty_qq,
    m1 = m1,
    m2 = m2
  )
}

#' Auxilary function for cond_prob
#' @author Giovanni Laudanno, Bart Haegeman
#' @noRd
cond_prob_rhs1 <- function(
  qvec,
  lambda,
  mu,
  nu,
  k,
  nu_matrix,
  m1,
  m2,
  empty_qq
) {
  lq2 <- length(qvec)
  lq <- sqrt(lq2)

  mm <- 2:(lq + 1)
  mm_plus_one <- mm + 1
  mm_minus_one <- mm - 1

  qq <- matrix(qvec, nrow = lq, ncol = lq)
  qq2 <- empty_qq
  qq2[mm, mm] <- qq

  dq1 <- (2 * k + m1 - 1) * qq2[mm, mm_minus_one] +
    (2 * k + m2 - 1) * qq2[mm_minus_one, mm] -
    (2 * k + m1 + m2) * qq # ok

  dq2 <- (m1 + 1) * qq2[mm, mm_plus_one] +
    (m2 + 1) * qq2[mm_plus_one, mm] -
    (2 * k + m1 + m2) * qq # ok

  dq3 <- nu_matrix %*% qq %*% t(nu_matrix) - qq # first cc is m1, t(cc) is m2

  dq <- lambda * dq1 + mu * dq2 + nu * dq3

  dq <- matrix(dq, nrow = lq2, ncol = 1)
  dq
}

#' Auxilary function for cond_prob
#' @author Giovanni Laudanno, Bart Haegeman
#' @noRd
cond_prob_rhs2 <- function(t, x, parms) {
  list(cond_prob_rhs1(
    qvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    k = parms$k,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_qq = parms$empty_qq
  ))
}

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the conditional probability
#' @author Giovanni Laudanno, Bart Haegeman
#' @export
cond_prob <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
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
  matrices <- cond_prob_matrices(q = q, lq = lq)

  # integrate equations
  parms <- list()
  parms$lambda <- lambda
  parms$mu <- mu
  parms$nu <- nu
  parms$kk <- 1
  parms$nu_matrix <- matrices$nu_matrix
  parms$m1 <- matrices$m1
  parms$m2 <- matrices$m2
  parms$empty_qq <- matrices$empty_qq
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
  p_m1_m2 <- q_m1_m2 / ((m1 + 1) * (m2 + 1)) # nolint lintr doesn't know math
  pc <- sum(p_m1_m2)

  if (!(pc >= 0 && pc <= 1)) {
    if (debug_mode == TRUE) {
      plot(
        p_m1_m2,
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
