#' Creates log nu matrix for the Q-equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_qeq_log_nu_mat <- function(lx) {
  nvec <- 1:lx - 1
  log_nu_mat <- matrix(-Inf, nrow = lx, ncol = lx)
  for (m1 in nvec) {
    for (a1 in 0:floor((m1 + 1) / 2)) { # nolint lintrbot is math's enemy
      aux <- log(m1 + 1) +
        lfactorial(m1 - a1) -
        lfactorial(m1 - 2 * a1 + 1) -
        lfactorial(a1)
      log_nu_mat[m1 + 1, m1 - a1 + 1] <- aux
    }
  }
  rownames(log_nu_mat) <- paste0("m1=", nvec)
  colnames(log_nu_mat) <- paste0("n1=", nvec)
  log_nu_mat
}

#' Creates log q-matrix for the Q-equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_qeq_log_q_mat <- function(lx, q) {
  nvec <- 1:lx - 1
  log_q_mat <- matrix(-Inf, nrow = lx, ncol = lx)
  for (m1 in nvec) {
    for (a1 in 0:floor((m1 + 1) / 2)) { # nolint lintrbot is math's enemy
      aux <- a1 * log(q) + (m1 + 1 - 2 * a1) * log(1 - q)
      log_q_mat[m1 + 1, m1 - a1 + 1] <- aux
    }
  }
  rownames(log_q_mat) <- paste0("m1=", nvec)
  colnames(log_q_mat) <- paste0("n1=", nvec)
  log_q_mat
}

#' Creates log nu matrix for the P-equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_peq_log_nu_mat <- function(lx) {
  nvec <- 1:lx - 1
  log_nu_mat <- matrix(-Inf, nrow = lx, ncol = lx)
  for (m1 in nvec) {
    for (n1 in 0:m1) {
      aux <- lchoose(n1, max(0, m1 - n1))
      log_nu_mat[m1 + 1, n1 + 1] <- aux
    }
  }
  rownames(log_nu_mat) <- paste0("m1=", nvec)
  colnames(log_nu_mat) <- paste0("n1=", nvec)
  log_nu_mat
}

#' Creates log q matrix for the P-equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_peq_log_q_mat <- function(lx, q) {
  nvec <- 1:lx - 1
  log_q_mat <- matrix(-Inf, nrow = lx, ncol = lx)
  for (m1 in nvec) {
    for (n1 in 0:m1) {
      aux <- (m1 - n1) * log(q) + (2 * n1 - m1) * log(1 - q)
      log_q_mat[m1 + 1, n1 + 1] <- aux
    }
  }
  rownames(log_q_mat) <- paste0("m1=", nvec)
  colnames(log_q_mat) <- paste0("n1=", nvec)
  log_q_mat
}

#' Creates log nu matrix for the specified equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_log_nu_mat <- function(lx, eq) {
  if (!(eq == "p_eq" || eq == "q_eq")) {
    stop("It is either Q- or P-equation!")
  }
  if (eq == "p_eq") {
    log_nu_mat <- condprob_peq_log_nu_mat(lx)
  }
  if (eq == "q_eq") {
    log_nu_mat <- condprob_qeq_log_nu_mat(lx)
  }
  log_nu_mat
}

#' Creates log q matrix for the specified equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_log_q_mat <- function(lx, q, eq) {
  if (!(eq == "p_eq" || eq == "q_eq")) {
    stop("It is either Q- or P-equation!")
  }
  if (eq == "p_eq") {
    log_q_mat <- condprob_peq_log_q_mat(lx = lx, q = q)
  }
  if (eq == "q_eq") {
    log_q_mat <- condprob_qeq_log_q_mat(lx = lx, q = q)
  }
  log_q_mat
}

#' Creates vector of parameters for the FORTRAN condprob code
#' @return a vector composed of (in order): lambda, mu, nu, components of the
#'  nu-q matrix, components of the matrix of columns (m1_mat), components of the
#'  matrix of rows (m2_mat), components of the empty framed matrix (empty_mat).
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_parmsvec <- function(
  pars,
  log_nu_mat,
  log_q_mat,
  lx,
  eq
) {
  lx2 <- lx ^ 2
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]

  log_nu_mat <- log_nu_mat[1:lx, 1:lx]
  log_q_mat <- log_q_mat[1:lx, 1:lx]
  nu_q_mat <- exp(log_nu_mat + log_q_mat)

  parmsvec <- c(
    lambda,
    mu,
    nu,
    matrix(nu_q_mat, nrow = lx2, ncol = 1)
  )
  parmsvec
}

# Differentials ----

#' dp total
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_dp <- function(
  pvec,
  parmsvec
) {
  lx2 <- length(pvec)
  lx <- sqrt(lx2)
  pp <- matrix(pvec, lx, lx)
  mm <- 2:(lx + 1)

  lambda <- parmsvec[1]
  mu <- parmsvec[2]
  nu <- parmsvec[3]
  nu_q_mat <- parmsvec[(3 + 1):(3 + lx2)]
  dim(nu_q_mat) <- c(lx, lx)

  pp2 <- matrix(0, lx + 2, lx + 2)
  pp2[mm, mm] <- pp
  mvec <- 1:lx - 1

  dp_lambda <- matrix(0, lx, lx)
  for (mm2 in mvec) {
    for (mm1 in mvec) {
      i <- mm2 + 1
      j <- mm1 + 1
      dp_lambda[i, j] <-
        (mm1 - 1) * pp2[i + 1, j] +
        (mm2 - 1) * pp2[i, j + 1] -
        (mm1 + mm2) * pp2[i + 1, j + 1]
    }
  }

  dp_mu <- matrix(0, lx, lx)
  for (mm2 in mvec) {
    for (mm1 in mvec) {
      i <- mm2 + 1
      j <- mm1 + 1
      dp_mu[i, j] <-
        (mm1 + 1) * pp2[i + 1, j + 2] +
        (mm2 + 1) * pp2[i + 2, j + 1] -
        (mm1 + mm2) * pp2[i + 1, j + 1]
    }
  }

  nu_q_mat2 <- t(nu_q_mat)
  dp_nu <- aux1 <- aux2 <-  matrix(0, lx, lx)
  for (m1 in 1:lx) {
    for (n2 in 1:lx) {
      sum1 <- 0
      for (n1 in 1:lx) {
        sum1 <- sum1 + nu_q_mat[m1, n1] * pp[n1, n2]
      }
      aux1[m1, n2] <- sum1
    }
  }
  for (m1 in 1:lx) {
    for (m2 in 1:lx) {
      sum1 <- 0
      for (n2 in 1:lx) {
        sum1 <- sum1 + aux1[m1, n2] * nu_q_mat2[n2, m2]
      }
      aux2[m1, m2] <- sum1
    }
  }
  dp_nu <- aux2 - pp

  dp <- lambda * dp_lambda + mu * dp_mu + nu * dp_nu
  dim(dp) <- c(lx2, 1)
  return(dp)
}

#' dq total
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_dq <- function(
  qvec,
  parmsvec
) {
  lx2 <- length(qvec)
  lx <- sqrt(lx2)
  qq <- matrix(qvec, lx, lx)
  mm <- 2:(lx + 1)

  lambda <- parmsvec[1]
  mu <- parmsvec[2]
  nu <- parmsvec[3]
  nu_q_mat <- parmsvec[(3 + 1):(3 + lx2)]
  dim(nu_q_mat) <- c(lx, lx)

  qq2 <- matrix(0, lx + 2, lx + 2)
  qq2[mm, mm] <- qq
  mvec <- 1:lx - 1

  dq_lambda <- matrix(0, lx, lx)
  for (mm2 in mvec) {
    for (mm1 in mvec) {
      i <- mm2 + 1
      j <- mm1 + 1
      dq_lambda[i, j] <-
        (2 + mm1 - 1) * qq2[i + 1, j] +
        (2 + mm2 - 1) * qq2[i, j + 1] -
        (2 + mm1 + mm2) * qq2[i + 1, j + 1]
    }
  }

  dq_mu <- matrix(0, lx, lx)
  for (mm2 in mvec) {
    for (mm1 in mvec) {
      i <- mm2 + 1
      j <- mm1 + 1
      dq_mu[i, j] <-
        (mm1 + 1) * qq2[i + 1, j + 2] +
        (mm2 + 1) * qq2[i + 2, j + 1] -
        (2 + mm1 + mm2) * qq2[i + 1, j + 1]
    }
  }

  nu_q_mat2 <- t(nu_q_mat)
  dq_nu <- aux1 <- aux2 <-  matrix(0, lx, lx)
  for (m1 in 1:lx) {
    for (n2 in 1:lx) {
      sum1 <- 0
      for (n1 in 1:lx) {
        sum1 <- sum1 + nu_q_mat[m1, n1] * qq[n1, n2]
      }
      aux1[m1, n2] <- sum1
    }
  }
  for (m1 in 1:lx) {
    for (m2 in 1:lx) {
      sum1 <- 0
      for (n2 in 1:lx) {
        sum1 <- sum1 + aux1[m1, n2] * nu_q_mat2[n2, m2]
      }
      aux2[m1, m2] <- sum1
    }
  }
  dq_nu <- aux2 - qq

  dq <- lambda * dq_lambda + mu * dq_mu + nu * dq_nu
  dim(dq) <- c(lx2, 1)
  return(dq)
}
