#' Creates m1 matrix for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_m1_mat <- function(lx_top) {
  nu_mat <- matrix(0, nrow = lx_top, ncol = lx_top)
  m1_mat <- col(nu_mat) - 1
  m1_mat
}

#' Creates m2 matrix for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_m2_mat <- function(m1_mat) {
  m2_mat <- t(m1_mat)
  m2_mat
}

#' Creates empty matrix for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_empty_mat <- function(lx_top) {
  empty_mat <- matrix(0, nrow = (lx_top + 2), ncol = (lx_top + 2))
  empty_mat
}

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
  m1_mat,
  m2_mat,
  empty_mat,
  log_nu_mat,
  log_q_mat,
  pars,
  lx,
  eq
) {
  lx2 <- lx ^ 2
  lxframe <- lx + 2
  lxframe2 <- lxframe ^ 2
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]

  log_nu_mat <- log_nu_mat[1:lx, 1:lx]
  log_q_mat <- log_q_mat[1:lx, 1:lx]
  nu_q_mat <- exp(log_nu_mat + log_q_mat)
  m1_mat <- m1_mat[1:lx, 1:lx]
  m2_mat <- m2_mat[1:lx, 1:lx]
  empty_mat <- empty_mat[1:lxframe, 1:lxframe]

  parmsvec <- c(
    lambda,
    mu,
    nu,
    matrix(nu_q_mat, nrow = lx2, ncol = 1),
    matrix(m1_mat, nrow = lx2, ncol = 1),
    matrix(m2_mat, nrow = lx2, ncol = 1),
    matrix(empty_mat, nrow = lxframe2, ncol = 1)
  )
  parmsvec
}
