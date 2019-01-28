#' @title Transition matrix for cond = 2
#' @description Creates the A matrix,
#' used for the cond = 2 calculation
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @author Giovanni Laudanno
#' @export
create_a_cond_2 <- function(pars, k, lx) {
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]

  testit::assert(lx < 2 ^ 31)
  nvec <- 0:lx
  m_visible <- 2
  mvec <- 0:(m_visible - 1)
  nu_terms <- rep(1, lx + 1)
  for (m in mvec) {
    nu_terms <- nu_terms -
      choose(k + nvec, m) * q ^ m * (1 - q) ^ (nvec + k - m)
  }

  m <- matrix(0, lx + 1, lx + 1)
  diag(m) <- (-nu) * nu_terms - (lambda + mu) * (k + nvec)
  m[row(m) == col(m) - 1] <- mu * nvec[2:(lx + 1)]
  m[row(m) == col(m) + 1] <- m[row(m) == col(m) + 1] +
    lambda * (nvec[1:(lx)] + 2 * k)
  m[length(nvec), length(nvec)] <- (-mu) * (k + nvec[length(nvec)])

  m
}

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the pc. If \code{is.nan(log(pc))} the log-likelihood estimation
#'   by \link{mbd_loglik} is -Inf
#' @author Giovanni Laudanno
#' @noRd
calculate_conditional_prob <- function(
  brts,
  pars,
  cond,
  lx = 1000,
  n_0 = 2,
  tips_interval = c(0, Inf),
  methode = "lsodes",
  abstol = 1e-16,
  reltol = 1e-10
) {
  if (!(cond %in% mbd_conds())) {
    stop("This conditioning is not implemented.")
  }
  pc0 <- pc1 <- pc2 <- 1
  if (cond == 0) {
    return(pc0)
  }
  total_time <- max(abs(brts))
  m <- 0:lx
  one_over_cm <- (3 * (m + 1)) / (m + 3)
  one_over_qm_binom <- 1 / choose(m + n_0, n_0)
  q_i <- c(1, rep(0, lx))
  testit::assert(length(one_over_cm) == length(m))
  testit::assert(length(one_over_qm_binom) == length(m))
  testit::assert(length(q_i) == length(m))
  # creating a_matrix
  matrix_a <- create_a(pars = pars, k = n_0, lx = lx) # nolint internal function
  # integrating the starting q_vector to t_p
  a2_v1 <- a_operator(
    q_vector = q_i,
    transition_matrix = matrix_a,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    abstol = abstol,
    reltol = reltol
  )

  total_product <- a2_v1 * one_over_cm * one_over_qm_binom
  missingspecies_min <- max(tips_interval[1] - n_0, 0)
  missingspecies_max <- min(tips_interval[2] - n_0, lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  pc1 <- sum(total_product[tips_components[1]:tips_components[2]])

  if (cond == 2) {
    # integrating the starting q_vector to t_p
    a2_v1_2 <- a_operator(
      q_vector = q_i,
      transition_matrix = create_a_cond_2(
        pars = pars,
        k = n_0,
        lx = lx
      ),
      time_interval = total_time,
      precision = 250L,
      methode = methode,
      abstol = abstol,
      reltol = reltol
    )
    total_product_2 <- a2_v1_2 * one_over_cm * one_over_qm_binom
    missingspecies_min <- max(tips_interval[1] - n_0, 0)
    missingspecies_max <- min(tips_interval[2] - n_0, lx)
    # +1 is because of the zero-th component
    tips_components <- 1 + c(missingspecies_min, missingspecies_max)
    pc2 <- pc1 - sum(total_product_2[tips_components[1]:tips_components[2]])
  }

  pc1 * (cond == 1) + pc2 * (cond == 2)
}
