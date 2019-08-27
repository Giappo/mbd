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
#' @return the conditional probability
#' @author Giovanni Laudanno
#' @export
calculate_conditional_prob <- function(
  brts,
  pars,
  cond,
  lx = 1000,
  n_0 = 2,
  tips_interval = c(n_0 * (cond > 0), Inf),
  methode = "lsodes",
  abstol = 1e-16,
  reltol = 1e-10,
  debug_mode = FALSE
) {
  check_cond(cond = cond, tips_interval = tips_interval, n_0 = n_0)
  if (cond == 0) {
    return(1)
  }
  total_time <- max(abs(brts))
  m <- 0:lx
  one_over_cm <- (3 * (m + 1)) / (m + 3)
  one_over_qm_binom <- 1 / choose(m + n_0, n_0)
  q_i <- c(1, rep(0, lx)); names(q_i) <- paste0("Q", 0:lx)
  testit::assert(length(one_over_cm) == length(m))
  testit::assert(length(one_over_qm_binom) == length(m))
  testit::assert(length(q_i) == length(m))
  # creating a_matrix
  matrix_a <- create_a(pars = pars, k = n_0, lx = lx) # nolint internal function
  # integrating the starting q_vector to t_p
  q_t <- mbd_solve(
    q_vector = q_i,
    matrix_a = matrix_a,
    time_interval = total_time,
    debug_mode = debug_mode
  )
  names(q_t) <- names(q_i)

  total_product <- q_t * one_over_cm * one_over_qm_binom
  missingspecies_min <- max(tips_interval[1] - n_0, 0)
  missingspecies_max <- min(tips_interval[2] - n_0, lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  pc <- sum(total_product[tips_components[1]:tips_components[2]])

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

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the pc. If \code{is.nan(log(pc))} the log-likelihood estimation
#'   by \link{mbd_loglik} is -Inf
#' @author Giovanni Laudanno
#' @noRd
calculate_conditional_prob2 <- function(
  brts,
  pars,
  cond,
  lx = 1000,
  n_0 = 2,
  tips_interval = c(n_0 * (cond > 0), Inf),
  methode = "lsodes",
  abstol = 1e-16,
  reltol = 1e-10,
  debug_mode = FALSE
) {
  check_cond(cond = cond, tips_interval = tips_interval, n_0 = n_0)
  if (cond == 0) {
    return(1)
  }
  pc0 <- pc1 <- pc2 <- 1
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
  names(a2_v1) <- paste0("Q", 0:lx)

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

  pc <- pc0 * (cond == 0) + pc1 * (cond == 1) + pc2 * (cond == 2)
  # small numerical errors can occur for short trees
  if (pc > 1 && pc < 1.01) {
    pc <- 1
  }
  if (!(pc >= 0 && pc <= 1)) {
    print(pc) # debug
    if (debug_mode == FALSE) {
      stop("problems: pc is wrong!") # debug
    }
  } # debug

  pc
}
