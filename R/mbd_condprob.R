# parmsvec ----

#' Creates vector of parameters for the conditional probability calculated with
#' the P-equation
#' @return nu_matrix for the P equation
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_nu_matrix_p <- function(
  pars,
  lx,
  absorb
) {
  mbd::check_lx(lx)

  nu <- pars[3]
  q <- pars[4]
  nvec <- 1:lx - 1

  log_q_mat <- log_nu_mat <- matrix(-Inf, nrow = lx, ncol = lx)

  if (nu > 0 && q > 0) {
    for (m1 in nvec) {
      for (n1 in 0:m1) {
        # nu component
        log_nu_mat[m1 + 1, n1 + 1] <-
          lchoose(n1, max(0, m1 - n1))
        # q component
        log_q_mat[m1 + 1, n1 + 1] <-
          (m1 - n1) * log(q) +
          (2 * n1 - m1) * log(1 - q)
      }
    }
    nu_q_mat <- exp(log_nu_mat + log_q_mat)

    # absorbing state?
    if (absorb == TRUE) {
      mm <- lx - 1
      nvec <- 1:lx - 1
      nvec2 <- nvec[2 * nvec - mm >= 0]
      for (n in nvec2) {
        k <- 0:(2 * n - mm)
        nu_q_mat[lx, n + 1] <- sum(
          choose(n, mm - n + k) *
            q ^ (mm - n + k) *
            (1 - q) ^ (n - (mm - n + k))
        )
      }
      testit::assert(all(
        abs(colSums(nu_q_mat) - 1) < 1e-5
      ))
    }

  } else {
    nu_q_mat <- diag(rep(1, lx))
  }
  rownames(nu_q_mat) <- paste0("m1=", nvec)
  colnames(nu_q_mat) <- paste0("n1=", nvec)

  nu_q_mat
}

#' Creates vector of parameters for the conditional probability calculated with
#' the P-equation
#' @return nu_matrix for the P equation
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_nu_matrix_q <- function(
  pars,
  lx,
  absorb
) {
  mbd::check_lx(lx)

  nu <- pars[3]
  q <- pars[4]
  nvec <- 1:lx - 1

  log_q_mat <- log_nu_mat <- matrix(-Inf, nrow = lx, ncol = lx)

  if (nu > 0 && q > 0) {
    for (m1 in nvec) {
      for (a1 in 0:floor((m1 + 1) / 2)) { # nolint lintrbot is math's enemy
        # nu component
        log_nu_mat[m1 + 1, m1 - a1 + 1] <-
          log(m1 + 1) +
          lfactorial(m1 - a1) -
          lfactorial(m1 - 2 * a1 + 1) -
          lfactorial(a1)
        # q component
        log_q_mat[m1 + 1, m1 - a1 + 1] <-
          a1 * log(q) +
          (m1 + 1 - 2 * a1) * log(1 - q)
      }
    }
    nu_q_mat <- exp(log_nu_mat + log_q_mat)

    # absorbing state?
    if (absorb == TRUE) {
      mm <- lx - 1
      nvec2 <- nvec[2 * nvec - mm + 1 > 0]
      for (n in nvec2) {
        k <- 0:(2 * n - mm + 1)
        nu_q_mat[lx, n + 1] <- sum(
          (
            choose(n, mm - n + k) + 2 * choose(n, mm - n + k - 1)
          ) *
            q ^ (mm - n + k) *
            (1 - q) ^ (n + 1 - (mm - n + k))
        )
      }
    }
  } else {
    nu_q_mat <- diag(rep(1, lx))
  }
  rownames(nu_q_mat) <- paste0("m1=", nvec)
  colnames(nu_q_mat) <- paste0("n1=", nvec)

  nu_q_mat
}

#' Creates vector of parameters for the conditional probability
#' @return a vector composed of (in order): lambda, mu, nu, components of the
#'  nu-q matrix, components of the matrix of columns (m1_mat), components of the
#'  matrix of rows (m2_mat), components of the empty framed matrix (empty_mat).
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_parmsvec <- function(
  pars,
  eq,
  lx,
  absorb,
  fortran = TRUE
) {
  mbd::check_lx(lx)

  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]

  if (eq == "q_eq") {
    nu_q_mat <- mbd::condprob_nu_matrix_q(pars = pars, lx = lx, absorb = absorb)
  }
  if (eq == "p_eq") {
    nu_q_mat <- mbd::condprob_nu_matrix_p(pars = pars, lx = lx, absorb = absorb)
  }
  if (eq != "q_eq" && eq != "p_eq") {
    lx <- 10
    fortran <- FALSE
    nu_q_mat <- mbd::condprob_nu_matrix_p(pars = pars, lx = lx, absorb = absorb)
  }

  if (fortran == TRUE) {
    dim(nu_q_mat) <- c(lx ^ 2, 1)
    parmsvec <- c(
      lambda,
      mu,
      nu,
      nu_q_mat
    )
  } else {
    empty_mat <- matrix(0, nrow = (lx + 2), ncol = (lx + 2))
    m1 <- col(nu_q_mat) - 1
    m2 <- row(nu_q_mat) - 1
    parmsvec <- list(
      lambda = lambda,
      mu = mu,
      nu = nu,
      nu_matrix = nu_q_mat,
      empty_mat = empty_mat,
      m1 = m1,
      m2 = m2
    )
  }
  parmsvec
}

# condprob_p -----

#' Calculates the lambda component of dp
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_lambda <- function(
  pp,
  pp2,
  m1,
  m2,
  mm
) {
  mm_minus_one <- mm - 1
  dp1 <- (m1 - 1) * pp2[mm, mm_minus_one] +
    (m2 - 1) * pp2[mm_minus_one, mm] -
    (m1 + m2) * pp # ok
  dp1
}

#' Calculates the mu component of dp
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_mu <- function(
  pp,
  pp2,
  m1,
  m2,
  mm
) {
  mm_plus_one <- mm + 1
  dp2 <- (m1 + 1) * pp2[mm, mm_plus_one] +
    (m2 + 1) * pp2[mm_plus_one, mm] -
    (m1 + m2) * pp # ok
  dp2
}

#' Calculates the nu component of dp
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_nu <- function(
  pp,
  nu_matrix
) {
  dp3 <- nu_matrix %*% pp %*% t(nu_matrix) - pp
  dp3
}

#' Auxilary function for cond_prob_p, computing rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp <- function(
  pvec,
  lambda,
  mu,
  nu,
  nu_matrix,
  m1,
  m2,
  empty_mat,
  t
) {
  lx <- sqrt(length(pvec))

  mm <- 1 + 1:lx

  pp <- matrix(pvec, nrow = lx, ncol = lx)
  pp2 <- empty_mat
  pp2[mm, mm] <- pp

  dp_lambda <-
    mbd::condprob_dp_lambda(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_mu <-
    mbd::condprob_dp_mu(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_nu <-
    mbd::condprob_dp_nu(pp = pp, nu_matrix = nu_matrix)
  dp <- lambda * dp_lambda + mu * dp_mu + nu * dp_nu
  dim(dp) <- c(lx ^ 2, 1)
  dp
}

#' Auxilary function for cond_prob_p, computing rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_rhs <- function(t, x, parms) {
  list(mbd::condprob_dp(
    pvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_mat = parms$empty_mat,
    t = t
  ))
}

#' Calculates the lambda component of dp
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_absorb_lambda <- function(
  pp,
  pp2,
  m1,
  m2,
  mm
) {
  mm_minus_one <- mm - 1
  m1a <- m1
  m1a[, lx] <- 0
  m2a <- t(m1a)

  dp1 <- (m1 - 1) * pp2[mm, mm_minus_one] +
    (m2 - 1) * pp2[mm_minus_one, mm] -
    (m1a + m2a) * pp
  dp1
}

#' Calculates the mu component of dp
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_absorb_mu <- function(
  pp,
  pp2,
  m1,
  m2,
  mm
) {
  mm_plus_one <- mm + 1
  m1a <- m1
  m1a[, lx] <- 0
  m2a <- t(m1a)

  dp2 <- (m1 + 1) * pp2[mm, mm_plus_one] +
    (m2 + 1) * pp2[mm_plus_one, mm] -
    (m1a + m2a) * pp
  dp2
}

#' Auxilary function for cond_prob_p, computing rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_absorb <- function(
  pvec,
  lambda,
  mu,
  nu,
  nu_matrix,
  m1,
  m2,
  empty_mat,
  t
) {
  lx <- sqrt(length(pvec))

  mm <- 1 + 1:lx

  pp <- matrix(pvec, nrow = lx, ncol = lx)
  pp2 <- empty_mat
  pp2[mm, mm] <- pp

  dp_lambda <- mbd::condprob_dp_absorb_lambda(
    pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm
  )
  dp_mu <- mbd::condprob_dp_absorb_mu(
    pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm
  )
  dp_nu <- mbd::condprob_dp_nu(pp = pp, nu_matrix = nu_matrix)
  dp <- lambda * dp_lambda + mu * dp_mu + nu * dp_nu
  dim(dp) <- c(lx ^ 2, 1)
  dp
}

#' Auxilary function for cond_prob_p, computing rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dp_absorb_rhs <- function(t, x, parms) {
  list(mbd::condprob_dp_absorb(
    pvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_mat = parms$empty_mat,
    t = t
  ))
}

#' Conditional probability P(n1, n2) for all the states n1 and n2
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_p_n1_n2 <- function(
  brts,
  parmsvec,
  lx,
  rhs_function = mbd::condprob_dp_rhs
) {
  mbd::check_lx(lx)

  age <- max(abs(brts)) # time between crown age and present

  # Define starting vector
  p_0 <- matrix(0, nrow = lx, ncol = lx)
  p_0[2, 2] <- 1
  dim(p_0) <- c(lx ^ 2, 1)

  # Integrate
  p_out <- mbd::mbd_solve(
    vector = p_0,
    time_interval = age,
    func = rhs_function,
    parms = parmsvec
  )

  p_n1_n2 <- matrix(p_out, nrow = lx, ncol = lx)
  p_n1_n2
}

# condprob_q -----

#' Calculates the lambda component of dq
#' @inheritParams default_params_doc
#' @export
condprob_dq_lambda <- function(
  qq,
  qq2,
  m1,
  m2,
  mm
) {
  k <- 1
  mm_minus_one <- mm - 1
  dq1 <- (2 * k + m1 - 1) * qq2[mm, mm_minus_one] +
    (2 * k + m2 - 1) * qq2[mm_minus_one, mm] -
    (2 * k + m1 + m2) * qq # ok
  dq1
}

#' Calculates the mu component of dq
#' @inheritParams default_params_doc
#' @export
condprob_dq_mu <- function(
  qq,
  qq2,
  m1,
  m2,
  mm
) {
  k <- 1
  mm_plus_one <- mm + 1
  dq2 <- (m1 + 1) * qq2[mm, mm_plus_one] +
    (m2 + 1) * qq2[mm_plus_one, mm] -
    (2 * k + m1 + m2) * qq # ok
  dq2
}

#' Calculates the nu component of dq
#' @inheritParams default_params_doc
#' @export
condprob_dq_nu <- function(
  qq,
  nu_matrix
) {

  #N_{m1, n1} Q_{n1, n2} N^T_{n2, m2}
  dq3 <- nu_matrix %*% qq %*% t(nu_matrix) - qq
  dq3
}

#' Auxilary function for cond_prob_q, creating rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dq <- function(
  qvec,
  lambda,
  mu,
  nu,
  nu_matrix,
  m1,
  m2,
  empty_mat,
  t
) {
  lx <- sqrt(length(qvec))

  mm <- 2:(lx + 1)

  qq <- matrix(qvec, nrow = lx, ncol = lx)
  qq2 <- empty_mat
  qq2[mm, mm] <- qq

  dq_lambda <- mbd::condprob_dq_lambda(
    qq = qq, qq2 = qq2, m1 = m1, m2 = m2, mm = mm
  )
  dq_mu <- mbd::condprob_dq_mu(
    qq = qq, qq2 = qq2, m1 = m1, m2 = m2, mm = mm
  )
  dq_nu <- mbd::condprob_dq_nu(
    qq = qq, nu_matrix = nu_matrix
  )
  dq <- lambda * dq_lambda + mu * dq_mu + nu * dq_nu
  dq <- matrix(dq, nrow = lx ^ 2, ncol = 1)
  dq
}

#' Auxilary function for cond_prob_q, creating rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dq_rhs <- function(t, x, parms) {
  list(mbd::condprob_dq(
    qvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_mat = parms$empty_mat,
    t = t
  ))
}

#' Calculates the lambda component of dq
#' @inheritParams default_params_doc
#' @export
condprob_dq_absorb_lambda <- function(
  qq,
  qq2,
  m1,
  m2,
  mm
) {
  k <- 1
  mm_minus_one <- mm - 1
  m1a <- m1
  m1a[, lx] <- - k
  m2a <- t(m1a)

  dq1 <- (2 * k + m1 - 1) * qq2[mm, mm_minus_one] +
    (2 * k + m2 - 1) * qq2[mm_minus_one, mm] -
    (2 * k + m1a + m2a) * qq # ok
  dq1
}

#' Calculates the mu component of dq
#' @inheritParams default_params_doc
#' @export
condprob_dq_absorb_mu <- function(
  qq,
  qq2,
  m1,
  m2,
  mm
) {
  k <- 1
  mm_plus_one <- mm + 1
  m1a <- m1
  m1a[, lx] <- - k
  m2a <- t(m1a)

  dq2 <- (m1 + 1) * qq2[mm, mm_plus_one] +
    (m2 + 1) * qq2[mm_plus_one, mm] -
    (2 * k + m1a + m2a) * qq # ok
  dq2
}

#' Auxilary function for cond_prob_q, creating rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dq_absorb <- function(
  qvec,
  lambda,
  mu,
  nu,
  nu_matrix,
  m1,
  m2,
  empty_mat,
  t
) {
  lx <- sqrt(length(qvec))

  mm <- 2:(lx + 1)

  qq <- matrix(qvec, nrow = lx, ncol = lx)
  qq2 <- empty_mat
  qq2[mm, mm] <- qq

  dq_lambda <- mbd::condprob_dq_lambda(
    qq = qq, qq2 = qq2, m1 = m1, m2 = m2, mm = mm
  )
  dq_mu <- mbd::condprob_dq_mu(
    qq = qq, qq2 = qq2, m1 = m1, m2 = m2, mm = mm
  )
  dq_nu <- mbd::condprob_dq_nu(
    qq = qq, nu_matrix = nu_matrix
  )
  dq <- lambda * dq_lambda + mu * dq_mu + nu * dq_nu
  dq <- matrix(dq, nrow = lx ^ 2, ncol = 1)
  dq
}

#' Auxilary function for cond_prob_q, creating rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
condprob_dq_absorb_rhs <- function(t, x, parms) {
  list(mbd::condprob_dq_absorb(
    qvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_mat = parms$empty_mat,
    t = t
  ))
}

#' Conditional probability Q(m1, m2) for all the states m1 and m2
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_q_m1_m2 <- function(
  brts,
  parmsvec,
  lx,
  rhs_function = mbd::condprob_dq_rhs
) {
  mbd::check_lx(lx)

  nu_q_mat <- parmsvec[3 + 1:(lx ^ 2)]
  dim(nu_q_mat) <- c(lx, lx)
  age <- max(abs(brts)) # time between crown age and present

  # Define starting vector
  q_0 <- matrix(0, nrow = lx, ncol = lx)
  q_0[1, 1] <- 1
  dim(q_0) <- c(lx ^ 2, 1)

  # Integrate
  q_out <- mbd::mbd_solve(
    vector = q_0,
    time_interval = age,
    func = rhs_function,
    parms = parmsvec
  )

  q_m1_m2 <- matrix(q_out, nrow = lx, ncol = lx)
  m1 <- col(nu_q_mat) - 1
  m2 <- t(m1)
  p_m1_m2 <- q_m1_m2 / ((m1 + 1) * (m2 + 1)) # nolint lintr has issues with math
  p_m1_m2
}

# conditional probability ----

#' Conditional probability calculated with the P-approach
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_p <- function(
  brts,
  parmsvec,
  fortran = TRUE,
  lx,
  absorb
) {
  mbd::check_lx(lx)

  if (fortran == TRUE) {
    rhs_function <- "mbd_runmodpcp"
  } else {
    if (absorb == FALSE) {
      rhs_function <- mbd::condprob_dp_rhs
    } else {
      rhs_function <- mbd::condprob_dp_absorb_rhs
    }
  }
  p_n1_n2 <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    rhs_function = rhs_function,
    lx = lx
  )
  pc <- 1 + p_n1_n2[1, 1] - sum(p_n1_n2[, 1]) - sum(p_n1_n2[1, ])
  mbd::check_pc(pc = pc)
  pc
}

#' Conditional probability calculated with the Q-approach
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_q <- function(
  brts,
  parmsvec,
  fortran = TRUE,
  lx,
  absorb
) {
  mbd::check_lx(lx)

  if (fortran == TRUE) {
    rhs_function <- "mbd_runmodpcq"
  } else {
    if (absorb == FALSE) {
      rhs_function <- mbd::condprob_dq_rhs
    } else {
      rhs_function <- mbd::condprob_dq_absorb_rhs
    }
  }
  q_m1_m2 <- mbd::condprob_q_m1_m2(
    brts = brts,
    parmsvec = parmsvec,
    rhs_function = rhs_function,
    lx = lx,
    absorb = absorb
  )
  pc <- sum(q_m1_m2)
  mbd::check_pc(pc = pc)
  pc
}

#' Estimates conditional probability using simulations
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @export
condprob_sim <- function(
  brts,
  parmsvec,
  lx = 1e4,
  saveit = TRUE,
  starting_seed = 1
) {
  mbd::check_lx(lx)

  # pars
  n_sims <- lx
  age <- max(brts)
  n_0 <- 2
  lambda <- parmsvec$lambda
  mu <- parmsvec$mu
  nu <- parmsvec$nu
  q <- 1 - parmsvec$nu_matrix[2, 2]
  pars <- c(lambda, mu, nu, q)

  # folder and files structure
  pars_filename <- mbd::get_pars_filename(
    pars = pars,
    age = age
  )
  sim_filename <- paste0(pars_filename, "-pc_sim.Rdata")
  taxa_plot_filename <- paste0(pars_filename, "-taxa_plot.png")
  delete_file <- FALSE
  if (file.exists(sim_filename)) {
    load(sim_filename)
    if (measure$n_sims >= n_sims) {
      return(measure$pc_sim)
    } else {
      delete_file <- TRUE
    }
  }

  # calculate pc_sim
  score <- 0
  n_tips <- rep(0, n_sims)
  seed_interval <- starting_seed:(starting_seed + n_sims - 1)
  for (seed in seed_interval) {
    sim <- mbd::mbd_sim(
      pars = pars,
      n_0 = n_0,
      cond = 0,
      age = age,
      seed = seed
    )

    l_matrix <- sim$l_matrix
    alive <- l_matrix[l_matrix[, 4] == -1, ]
    alive <- matrix(alive, ncol = 4)
    n_tips[seed] <- nrow(alive)
    crown_species_dead <- (length(unique(sign(alive[, 3]))) != n_0)
    crown_survival <- !crown_species_dead

    score <- score + 1 * crown_survival
  }
  pc_sim <- score / n_sims

  if (saveit == TRUE) {
    measure <- list(
      pc_sim = pc_sim,
      n_sims = n_sims,
      n_tips = n_tips
    )

    if (delete_file == TRUE) {
      file.remove(sim_filename)
      file.remove(taxa_plot_filename)
    }

    # save pc_sim
    save(
      measure,
      file = sim_filename
    )

    # save plot taxa
    taxa_plot <- ggplot2::ggplot(
      data = data.frame(measure), ggplot2::aes(x = n_tips)
    ) +
      ggplot2::geom_histogram(bins = 30) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(
        paste0(
          mbd::get_param_names(),
          " = ",
          mbd::to_string2(signif(pars, digits = 3)),
          collapse = ", "
        ),
        subtitle = paste0("crown_age = ", max(brts), ", n_sims = ", n_sims)
      )
    ggplot2::ggsave(filename = taxa_plot_filename, plot = taxa_plot)

  }
  mbd::check_pc(pc = pc_sim)
  pc_sim
}

#' Conditional probability using Nee et al's approach
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_nee <- function(
  brts,
  parmsvec
) {

  n_0 <- 2
  lambda <- parmsvec$lambda
  mu <- parmsvec$mu
  nu <- parmsvec$nu
  q <- 1 - parmsvec$nu_matrix[2, 2]
  pars <- c(lambda, mu, nu, q)

  nee_pars <- mbd::get_nee_pars(pars = pars)
  age <- max(abs(brts))

  pc <- mbd::p_t(
    lambda = nee_pars[1],
    mu = nee_pars[2],
    t = age
  ) ^ n_0
  pc
}

#' Calculate conditional probability: loglik mode
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob <- function(
  brts,
  parmsvec,
  lx,
  eq,
  absorb = TRUE,
  fortran = TRUE
) {

  if (eq == "sim") {
    return(mbd::condprob_sim(
      brts = brts,
      parmsvec = parmsvec,
      lx = lx
    ))
  }
  if (eq == "p_eq") {
    return(mbd::condprob_p(
      brts = brts,
      parmsvec = parmsvec,
      lx = lx,
      absorb = absorb,
      fortran = fortran
    ))
  }
  if (eq == "q_eq") {
    return(mbd::condprob_q(
      brts = brts,
      parmsvec = parmsvec,
      lx = lx,
      absorb = absorb,
      fortran = fortran
    ))
  }
  if (eq == "nee") {
    return(mbd::condprob_nee(
      brts = brts,
      parmsvec = parmsvec
    ))
  }
  stop("'eq' is not valid")
}

#' Calculate conditional probability: user mode
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
calculate_condprob <- function(
  pars,
  brts,
  lx,
  eq = "nee",
  absorb = TRUE,
  fortran = TRUE
) {

  if (eq == "nee") {
    lx <- 10
  }

  mbd::condprob(
    brts = brts,
    fortran = fortran,
    lx = lx,
    eq = eq,
    parmsvec = mbd::condprob_parmsvec(
      pars = pars,
      eq = eq,
      lx = lx,
      absorb = absorb,
      fortran = fortran
    )
  )
}

# condprob utils ----

#' Maximum allowed value for lx in condprob
#' @export
max_lx_condprob <- function() {
  max_lx <- 70
  max_lx
}

#' Minimum allowed value for lx in condprob
#' @export
min_lx_condprob <- function() {
  min_lx <- 20
  min_lx
}

#' Fix the value for lx for condprob in \link{mbd_loglik}
#' @inheritParams default_params_doc
#' @export
get_lx_condprob <- function(pars, n_0, age, lx) {
  rm(pars)
  rm(n_0)
  rm(age)

  lx_condprob <- max(
    mbd::min_lx_condprob(),
    min(
      mbd::max_lx_condprob(),
      ceiling(lx / 2)
    )
  )
  lx_condprob
}
