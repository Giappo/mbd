# Build Parmsvec ----

#' Creates log nu matrix for the Q-equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_qeq_log_nu_mat <- function(lx) {
  nvec <- 1:lx - 1
  log_nu_mat <- matrix(-Inf, nrow = lx, ncol = lx)
  for (m1 in nvec) {
    for (a1 in 0:floor((m1 + 1) / 2)) { # nolint lintrbot is math's enemy
      log_nu_mat[m1 + 1, m1 - a1 + 1] <- log(m1 + 1) +
        lfactorial(m1 - a1) -
        lfactorial(m1 - 2 * a1 + 1) -
        lfactorial(a1)
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
  if (q != 0) {
    for (m1 in nvec) {
      for (a1 in 0:floor((m1 + 1) / 2)) { # nolint lintrbot is math's enemy
        log_q_mat[m1 + 1, m1 - a1 + 1] <-
          a1 * log(q) + (m1 + 1 - 2 * a1) * log(1 - q)
      }
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
      log_nu_mat[m1 + 1, n1 + 1] <- lchoose(n1, max(0, m1 - n1))
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
  if (q != 0) {
    for (m1 in nvec) {
      for (n1 in 0:m1) {
        log_q_mat[m1 + 1, n1 + 1] <-
          (m1 - n1) * log(q) + (2 * n1 - m1) * log(1 - q)
      }
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
  mbd::check_condprob_eq(eq = eq)
  if (eq == "p_eq" || eq == "sim") {
    log_nu_mat <- mbd::condprob_peq_log_nu_mat(lx)
  }
  if (eq == "q_eq") {
    log_nu_mat <- mbd::condprob_qeq_log_nu_mat(lx)
  }
  log_nu_mat
}

#' Creates log q matrix for the specified equation for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_log_q_mat <- function(lx, q, eq) {
  mbd::check_condprob_eq(eq = eq)
  if (eq == "p_eq" || eq == "sim") {
    log_q_mat <- mbd::condprob_peq_log_q_mat(lx = lx, q = q)
  }
  if (eq == "q_eq") {
    log_q_mat <- mbd::condprob_qeq_log_q_mat(lx = lx, q = q)
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
  fortran = TRUE
) {
  lx2 <- lx ^ 2
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]

  if (nu > 0 && q > 0) {
    log_nu_mat <- log_nu_mat[1:lx, 1:lx]
    log_q_mat <- log_q_mat[1:lx, 1:lx]
    nu_q_mat <- exp(log_nu_mat + log_q_mat)
  } else {
    nu_q_mat <- diag(rep(1, lx))
  }
  if (fortran == TRUE) {
    dim(nu_q_mat) <- c(lx2, 1)
    parmsvec <- c(
      lambda,
      mu,
      nu,
      nu_q_mat
    )
  } else {
    lambda <- pars[1]
    mu <- pars[2]
    nu <- pars[3]
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

#' Quickly builds parmsvec
#' @inheritParams default_params_doc
#' @export
create_fast_parmsvec <- function(
  pars,
  lx,
  eq,
  fortran
) {
  if (eq == "sim") {
    lx <- 10
    fortran <- FALSE
  }
  parmsvec <- mbd::condprob_parmsvec(
    pars = pars,
    log_nu_mat = mbd::condprob_log_nu_mat(lx = lx, eq = eq),
    log_q_mat = mbd::condprob_log_q_mat(lx = lx, q = pars[4], eq = eq),
    lx = lx,
    fortran = fortran
  )
  parmsvec
}

# Differentials in R ----

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
  lx2 <- length(pvec)
  lx <- sqrt(lx2)

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
  dim(dp) <- c(lx2, 1)
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
  lx2 <- length(qvec)
  lx <- sqrt(lx2)

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
  dq <- matrix(dq, nrow = lx2, ncol = 1)
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

# Compute probability distributions ----

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

  lx2 <- lx ^ 2
  age <- max(abs(brts)) # time between crown age and present

  # Define starting vector
  p_0 <- matrix(0, nrow = lx, ncol = lx)
  p_0[2, 2] <- 1
  dim(p_0) <- c(lx2, 1)

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

  lx2 <- lx ^ 2
  nu_q_mat <- parmsvec[3 + 1:lx2]
  dim(nu_q_mat) <- c(lx, lx)
  age <- max(abs(brts)) # time between crown age and present

  # Define starting vector
  q_0 <- matrix(0, nrow = lx, ncol = lx)
  q_0[1, 1] <- 1
  dim(q_0) <- c(lx2, 1)

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

# Conditional probability ----

#' Conditional probability calculated with the P-approach
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_p <- function(
  brts,
  parmsvec,
  fortran = TRUE,
  lx
) {
  if (fortran == TRUE) {
    rhs_function <- "mbd_runmodpcp"
  } else {
    rhs_function <- mbd::condprob_dp_rhs
  }
  p_n1_n2 <- mbd::condprob_p_n1_n2(
    brts = brts,
    parmsvec = parmsvec,
    rhs_function = rhs_function,
    lx = lx
  )
  pc <- 1 + p_n1_n2[1, 1] - sum(p_n1_n2[, 1]) - sum(p_n1_n2[1, ])
  pc <- DDD::roundn(pc, digits = 8)
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
  lx
) {
  if (fortran == TRUE) {
    rhs_function <- "mbd_runmodpcq"
  } else {
    rhs_function <- mbd::condprob_dq_rhs
  }
  q_m1_m2 <- mbd::condprob_q_m1_m2(
    brts = brts,
    parmsvec = parmsvec,
    rhs_function = rhs_function,
    lx = lx
  )
  pc <- sum(q_m1_m2)
  pc <- DDD::roundn(pc, digits = 8)
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
  saveit = TRUE
) {
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
  for (seed in 1:n_sims) {
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
        paste0(mbd::get_param_names(), " = ", pars, collapse = ", "),
        subtitle = paste0("crown_age = ", max(brts), ", n_sims = ", n_sims)
      )
    ggplot2::ggsave(filename = taxa_plot_filename, plot = taxa_plot)

  }
  pc <- DDD::roundn(pc_sim, digits = 8)
  mbd::check_pc(pc = pc)
  pc
}

#' Calculate conditional probability: loglik mode
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob <- function(
  brts,
  fortran = TRUE,
  lx,
  eq,
  parmsvec
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
      fortran = fortran,
      lx = lx
    ))
  }
  if (eq == "q_eq") {
    return(mbd::condprob_q(
      brts = brts,
      parmsvec = parmsvec,
      fortran = fortran,
      lx = lx
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
  eq,
  fortran = TRUE
) {
  mbd::condprob(
    brts = brts,
    fortran = fortran,
    lx = lx,
    eq = eq,
    parmsvec = mbd::create_fast_parmsvec(
      pars = pars,
      lx = lx,
      eq = eq,
      fortran = fortran
    )
  )
}

# Condprob selector ----

#' Selects the best eq for condprob
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
condprob_select_eq <- function(
  pars,
  fortran = TRUE
) {
  lx_stepsize <- 100
  lx_max <- 150
  threshold <- 0.02

  res <- list()
  for (eq in c("p_eq", "q_eq")) {
    lx <- 10
    log_nu_mat <- mbd::condprob_log_nu_mat(lx = lx_max, eq = eq)
    log_q_mat <- mbd::condprob_log_q_mat(lx = lx_max, q = pars[4], eq = eq)
    coeff <- 10
    pc1 <- mbd::condprob(
      brts = brts,
      fortran = fortran,
      lx = lx,
      eq = eq,
      parmsvec = mbd::condprob_parmsvec(
        pars = pars,
        log_nu_mat = log_nu_mat,
        log_q_mat = log_q_mat,
        lx = lx,
        fortran = fortran
      )
    )
    while(abs(coeff) > threshold) {
      delta_lx <- ceiling(lx_stepsize / lx)
      lx <- lx + delta_lx
      pc2 <- mbd::condprob(
        brts = brts,
        fortran = fortran,
        lx = lx,
        eq = eq,
        parmsvec = mbd::condprob_parmsvec(
          pars = pars,
          log_nu_mat = log_nu_mat,
          log_q_mat = log_q_mat,
          lx = lx,
          fortran = fortran
        )
      )
      coeff <- (pc2 - pc1) / delta_lx
      pc1 <- pc2
    }
    res[[eq]] <- list(
      coeff = coeff,
      lx1 = lx - delta_lx,
      lx2 = lx,
      y2 = pc2
    )
  }
  # calculate the intersection point of linear extrapolations
  p3 <- (res$q_eq$y2 * res$p_eq$coeff - res$p_eq$y2 * res$q_eq$coeff)
  p3 <- p3 / (res$p_eq$coeff - res$q_eq$coeff)
  testit::assert(p3 > res$q_eq$y2)
  testit::assert(p3 < res$p_eq$y2)
  if (p3 < 0.5) {
    eq <- "q_eq"
  } else {
    eq <- "p_eq"
  }
  eq
}
