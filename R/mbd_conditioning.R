# cond_prob_q -----

#' Calculates the lambda component of dq
#' @inheritParams default_params_doc
#' @export
cond_prob_dq_lambda <- function(
  qq,
  qq2,
  k,
  m1,
  m2,
  mm
) {
  mm_minus_one <- mm - 1
  dq1 <- (2 * k + m1 - 1) * qq2[mm, mm_minus_one] +
    (2 * k + m2 - 1) * qq2[mm_minus_one, mm] -
    (2 * k + m1 + m2) * qq # ok
  dq1
}

#' Calculates the mu component of dq
#' @inheritParams default_params_doc
#' @export
cond_prob_dq_mu <- function(
  qq,
  qq2,
  k,
  m1,
  m2,
  mm
) {
  mm_plus_one <- mm + 1
  dq2 <- (m1 + 1) * qq2[mm, mm_plus_one] +
    (m2 + 1) * qq2[mm_plus_one, mm] -
    (2 * k + m1 + m2) * qq # ok
  dq2
}

#' Calculates the nu component of dq
#' @inheritParams default_params_doc
#' @export
cond_prob_dq_nu <- function(
  qq,
  nu_matrix
) {

  #N_{m1, n1} Q_{n1, n2} N^T_{n2, m2}
  dq3 <- nu_matrix %*% qq %*% t(nu_matrix) - qq
  dq3
}

#' Auxilary function for cond_prob, creating useful matrices
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno, Bart Haegeman
#' @export
cond_prob_q_matrices <- function(
  q,
  lx
) {
  nvec <- 1:lx - 1
  nu_matrix <- matrix(0, nrow = lx, ncol = lx)
  if (q != 0) {
    for (m1 in 0:(lx - 1)) {
      for (a1 in 0:floor((m1 + 1) / 2)) { # nolint lintrbot is math's enemy
        aux <- log(m1 + 1) +
          lfactorial(m1 - a1) -
          lfactorial(m1 - 2 * a1 + 1) -
          lfactorial(a1)
        aux <- aux + a1 * log(q) + (m1 + 1 - 2 * a1) * log(1 - q)
        nu_matrix[m1 + 1, m1 - a1 + 1] <- exp(aux)
      }
    }
  }
  rownames(nu_matrix) <- paste0("m1=", nvec)
  colnames(nu_matrix) <- paste0("n1=", nvec)

  empty_qq <- matrix(0, nrow = (lx + 2), ncol = (lx + 2))
  m1 <- col(nu_matrix) - 1
  m2 <- row(nu_matrix) - 1
  list(
    nu_matrix = nu_matrix,
    empty_qq = empty_qq,
    m1 = m1,
    m2 = m2
  )
}

#' Auxilary function for cond_prob_q, creating rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
cond_prob_q_rhs1 <- function(
  qvec,
  lambda,
  mu,
  nu,
  k,
  nu_matrix,
  m1,
  m2,
  empty_qq,
  t
) {
  lx2 <- length(qvec)
  lx <- sqrt(lx2)

  mm <- 2:(lx + 1)

  qq <- matrix(qvec, nrow = lx, ncol = lx)
  qq2 <- empty_qq
  qq2[mm, mm] <- qq

  dq_lambda <-
    cond_prob_dq_lambda(qq = qq, qq2 = qq2, k = k, m1 = m1, m2 = m2, mm = mm)
  dq_mu <-
    cond_prob_dq_mu(qq = qq, qq2 = qq2, k = k, m1 = m1, m2 = m2, mm = mm)
  dq_nu <-
    cond_prob_dq_nu(qq = qq, nu_matrix = nu_matrix)
  dq <- lambda * dq_lambda + mu * dq_mu + nu * dq_nu
  dq <- matrix(dq, nrow = lx2, ncol = 1)
  dq
}

#' Auxilary function for cond_prob_q, creating rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
cond_prob_q_rhs2 <- function(t, x, parms) {
  list(cond_prob_q_rhs1(
    qvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    k = parms$k,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_qq = parms$empty_qq,
    t = t
  ))
}

#' Solution of the integration of the ODE for Q_{m1, m2}
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
prob_cond_get_q_m1_m2 <- function(
  pars,
  brts,
  matrices,
  rhs_function = cond_prob_q_rhs2
) {
  tt <- max(abs(brts)) # time between crown age and present
  lx <- ncol(matrices$nu_matrix)

  parms <- list()
  parms$lambda <- pars[1]
  parms$mu <- pars[2]
  parms$nu <- pars[3]
  parms$kk <- 1
  parms$nu_matrix <- matrices$nu_matrix
  parms$m1 <- matrices$m1
  parms$m2 <- matrices$m2
  parms$empty_qq <- matrices$empty_qq
  q_0 <- c(y = c(1, rep(0, lx ^ 2 - 1)))

  ode_out <- mbd_solve(
    vector = q_0,
    func = rhs_function,
    time_interval = tt,
    parms = parms
  )
  q_m1_m2 <- matrix(ode_out, nrow = lx, ncol = lx)
  q_m1_m2
}

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the conditional probability
#' @author Giovanni Laudanno, Bart Haegeman
#' @export
cond_prob_q <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  lx = 30,
  tips_interval = c(n_0 * (cond > 0), Inf),
  debug_mode = FALSE
) {
  if (n_0 != 2) {
    stop("This works only for n_0 == 2.")
  }
  if (cond == 0) {
    return(1)
  }

  # construct auxiliary matrix
  matrices <- cond_prob_q_matrices(q = pars[4], lx = lx)

  # integrate equations
  q_m1_m2 <- prob_cond_get_q_m1_m2(
    pars = pars,
    brts = brts,
    matrices = matrices,
    rhs_function = cond_prob_q_rhs2
  )

  # compute conditioning probability
  m1 <- col(q_m1_m2) - 1
  m2 <- row(q_m1_m2) - 1
  p_m1_m2 <- q_m1_m2 / ((m1 + 1) * (m2 + 1)) # nolint lintr has issues with math
  pc <- sum(p_m1_m2)

  pc <- DDD::roundn(pc, digits = 8)
  if (!(pc >= 0 && pc <= 1)) {
    if (debug_mode == TRUE) {
      graphics::plot(
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

# cond_prob_p -----

#' Calculates the lambda component of dp
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
cond_prob_dp_lambda <- function(
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
cond_prob_dp_mu <- function(
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
cond_prob_dp_nu <- function(
  pp,
  nu_matrix
) {
  #N_{m1, n1} P_{n1, n2} N^T_{n2, m2}
  dp3 <- nu_matrix %*% pp %*% t(nu_matrix) - pp
  dp3
}

#' Auxilary function for cond_prob_p, computing useful matrices
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
cond_prob_p_matrices <- function(
  q,
  lx
) {
  nvec <- 1:lx - 1
  nu_matrix <- matrix(0, nrow = lx, ncol = lx)
  if (q != 0) {
    # alt: for (m1 in 0:(lx - 1)) {
    # alt: for (n1 in 0:m1) {
    for (n1 in nvec) {
      for (m1 in n1:(lx - 1)) {
        aux <- lchoose(n1, max(0, m1 - n1))
        aux <- aux + (m1 - n1) * log(q) + (2 * n1 - m1) * log(1 - q)
        nu_matrix[m1 + 1, n1 + 1] <- exp(aux)
      }
    }
  }
  rownames(nu_matrix) <- paste0("m1=", nvec)
  colnames(nu_matrix) <- paste0("n1=", nvec)

  empty_pp <- matrix(0, nrow = (lx + 2), ncol = (lx + 2))
  m1 <- col(nu_matrix) - 1
  m2 <- row(nu_matrix) - 1
  list(
    nu_matrix = nu_matrix,
    empty_pp = empty_pp,
    m1 = m1,
    m2 = m2
  )
}

#' Auxilary function for cond_prob_p, computing rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
cond_prob_p_rhs1 <- function(
  pvec,
  lambda,
  mu,
  nu,
  nu_matrix,
  m1,
  m2,
  empty_pp,
  t
) {
  lx2 <- length(pvec)
  lx <- sqrt(lx2)

  mm <- 2:(lx + 1)

  pp <- matrix(pvec, nrow = lx, ncol = lx)
  pp2 <- empty_pp
  pp2[mm, mm] <- pp

  dp_lambda <-
    cond_prob_dp_lambda(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_mu <-
    cond_prob_dp_mu(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_nu <-
    cond_prob_dp_nu(pp = pp, nu_matrix = nu_matrix)
  dp <- lambda * dp_lambda + mu * dp_mu + nu * dp_nu
  dp <- matrix(dp, nrow = lx2, ncol = 1)
  dp
}

#' Auxilary function for cond_prob_p, computing rhs
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
cond_prob_p_rhs2 <- function(t, x, parms) {
  list(cond_prob_p_rhs1(
    pvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_pp = parms$empty_pp,
    t = t
  ))
}

#' Solution of the integration of the ODE for P_{n1, n2}
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
prob_cond_get_p_m1_m2 <- function(
  pars,
  brts,
  matrices,
  rhs_function = cond_prob_p_rhs2
) {
  tt <- max(abs(brts)) # time between crown age and present
  lx <- ncol(matrices$nu_matrix)

  parms <- list()
  parms$lambda <- pars[1]
  parms$mu <- pars[2]
  parms$nu <- pars[3]
  parms$nu_matrix <- matrices$nu_matrix
  parms$m1 <- matrices$m1
  parms$m2 <- matrices$m2
  parms$empty_pp <- matrices$empty_pp
  p_0 <- matrix(0, nrow = lx, ncol = lx)
  p_0[2, 2] <- 1
  p_0 <- matrix(p_0, nrow = lx ^ 2, ncol = 1)

  ode_out <- mbd_solve(
    vector = p_0,
    func = rhs_function,
    time_interval = tt,
    parms = parms
  )

  p_m1_m2 <- matrix(ode_out, nrow = lx, ncol = lx)
  p_m1_m2
}

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @return the conditional probability
#' @author Giovanni Laudanno, Bart Haegeman
#' @inheritParams default_params_doc
#' @export
cond_prob_p <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  lx = 30,
  tips_interval = c(n_0 * (cond > 0), Inf),
  debug_mode = FALSE
) {
  if (n_0 != 2) {
    stop("This works only for n_0 == 2.")
  }
  if (cond == 0) {
    return(1)
  }
  rm(tips_interval)

  # construct auxiliary matrix
  matrices <- cond_prob_p_matrices(q = pars[4], lx = lx)

  # integrate equations
  p_m1_m2 <- prob_cond_get_p_m1_m2(
    pars = pars,
    brts = brts,
    matrices = matrices,
    rhs_function = cond_prob_p_rhs2
  )

  # compute conditioning probability
  pc <- 1 + p_m1_m2[1, 1] - sum(p_m1_m2[, 1]) - sum(p_m1_m2[1, ])

  pc <- DDD::roundn(pc, digits = 8)
  if (!(pc >= 0 && pc <= 1)) {
    if (debug_mode != TRUE) {
      stop("problems: pc is wrong!")
    }
  } # debug

  pc
}

# cond_prob_sim -----

#' Estimates conditional probability using simulations
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @export
cond_prob_sim <- function(
  pars,
  brts,
  cond = 1,
  n_0 = 2,
  lx = 30,
  tips_interval = c(n_0 * (cond > 0), Inf),
  n_sims = 1e4,
  saveit = TRUE
) {
  age <- max(brts)
  n_0 <- 2
  testit::assert(cond == 1)
  rm(tips_interval)

  # folder and files structure
  pars_filename <- get_pars_filename(
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
    sim <- mbd_sim(
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
        paste0(get_param_names(), " = ", pars, collapse = ", "),
        subtitle = paste0("crown_age = ", max(brts), ", n_sims = ", n_sims)
      )
    ggplot2::ggsave(filename = taxa_plot_filename, plot = taxa_plot)

  }
  pc_sim
}

##### cond_prob
#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the conditional probability
#' @author Giovanni Laudanno, Bart Haegeman
#' @export
cond_prob <- cond_prob_p

# cond_prob_p2 -----
# cond_prob_p2 defined with the zero term in nu

#' @noRd
cond_prob_dp_nu2 <- function(
  pp,
  nu_matrix
) {
  lx <- nrow(nu_matrix)
  diag_matrix <- matrix(0, nrow = lx, ncol = lx)
  diag(diag_matrix) <- diag(nu_matrix) # this is (1 - q) ^ m
  loss_term <-
    diag_matrix %*% pp %*% t(nu_matrix) +
    nu_matrix %*% pp %*% t(diag_matrix) -
    diag_matrix %*% pp %*% t(diag_matrix) -
    pp

  nu_matrix2 <- nu_matrix
  diag(nu_matrix2) <- c(
    1,
    rep(0, nrow(nu_matrix) - 1)
  )
  #N_{m1, n1} P_{n1, n2} N^T_{n2, m2}
  dp3 <- nu_matrix2 %*% pp %*% t(nu_matrix2) + loss_term
  dp3
}

#' Auxilary function for cond_prob_p
#' @author Giovanni Laudanno, Bart Haegeman
#' @noRd
cond_prob_p2_matrices <- function(
  q,
  lx
) {

  nu_matrix <- matrix(0, nrow = lx, ncol = lx)
  if (q != 0) {
    # alt: for (m1 in 0:(lx - 1)) {
    # alt: for (n1 in 0:m1) {
    for (n1 in 0:(lx - 1)) {
      for (m1 in n1:(lx - 1)) {
        aux <- lchoose(n1, max(0, m1 - n1))
        aux <- aux + (m1 - n1) * log(q) + (2 * n1 - m1) * log(1 - q)
        nu_matrix[m1 + 1, n1 + 1] <- exp(aux)
      }
    }
  }
  rownames(nu_matrix) <- paste0("m1=", 0:(lx - 1))
  colnames(nu_matrix) <- paste0("n1=", 0:(lx - 1))

  empty_pp <- matrix(0, nrow = (lx + 2), ncol = (lx + 2))
  m1 <- col(nu_matrix) - 1
  m2 <- row(nu_matrix) - 1
  list(
    nu_matrix = nu_matrix,
    empty_pp = empty_pp,
    m1 = m1,
    m2 = m2
  )
}

#' Auxilary function for cond_prob_p
#' @author Giovanni Laudanno, Bart Haegeman
#' @noRd
cond_prob_p2_rhs1 <- function(
  pvec,
  lambda,
  mu,
  nu,
  nu_matrix,
  m1,
  m2,
  empty_pp,
  t
) {
  lx2 <- length(pvec)
  lx <- sqrt(lx2)

  mm <- 2:(lx + 1)

  pp <- matrix(pvec, nrow = lx, ncol = lx)
  pp2 <- empty_pp
  pp2[mm, mm] <- pp

  dp_lambda <-
    cond_prob_dp_lambda(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_mu <-
    cond_prob_dp_mu(pp = pp, pp2 = pp2, m1 = m1, m2 = m2, mm = mm)
  dp_nu <-
    cond_prob_dp_nu2(pp = pp, nu_matrix = nu_matrix)
  dp <- lambda * dp_lambda + mu * dp_mu + nu * dp_nu
  dp <- matrix(dp, nrow = lx2, ncol = 1)
  dp
}

#' Auxilary function for cond_prob_p
#' @author Giovanni Laudanno, Bart Haegeman
#' @noRd
cond_prob_p2_rhs2 <- function(t, x, parms) {
  list(cond_prob_p2_rhs1(
    pvec = x,
    lambda = parms$lambda,
    mu = parms$mu,
    nu = parms$nu,
    nu_matrix = parms$nu_matrix,
    m1 = parms$m1,
    m2 = parms$m2,
    empty_pp = parms$empty_pp,
    t = t
  ))
}

#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @return the conditional probability
#' @author Giovanni Laudanno, Bart Haegeman
#' @noRd
cond_prob_p2 <- function(
  pars,
  brts,
  cond,
  n_0 = 2,
  lx = 30,
  tips_interval = c(n_0 * (cond > 0), Inf),
  debug_mode = FALSE
) {
  if (n_0 != 2) {
    stop("This works only for n_0 == 2.")
  }
  if (cond == 0) {
    return(1)
  }
  rm(tips_interval)
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  tt <- max(abs(brts)) # time between crown age and present
  times <- c(0, tt)

  # construct auxiliary matrix
  matrices <- cond_prob_p2_matrices(q = q, lx = lx)

  # integrate equations
  parms <- list()
  parms$lambda <- lambda
  parms$mu <- mu
  parms$nu <- nu
  parms$nu_matrix <- matrices$nu_matrix
  parms$m1 <- matrices$m1
  parms$m2 <- matrices$m2
  parms$empty_pp <- matrices$empty_pp
  p_0 <- matrix(0, nrow = lx, ncol = lx)
  p_0[2, 2] <- 1
  p_0 <- matrix(p_0, nrow = lx ^ 2, ncol = 1)

  ode_out <- deSolve::ode(
    y = p_0,
    times = times,
    func = cond_prob_p_rhs2,
    parms = parms,
    method = "lsoda",
    atol = 1e-100,
    rtol = 1e-6,
    tcrit = tt
  )[2, -1]
  p_m1_m2 <- matrix(ode_out, nrow = lx, ncol = lx)
  somma <- sum(p_m1_m2)
  testit::assert(somma >= 0.98)

  # compute conditioning probability
  pc <- 1 + p_m1_m2[1, 1] - sum(p_m1_m2[, 1]) - sum(p_m1_m2[1, ])

  if (!(pc >= 0 && pc <= 1)) {
    if (debug_mode != TRUE) {
      stop("problems: pc is wrong!")
    }
  } # debug

  pc
}
