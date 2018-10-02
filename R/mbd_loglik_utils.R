#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
hyper_a_hanno <- function(n_species, k, q) {
  # HG function: fast O(N), updated after Moulis meeting
  #this is the matrix builder: helps to create A and B operators
  #it produces the structure
  # q ^ (m - n) * (1 - q) ^ (k + 2 * n-m) *
  #  sum_j 2 ^ j choose(k, j) * choose(n, m - n - j)
  j <- 0:k
  a_1 <- (1 - q) ^ (k) * choose(k, j) * (2)^j
  n_species <- n_species + 1
  matrix_a <- diag(a_1[1], nrow = n_species + 2, ncol = n_species + 2)
  matrix_a[1:(k + 1), 1] <- a_1
  for (dst in 2:n_species) {
    src <- dst - 1
    s <- src:min(n_species, 2 * src + k - 1)
    matrix_a[s + 2, dst] <- matrix_a[s, src] + matrix_a[s + 1, src]
    m <- s - 1
    n <- src - 1;
    matrix_a[s, src] <- matrix_a[s, src] * q ^ (m - n) * (1 - q) ^ (2 * n - m)
  }
  matrix_a[n_species, n_species] <- matrix_a[n_species, n_species] *
    (1 - q) ^ (n_species - 1);
  matrix_a[1:n_species, 1:n_species]
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
create_a_zero <- function(
  max_number_of_species,
  lambda,
  mu,
  q,
  k,
  matrix_builder = hyper_a_hanno
){
  testit::assert(max_number_of_species < Inf)
  nvec <- 0:max_number_of_species
  my_matrix <- lambda * matrix_builder(
    n_species = max_number_of_species, k = k, q = q
  )

  #new version to avoid the dumpster problem at the end of the matrix
  diag(my_matrix) <= (-lambda) * (1 - (1 - q) ^ (k + nvec)) - mu * (nvec + k)

  my_matrix[row(my_matrix) == col(my_matrix) - 1] <- mu *
    nvec[2:(max_number_of_species + 1)]
  my_matrix
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
create_b_zero <- function(
  max_number_of_species,
  q,
  k,
  b,
  matrix_builder = hyper_a_hanno
) {
  #lambda * choose(k, b) * q^b  is going to be added in logB in the main script
  k2 <- k - b
  matrix_builder(n_species = max_number_of_species, k = k2, q = q)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_a <- function(
  lambda,
  mu,
  nu,
  q,
  k,
  max_number_of_species
) {
  testit::assert(max_number_of_species < Inf)
  nvec <- 0:max_number_of_species
  m <- create_a_zero(max_number_of_species = max_number_of_species,
                       lambda = nu, mu = mu, q = q, k = k)
  m[row(m) == col(m) + 1] <- m[row(m) == col(m) + 1] +
    lambda * (nvec[1:(max_number_of_species)] + 2 * k)

  # new version to avoid the dumpster problem at the end of the matrix
  m[row(m) == col(m)] <- m[row(m) == col(m)] - lambda * (nvec + k)

  m
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_b <- function(
  lambda,
  nu,
  q,
  k,
  b,
  max_number_of_species
) {
  m <- create_b_zero(
    max_number_of_species = max_number_of_species,
    q = q,
    k = k,
    b = b
  )
  lambda * k * diag(max_number_of_species + 1) *
    (b == 1) + nu * choose(k, b) * (q^b) * m
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_a_no_mbd <- function(
  lambda,
  mu,
  nu,
  q,
  k,
  max_number_of_species,
  minimum_multiple_births
) {
  testit::assert(max_number_of_species < Inf)
  nvec <- 0:max_number_of_species
  m <- create_a_zero(
    max_number_of_species = max_number_of_species,
    lambda = nu, mu = mu, q = q, k = k
  )
  m[row(m) == col(m) + 1] <- m[row(m) == col(m) + 1] + lambda *
    (nvec[1:(max_number_of_species)] + 2 * k)
  m[row(m) == col(m)] <- m[row(m) == col(m)] -
    c(lambda * (nvec[-length(nvec)] + k), 0)
  m[row(m) > col(m) + minimum_multiple_births] <- 0
  m
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
create_b_no_mbd <- function(
  lambda,
  nu,
  q,
  k,
  b,
  max_number_of_species,
  minimum_multiple_births
) {
  m <- create_b_zero(
    max_number_of_species = max_number_of_species,
    q = q, k = k, b = b
  )
  b_matrix <- lambda * k * diag(max_number_of_species + 1) *
    (b == 1) + nu * choose(k, b) * (q^b) * m
  b_matrix[row(b_matrix) > col(b_matrix) + minimum_multiple_births] <- 0
  b_matrix
}

#' The A operator is given by the integration of a set of differential equations
#' between two consecutive nodes. So, defined the set in the time interval 
#' [t_{i-1}, t_i], where k species are present in the phylogeny, as:
#' 
#' d  
#' --Q^k_m(t) = SUM_n(M^k_m,n * Q^k_n(t) 
#' dt  
#' 
#' where m, n, label the amount of unseen species in the phylogeny,
#' A is thus defined as:
#' 
#' A(t_i - t_{i-1}) = exp(M(t_k - t_{k-1})
#' @inheritParams default_params_doc
#' @noRd
a_operator <- function(
  q_matrix,
  transition_matrix,
  time_interval,
  precision = 50L,
  a_abstol = 1e-16,
  a_reltol = 1e-10,
  methode = "expo"
) {
  precision_limit <- 2000
  precision_step1 <- 40
  precision_step2 <- 50
  max_repetitions <- 10
  result <- rep(-1, length(q_matrix))
  bad_result <- 0

  testit::assert(methode != "sexpm")
  if (methode == "expo") {
    result_nan <- result_negative <- 1
    repetition <- 1
    while ((result_nan == 1 | result_negative == 1) &
      repetition < max_repetitions
    ) {
      result <- try(
        expoRkit::expv(
          v = q_matrix,
          x = transition_matrix,
          t = time_interval,
          m = precision
        ),
        silent = TRUE
      )

      result_nan <- (any(!is.numeric(result)) || any(is.nan(result)))
      if (result_nan) {
        precision <- precision - precision_step1
      } else {
        result_negative <- (any(result < 0))
        if (result_negative) {
          precision <- precision + precision_step2
          if (precision > precision_limit) {
            break
          }
        }
      }
      repetition <- repetition + 1
    }
  }

  bad_result <- (any(!is.numeric(result)) || any(is.nan(result)))
  if (!bad_result) {
    bad_result <- (any(result < 0))
  }

  if (methode == "lsoda" | bad_result) {
    times <- c(0, time_interval)
    ode_matrix <- transition_matrix
    R.utils::withTimeout(result <- deSolve::ode(
      y = q_matrix,
      times = times,
      func = mbd_loglik_rhs,
      parms = ode_matrix,
      atol = a_abstol,
      rtol = a_reltol)[2, -1],
      timeout = 1001
    )
  }

  result
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
mbd_loglik_rhs <- function(
  t,
  x,
  pars
) {
  #builds right hand side of the ODE set for multiple birth model
  with(as.list(x), {
    starting_vector <- x
    transition_matrix <- pars
    dx <- rep(0, length(starting_vector))
    dx <- drop(transition_matrix %*% starting_vector)
    out <- (dx)
    names(out) <- names(x)
    return(list(out))
  })
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
determine_k_limit <- function(
  pars,
  brts,
  lx,
  soc,
  methode,
  abstol = 1e-16,
  reltol = 1e-10
) {
  lambda <- pars[1]
  nu <- pars[3]
  q <- pars[4]
  mvec <- 0:lx
  q_i <- c(1, rep(0, lx))
  total_time <- max(abs(brts));
  t_zero <- create_a(
    lambda = lambda,
    mu = 0,
    nu = nu,
    q = q,
    k = soc,
    max_number_of_species = lx
  )
  p_m <- a_operator(
    q_matrix = q_i,
    transition_matrix = t_zero,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    a_abstol = abstol,
    a_reltol = reltol
  )
  soc + max(mvec[(mvec %in% which((cumsum(p_m / sum(p_m))) <= 0.95))])
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calc_cond_prob <- function(
  brts,
  pars,
  lx = 1000,
  soc = 2,
  tips_interval = c(0, Inf),
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10
) {
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts));

  m <- 0:lx; length(m)
  one_over_cm <- (3 * (m + 1)) / (m + 3); length(one_over_cm)
  one_over_qm_binom <- 1 / choose((m + soc), soc)
  #starting with k = 0 and m = 2 missing species
  q_i <- rep(0, lx + 1);  q_i[3] <- 1

  t_matrix <- create_a(
    lambda = lambda, mu = mu, nu = nu, q = q, k = 0,
    max_number_of_species = lx
  )

  a2_v1 <- a_operator(
    q_matrix = q_i,
    transition_matrix = t_matrix,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    a_abstol = abstol,
    a_reltol = reltol
  )

  total_product <- a2_v1 * one_over_cm * one_over_qm_binom

  #these are the components I want to exclude
  # (the one corresponding to 0 and 1 tips)
  tips_components <- 1 + 0:1
  1 - sum(total_product[tips_components])
}


#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calc_cond_prob0 <- function(
  brts,
  pars,
  lx = 1000,
  soc = 2,
  tips_interval = c(0, Inf),
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10
) {
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts));

  m <- 0:lx; length(m)
  one_over_cm <- (3 * (m + 1)) / (m + 3); length(one_over_cm)
  one_over_qm_binom <- 1 / choose((m + soc), soc)
  q_i <- c(1, rep(0, lx)); length(q_i)

  t_matrix <- create_a(lambda = lambda, mu = mu, nu = nu, q = q, k = soc,
                       max_number_of_species = lx)

  a2_v1 <- a_operator(
    q_matrix = q_i,
    transition_matrix = t_matrix,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    a_abstol = abstol,
    a_reltol = reltol
  )
  total_product <- a2_v1 * one_over_cm * one_over_qm_binom
  missingspecies_min <- max((tips_interval[1] - 2), 0)
  missingspecies_max <- min((tips_interval[2] - 2), lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  sum(total_product[tips_components[1]:tips_components[2]])
}


#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calc_cond_prob_zero_p_b <- function(
  brts,
  pars,
  lx = 200,
  soc = 2,
  tips_interval = c(0, Inf),
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10
) {
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts))
  if (mu != 0) {
    cat("mu is supposed to be equal zero to use this function")
    return(NA)
  }

  m <- 0:lx; length(m)
  one_over_cm <- (3 * (m + 1)) / (m + 3); length(one_over_cm)
  one_over_qm_binom <- 1 / choose((m + soc), soc)
  q_i <- c(1, rep(0, lx)); length(q_i)

  t_matrix <- create_a(
    lambda = lambda,
    mu = mu,
    nu = nu,
    q = q,
    k = soc,
    max_number_of_species = lx
  )

  a2_v1 <- a_operator(
    q_matrix = q_i,
    transition_matrix = t_matrix,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    a_abstol = abstol,
    a_reltol = reltol
  )

  total_product <- a2_v1 * one_over_cm * one_over_qm_binom
  missingspecies_min <- max((tips_interval[1] - 2), 0)
  missingspecies_max <- min((tips_interval[2] - 2), lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  opposite_tips_components <- (m + 1)[
    !((m + 1) %in% (tips_components[1]:tips_components[2]))
  ]
  1 - sum(total_product[opposite_tips_components])
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @noRd
find_best_lx_for_pc <- function(
  brts,
  pars,
  soc = 2,
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10,
  iterations = 20,
  interval_min = 500,
  interval_max = 1400
) {
  a <- iterations / 2
  interval_width <- interval_max - interval_min
  step1 <- floor(interval_width / a)

  lx_test <- rep(
    NA,
    length(lxvec <- seq(interval_min + step1, interval_max - step1, step1))
  )
  i <- 1
  right_lx_coord <- 0
  for (lx2 in lxvec) {
    lx_test[i] <- mbd::calc_cond_prob0(
      brts = brts,
      pars = c(pars[1], 0, pars[3], pars[4]),
      lx = lx2,
      soc = soc,
      tips_interval = c(0, Inf),
      methode = methode,
      abstol = abstol,
      reltol = reltol
    )
    if (!is.na(abs(lx_test[i]))) {
      if (abs(lx_test[i] - 1) < 0.01) {
        right_lx_coord <- i
        lx <- lxvec[right_lx_coord]
        break
      }
    }
    i <- i + 1
  }
  if (right_lx_coord == 0) {
    right_lx_coord <- which(
      abs(lx_test - 1) == min(abs(lx_test - 1), na.rm = TRUE)
    )
    lx <- lxvec[right_lx_coord]
  }

  lx_test2 <- rep(
    NA,
    length(lxvec2 <- floor(seq(lx - step1, lx + step1, 2 * step1 / a)))
  )
  j <- 1
  right_lx_coord2 <- 0
  for (lx2 in lxvec) {
    lx_test2[j] <- mbd::calc_cond_prob0(
      brts = brts,
      pars = c(pars[1], 0, pars[3], pars[4]),
      lx = lx2,
      soc = soc,
      tips_interval = c(0, Inf),
      methode = methode,
      abstol = abstol,
      reltol = reltol
    )
    if (!is.na(abs(lx_test2[i]))) {
      if (abs(lx_test2[j] - 1) < 0.01) {
        right_lx_coord2 <- j
        lx <- lxvec2[right_lx_coord2]
        break
      }
    }
    j <- j + 1
  }
  if (right_lx_coord2 == 0) {
    right_lx_coord2 <- which(
      abs(lx_test2 - 1) == min(abs(lx_test2 - 1), na.rm = TRUE)
    )
    lx <- lxvec2[right_lx_coord2]
  }

  return(lx)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
calc_cond_prob1 <- function(
  brts,
  pars,
  soc = 2,
  tips_interval = c(0, Inf),
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10
){

  lx <- find_best_lx_for_pc(brts = brts, pars = pars, soc = soc) # nolint internal function
  testit::assert(pars[2] == 0)
  pc <- mbd::calc_cond_prob_zero_p_b(
    brts = brts,
    pars = pars,
    lx = lx,
    soc = soc,
    tips_interval = tips_interval,
    methode = methode,
    abstol = abstol,
    reltol = reltol
  )
  list(pc = pc, lx = lx)
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
alpha_conditional_probability <- function(
  brts, pars, alpha, tips_interval = c(0, Inf),
  cond = 1, soc = 2, methode = "expo",
  abstol = 1e-16, reltol = 1e-10,
  minimum_multiple_births = 0
) {
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  min_tips <- tips_interval[1]
  max_tips <- tips_interval[2]
  min_tips <- max(min_tips, soc * cond) #check this
  init_n_lineages <- soc
  total_time <- max(abs(brts));
  births <- c(0, brts2time_intervals_and_births(brts)$births) # nolint internal function
  k_interval <- init_n_lineages + cumsum(births)
  max_k <- max(k_interval)
  # alpha is the proportionality factor between max_k
  # and the edge of the matrix
  max_number_of_species <- alpha * max_k
  testit::assert(max_number_of_species < Inf)

  if (!(cond == 1 | tips_interval[1] > 0 | tips_interval[2] < Inf)) {
    pc <- 1; a2_v1 <- c(1, rep(0, max_number_of_species))
  } else {
    m <- 0:max_number_of_species;
    one_over_cm <- (3 * (m + 1)) / (m + 3)
    one_over_qm_binom <- 1 / choose((m + init_n_lineages), init_n_lineages)
    # applying tips constrain
    tips_components <- (1 + min_tips):(1 + min(max_tips, max_number_of_species))
    if (cond == 1) {
      # I am already considering the starting species to survive.
      # I must not double count them!
      tips_components <- tips_components - init_n_lineages
    }

    q_i <- c(1, rep(0, max_number_of_species))
    mk_n_zero <- create_a(
      lambda = lambda,
      mu = mu,
      nu = nu,
      q = q,
      k = soc,
      max_number_of_species = max_number_of_species
    )
    a2_v1 <- a_operator(
      q_matrix = q_i,
      transition_matrix = mk_n_zero,
      time_interval = total_time,
      precision = 50L,
      methode = methode,
      a_abstol = abstol,
      a_reltol = reltol
    )
    if (methode != "sexpm") {
      # it removes some small negative values
      # that can occurr as bugs from the integration process
      a2_v1 <- negatives_correction(a2_v1, pars) # nolint internal function
    }

    #adjust for the required minimum amount of mbd
    if (minimum_multiple_births > 0) {
      mk_n_zero_no_mbd <- create_a_no_mbd(
        lambda = lambda,
        mu = mu,
        nu = nu,
        q = q,
        k = soc,
        max_number_of_species = max_number_of_species,
        minimum_multiple_births = minimum_multiple_births
      )
      a2_v1_no_mbd <- a_operator(
        q_matrix = q_i,
        transition_matrix = mk_n_zero_no_mbd,
        time_interval = total_time,
        precision = 50L,
        methode = methode,
        a_abstol = abstol,
        a_reltol = reltol
      )
      a2_v1 <- a2_v1 - a2_v1_no_mbd
    }

    total_product <- a2_v1 * one_over_cm * one_over_qm_binom
    pc <- sum(total_product[tips_components])
  }
  return(list(pc = pc, a2_v1 = a2_v1))
}

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
alpha_analysis <- function(
  brts,
  pars,
  tips_interval,
  cond,
  soc,
  alpha0,
  max_k,
  methode = "expo",
  abstol,
  reltol,
  minimum_multiple_births
) {
  delta_alpha <- 1
  count <- 0
  same_result_count <- 0
  pc_notanumber <- 1
  alpha <- alpha0
  while (pc_notanumber) {
    pc_1 <- alpha_conditional_probability(
      brts = brts,
      pars = pars,
      tips_interval = tips_interval,
      cond = cond,
      soc = soc,
      alpha = alpha,
      methode = methode,
      abstol = abstol,
      reltol = reltol,
      minimum_multiple_births = minimum_multiple_births
    )$pc
    pc_notanumber <- is.nan(pc_1)
    alpha <- alpha - pc_notanumber
  }
  while (delta_alpha != 0 && count < 100 && same_result_count < 5) {
    pc_2 <- alpha_conditional_probability(
      brts = brts,
      pars = pars,
      tips_interval = tips_interval,
      cond = cond,
      soc = soc,
      alpha = alpha + delta_alpha,
      methode = methode,
      abstol = abstol,
      reltol = reltol,
      minimum_multiple_births = minimum_multiple_births
    )$pc
    if (is.nan(pc_2)) {
      delta_alpha <- delta_alpha - 1
    } else if (pc_2 < pc_1) {
      delta_alpha <- delta_alpha + 1
      same_result_count <- same_result_count + 1
    } else {
      same_result_count <- 0
      delta_pc <- abs(pc_2 - pc_1) / pc_1;
      delta_alpha <- floor(10 * delta_pc)
      alpha <- alpha + delta_alpha
      pc_1 <- pc_2
    }

    count <- count + 1
  }
  if (max_k * alpha >= 2000) {
    #check to see whether alpha is too big to be handled without memory issues
    alpha <- floor(1500 / max_k)
    pc_1 <- alpha_conditional_probability(
      brts = brts,
      pars = pars,
      tips_interval = tips_interval,
      cond = cond,
      soc = soc,
      alpha = alpha,
      methode = methode,
      abstol = abstol,
      reltol = reltol,
      minimum_multiple_births = minimum_multiple_births
    )$pc
  }
  pc <- pc_1
  if (count >= 100) {
    alpha <- 10
  }
  if (pc <= 0 | pc == Inf | pc == -Inf) {
    pc <- 1
    print("there's a problem with pc")
  }
  return(list(pc = pc, alpha = alpha))
}
