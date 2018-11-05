# THESE SCRIPTS ARE NOT PART OF THE CORE MBD PACKAGE. THEY ARE ONLY EXPERIMENTAL: DON'T CALL THEM!

# Matrix Builders

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
create_a_no_mbd <- function(
  lambda,
  mu,
  nu,
  q,
  k,
  max_number_of_species,
  minimum_multiple_births
) {
  testit::assert(max_number_of_species < 2 ^ 31)
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

# This stuff might be useful

#' @title Internal mbd function
#' @description Internal mbd function.
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
determine_k_limit <- function(
  pars,
  brts,
  lx,
  soc,
  methode,
  abstol = 1e-16,
  reltol = 1e-10
) 
{
  lambda <- pars[1]
  nu <- pars[3]
  q <- pars[4]
  mvec <- 0:lx
  q_i <- c(1, rep(0, lx))
  total_time <- max(abs(brts));
  matrix_a <- create_a(
    pars = pars,
    k = soc,
    max_number_of_species = lx
  )
  p_m <- a_operator(
    q_vector = q_i,
    transition_matrix = matrix_a,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    a_abstol = abstol,
    a_reltol = reltol
  )
  soc + max(mvec[(mvec %in% which((cumsum(p_m / sum(p_m))) <= 0.95))])
}

# Alpha analysis

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
