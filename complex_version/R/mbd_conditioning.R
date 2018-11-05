#' Called by \link{mbd_loglik} if and only if conditioned on non-extinction
#' @inheritParams default_params_doc
#' @return the pc. If \code{is.nan(log(pc))} the log-likelihood estimation
#'   by \link{mbd_loglik} failed
#' @author Giovanni Laudanno
#' @noRd
calc_cond_prob <- function(
brts,
pars,
lx = 1000,
soc = 2,
tips_interval = c(0, Inf),
methode = "expo",
abstol = 1e-16,
reltol = 1e-10
) 
{
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts));
  
  m <- 0:lx; length(m)
  one_over_cm <- (3 * (m + 1)) / (m + 3); length(one_over_cm)
  one_over_qm_binom <- 1 / choose((m + soc), soc)
  q_i <- c(1, rep(0, lx)); length(q_i)
  
  matrix_a <- create_a(
    pars = pars, 
    k = soc,
    max_number_of_species = lx
  )
  
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
calc_cond_prob2 <- function(
  brts,
  pars,
  lx = 1000,
  soc = 2,
  tips_interval = c(0, Inf),
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10
) {
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4];
  total_time <- max(abs(brts))
  m <- 0:lx
  one_over_cm <- (3 * (m + 1)) / (m + 3)
  one_over_qm_binom <- 1 / choose((m + soc), soc)
  #starting with k = 0 and m = 2 missing species
  q_i <- rep(0, lx + 1)
  q_i[3] <- 1
  
  matrix_a <- create_a(
    pars = pars, 
    k = 0,
    max_number_of_species = lx
  )
  
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
  
  matrix_a <- create_a(
    pars = pars,
    k = soc,
    max_number_of_species = lx
  )
  
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
  missingspecies_min <- max((tips_interval[1] - 2), 0)
  missingspecies_max <- min((tips_interval[2] - 2), lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  opposite_tips_components <- (m + 1)[
    !((m + 1) %in% (tips_components[1]:tips_components[2]))
    ]
  1 - sum(total_product[opposite_tips_components])
}

#' @title Conditional probability (alpha)
#' @description Conditional probability calculator for the alpha approach
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
mbd_calc_alpha_cond_prob <- function(
  brts,
  pars,
  alpha,
  tips_interval = c(0, Inf),
  cond = 1,
  soc = 2,
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10,
  minimum_multiple_births = 0
) 
{
  lambda <- pars[1]
  mu     <- pars[2]
  nu     <- pars[3]
  q      <- pars[4]
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
  if (max_number_of_species >= 2^31) {
    stop(
      "Maximum number of species exceeds R's memory limits. ",
      "max_number_of_species = ", max_number_of_species,
      ", alpha = ", alpha,
      ", max_k = ", max_k
    )
  }
  
  if (!(cond == 1 | tips_interval[1] > 0 | tips_interval[2] < Inf)) {
    pc <- 1; a2_v1 <- c(1, rep(0, max_number_of_species))
  } else {
    m <- 0:max_number_of_species
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
      pars = pars,
      k = soc,
      max_number_of_species = max_number_of_species
    )
    a2_v1 <- a_operator(
      q_vector = q_i,
      transition_matrix = mk_n_zero,
      time_interval = total_time,
      precision = 50L,
      methode = methode,
      abstol = abstol,
      reltol = reltol
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
        q_vector = q_i,
        transition_matrix = mk_n_zero_no_mbd,
        time_interval = total_time,
        precision = 50L,
        methode = methode,
        abstol = abstol,
        reltol = reltol
      )
      a2_v1 <- a2_v1 - a2_v1_no_mbd
    }
    
    total_product <- a2_v1 * one_over_cm * one_over_qm_binom
    pc <- sum(total_product[tips_components])
  }
  return(list(pc = pc, a2_v1 = a2_v1))
}
