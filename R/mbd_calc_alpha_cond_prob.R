#' The alpha conditional probability is ...
#' @inheritParams default_params_doc
#' @noRd
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
