#' @author Giovanni Laudanno
#' @title Creates simulated trees under the multiple birth death process,
#'   including both sympatric and allopatric speciation
#' @description mbd_sim produces simulated trees allowing for
#'   three kind of events: sympatric speciation, multiple allopatric
#'   speciations and extinction.
#' @inheritParams default_params_doc
#' @examples
#' out <- mbd_sim(
#'   pars = c(0.6, 0.1, 0.4, 0.1),
#'   n_0 = 2,
#'   age = 10,
#'   cond = 1
#' )
#' graphics::plot(out$full_tree)
#' graphics::plot(out$reconstructed_tree)
#' out$l_matrix
#'
#' @export
mbd_sim <- function(
  pars,
  n_0 = 2,
  age = 10,
  cond = 1,
  seed = NA,
  tips_interval = c(n_0 * (cond > 0), Inf),
  brts_precision = 8
) {
  mbd::check_cond(cond = cond, tips_interval = tips_interval, n_0 = n_0)
  mbd::check_seed(seed = seed)
  if (
    mbd::check_pars(pars = pars) == "wrong"
  ) {
    stop("These parameters are wrong. Please check.")
  }
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  keep_the_sim <- 0
  while (keep_the_sim == 0) {
    total_count <- n_0 # all the species that ever appeared
    n_species <- n_0 # all the surviving species
    l_matrix <- mbd::initialize_l_matrix(age = age, n_0 = n_0)
    pool <- unique(l_matrix[, 3])[unique(l_matrix[, 3]) != 0]
    t <- age
    while (t > 0 && length(pool) > 0) {
      n_species <- length(pool)
      total_rate <- n_species * (lambda + mu) + nu
      if (total_rate > 0) {
        delta_t <- stats::rexp(1, rate = total_rate)
        outcome <- sample(
          c(-1, 1, 2),
          size = 1,
          prob = c(n_species * mu, n_species * lambda, nu)
        )
        delta_n <- (outcome == -1) * -1 +
          (outcome == 1) * 1 +
          (outcome == 2) * stats::rbinom(
            n = 1,
            size = n_species,
            prob = q
          )
        t <- t - delta_t

        if (delta_n > 0 & t > 0) {
          if (n_species > 1) {
            parents <- sample(pool, replace = FALSE, size = delta_n)
          } else {
            parents <- pool
          }
          new_interval <- (total_count + 1):(total_count + delta_n)
          new_ids <- abs(new_interval) * sign(parents)
          l_matrix[new_interval, 1] <- t
          l_matrix[new_interval, 2] <- parents
          l_matrix[new_interval, 3] <- new_ids
          pool <- c(pool, new_ids)
          total_count <- total_count + delta_n
        }
        if (delta_n < 0 & t > 0) {
          if (n_species > 1) {
            dead <- sample(pool, replace = FALSE, size = 1)
          } else {
            dead <- pool
          }
          l_matrix[abs(dead), 4] <- t
          pool <- pool[pool != dead]
        }
      } else {
        t <- 0
      }
    }
    l_matrix <- l_matrix[1:total_count, ]
    dim(l_matrix) <- c(total_count, 4)
    l_matrix[, 1] <- mbd::approximate_brts(
      brts = l_matrix[, 1],
      brts_precision = brts_precision
    )

    keep_the_sim <- mbd::evaluate_sim(
      l_matrix = l_matrix,
      cond = cond,
      n_0 = n_0,
      tips_interval = tips_interval,
      pool = pool
    )
  }
  info <- mbd::get_info_l_matrix(
    l_matrix = l_matrix,
    brts_precision = brts_precision
  )

  list(
    brts = info$brts,
    reconstructed_tree = info$reconstructed_tree,
    full_tree = info$full_tree,
    l_matrix = l_matrix
  )
}
