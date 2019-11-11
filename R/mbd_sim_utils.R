#' Extract information from the l-table
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
get_info_l_matrix <- function(
  l_matrix,
  brts_precision
) {
  n_survivors <- sum(l_matrix[, 4] == -1)
  if (n_survivors >= 2) {
    time_points <- unlist(unname(sort(
      DDD::L2brts(l_matrix, dropextinct = TRUE), decreasing = TRUE
    )))
    reconstructed_tree <- DDD::L2phylo(unname(l_matrix), dropextinct = TRUE)
  }
  if (n_survivors == 1) {
    time_points <- unname(l_matrix[which(l_matrix[, 4] == -1), 1])
    reconstructed_tree <- mbd::create_singleton_phylo(time_points)
  }
  if (n_survivors == 0) {
    time_points <- c()
    reconstructed_tree <- mbd::create_empty_phylo()
  }
  colnames(l_matrix) <- c("birth_time", "parent", "id", "death_time")
  brts <- sort(abs(as.numeric(time_points)), decreasing = TRUE)
  brts <- DDD::roundn(brts, digits = brts_precision)

  full_tree <- DDD::L2phylo(l_matrix, dropextinct = FALSE)
  list(
    brts = brts,
    reconstructed_tree = reconstructed_tree,
    full_tree = full_tree
  )
}

#' Establish whether the outcome of the simulation is in accordance with the
#'  required conditioning
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
evaluate_sim <- function(
  l_matrix,
  cond,
  n_0,
  tips_interval,
  pool
) {

  tips <- -1
  crown_species_dead <- cond

  # count extants species
  alive <- l_matrix[l_matrix[, 4] == -1, ]
  alive <- matrix(alive, ncol = 4)
  tips <- nrow(alive)
  testit::assert(length(pool) == tips)

  # tips check
  tips_condition <- tips >= tips_interval[1] & tips <= tips_interval[2]

  # survival of crown check
  crown_species_dead <- length(unique(sign(alive[, 3]))) != n_0
  crown_survival <- !crown_species_dead

  # should i keep this simulation?
  keep_the_sim <- (cond == 0) * 1 +
    (cond == 1) * (crown_survival && tips_condition)
  keep_the_sim
}

#' Initialize the l-table
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
initialize_l_matrix <- function(
  age,
  n_0
) {
  l_matrix <- matrix(0, nrow = 1e+06, 4)
  l_matrix[, 4] <- -1
  l_matrix[, 3] <- 0
  l_matrix[1, 1:4] <- c(age, 0, -1, -1)
  if (n_0 == 2) {
    l_matrix[2, 1:4] <- c(age, -1, 2, -1)
  }
  l_matrix
}
