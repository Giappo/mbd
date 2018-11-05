#' @author Giovanni Laudanno
#' @title Creates simulated trees under the multiple birth death process,
#'   including both sympatric and allopatric speciation
#' @description mbd_sim produces simulated trees allowing for
#'   three kind of events: sympatric speciation, multiple allopatric
#'   speciations and extinction.
#' @inheritParams default_params_doc
#' @examples
#' out <- mbd_sim(
#'   pars = c(0.6, 0.1, 0.4, 0.1), n_0 = 2, age = 10, cond = 1,
#'   tips_interval = c(0, Inf)
#' )
#' graphics::plot(out$tas)
#' graphics::plot(out$tes)
#' out$L
#'
#' @export
mbd_sim <- function(
  pars,
  n_0 = 2,
  age = 10,
  cond = 1,
  tips_interval = c(0, Inf),
  minimum_multiple_births = 0
) 
{
  if (length(pars) != 4) {
    stop("'pars' must have four parameters")
  }
  if (pars[1] < 0.0) {
    stop("The sympatric speciation rate 'pars[1]' must be positive")
  }
  if (pars[2] < 0.0) {
    stop("The extinction rate 'pars[2]' must be positive")
  }
  if (pars[3] < 0.0) {
    stop(
      "The multiple allopatric speciation trigger rate ",
      "'pars[3]' must be positive"
    )
  }
  if (pars[4] < 0.0) {
    stop("The single-lineage speciation probability 'pars[4]' must be positive")
  }
  if (tips_interval[2] < tips_interval[1]) {
    stop(
      "'tips_interval' must contain two values, ",
      "of which the second is larger"
    )
  }
  if (any(tips_interval < 0)) {
    stop("'tips_interval' must contain two positive values")
  }
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  init_n_lineages <- n_0
  tips <- -1; crown_species_dead <- cond; multiple_births_check <- 0;
  keep_the_sim <- 0
  while (keep_the_sim == 0 | multiple_births_check == 0) {
    total_count <- init_n_lineages
    pool <- 1:init_n_lineages
    while (total_count == init_n_lineages | length(pool) < init_n_lineages) {
      total_count <- init_n_lineages
      n_species <- init_n_lineages
      pool <- c(-1, 2)
      t <- age
      l_matrix <- matrix(0, nrow = 1e6, 4)
      l_matrix[, 4] <- -1
      l_matrix[, 3] <- 0
      l_matrix[1, 1:4] <- c(t, 0, -1, -1)
      l_matrix[2, 1:4] <- c(t, -1, 2, -1)
      while (t > 0) {
        n_species <- length(pool)
        total_rate <- n_species * (lambda + mu) + nu
        if (total_rate > 0) {
          delta_t <- stats::rexp(1, rate = total_rate)
          outcome <- sample(
            c(-1, 1, 2), size = 1,
            prob = c(n_species * mu, n_species * lambda, nu)
          )
          delta_n <- -1 * (outcome == -1) + 1 * (outcome == 1) +
            (outcome == 2) * stats::rbinom(n = 1, size = n_species, prob = q)
          t <- t - delta_t
          
          if (delta_n > 0 & t > 0) {
            if (n_species > 1) {
              parents <- sample(pool, replace = FALSE, size = delta_n)
            } else {
              parents <- pool
            }
            new_interval <- (total_count + 1):(total_count + delta_n)
            #-(delta_n:1)* 1e-5 add this if you need separate time points
            l_matrix[new_interval, 1] <- t
            l_matrix[new_interval, 2] <- parents
            l_matrix[new_interval, 3] <- abs(new_interval) * sign(parents)
            
            pool <- c(pool, abs(new_interval) * sign(parents))
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
    }
    l_matrix <- l_matrix[(1:total_count), ]
    #tips check
    tips <- length(l_matrix[, 4][l_matrix[, 4] == -1])
    #survival of crown check
    alive <- l_matrix[l_matrix[, 4] == -1, ]
    alive <- matrix(alive, ncol = 4)
    #if cond == 0 they will always look like they're alive, because I don't care
    crown_species_dead <- (length(unique(sign(alive[, 3]))) != 2) * cond
    #multiple births check
    births_rec_tree <- unlist(
      unname(
        sort(DDD::L2brts(l_matrix, dropextinct = TRUE), decreasing = TRUE)
      )
    )
    births_full_tree <- unlist(
      unname(
        sort(DDD::L2brts(l_matrix, dropextinct = FALSE), decreasing = TRUE)
      )
    )
    multiple_births_rec_tree <- sum(duplicated(births_rec_tree))
    multi_births_full_tree <- sum(duplicated(births_full_tree))
    #should i consider the full tree or the reconstructed one???
    multiple_births_check <- multiple_births_rec_tree >= minimum_multiple_births
    
    #should i keep this simulation?
    keep_the_sim <- (!crown_species_dead) &
      (tips >= tips_interval[1] &
         tips <= tips_interval[2]
      )
  }
  time_points <- unlist(
    unname(sort(DDD::L2brts(l_matrix, dropextinct = TRUE), decreasing = TRUE))
  )
  
  colnames(l_matrix) <- c("birth_time",
                          "parent",
                          "id",
                          "death_time")
  
  brts <- sort(abs(as.numeric(time_points)), decreasing = TRUE)
  reconstructed_tree <- DDD::L2phylo(unname(l_matrix), dropextinct = TRUE)
  full_tree <- DDD::L2phylo(l_matrix, dropextinct = FALSE)
  list(
    brts = brts,
    reconstructed_tree = reconstructed_tree,
    full_tree = full_tree,
    l_matrix = l_matrix,
    minimum_multiple_births = multi_births_full_tree
  )
}

# mbd_sim_dataset---------------------------------
# moved to razzo

# mbd_sim_dataset0---------------------------------
# moved to razzo
