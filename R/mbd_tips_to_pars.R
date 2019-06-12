#' @author Giovanni Laudanno
#' @title Calculate the difference in cumulative probabilities for a given
#'   percentile
#' @description Given the number of tips, it returns the difference between the
#'   cumulative probability to have less than that number of tips, and the
#'   required cumulative probability.
#' @inheritParams default_params_doc
#' @param percentile The area under the curve (aka probability mass) of the
#'   distribution up to the given number of tips.
#' @param tips The number of tips that you don't want to exceed when you
#'   simulate a tree (with a percent confidence of "percentile").
#' @return The parameters
#' @export
calc_diff_percentile <- function(
  pars,
  n_0 = 2,
  tips = 100,
  lx = 2 * tips,
  age = 10,
  percentile = 0.99,
  methode = "lsodes",
  abstol = 1e-16,
  reltol = 1e-10
) {
  q_i <- c(1, rep(0, lx))
  matrix_a <- create_a(pars = pars, k = n_0, lx = lx) # nolint internal function
  distribution <- a_operator( # integrating the starting q_vector to t_p
    q_vector = q_i,
    transition_matrix = matrix_a,
    time_interval = age,
    precision = 250L,
    methode = methode,
    abstol = abstol,
    reltol = reltol
  )
  distribution <- distribution / sum(distribution)
  # plot(distribution)
  cumulative <- cumsum(distribution)
  n_additional <- tips - n_0
  cum_prob_tips <- cumulative[n_additional + 1] # +1 is for the 0-th component
  diff <- abs(cum_prob_tips - percentile)
  diff
}

#' @author Giovanni Laudanno
#' @title Find the best parameters to have a fixed probability to simulate
#'   trees with less than the given number of tips
#' @description It returns the parameters that minimize the difference between
#'   the requested percentile and the cumulative probability for the given
#'   number of tips.
#' @inheritParams default_params_doc
#' @param percentile The area under the curve (aka probability mass) of the
#'   distribution up to the given number of tips.
#' @param tips The number of tips that you don't want to exceed when you
#'   simulate a tree (with a percent confidence of "percentile").
#' @return The parameters
#' @export
mbd_tips_to_pars_single <- function(
  start_pars = c(0.2, 0.15, 1.5, 0.15),
  tips = 100,
  percentile = 0.99,
  age = 10,
  n_0 = 2,
  optim_ids = c(FALSE, FALSE, TRUE, TRUE),
  methode = "ode45", #methode <- "lsodes"
  verbose = TRUE
) {

  optim_fun <- function(tr_optim_pars) {
    pars2 <- rep(0, length(start_pars))
    optim_pars <- pars_transform_back(tr_optim_pars) # nolint internal function
    pars2[optim_ids] <- optim_pars
    pars2[!optim_ids] <- start_pars[!optim_ids]
    
    if (
      are_these_parameters_wrong(
        brts = age,
        pars = pars2,
        safety_threshold = 1e-4,
        n_0 = n_0
        )
      ) {
      return(100) 
    }
    optim_out <- calc_diff_percentile(
      pars = pars2,
      n_0 = n_0,
      tips = tips,
      age = age,
      percentile = percentile
    )
    if (verbose == TRUE) {
      printed_values <- paste0(
        signif(c(pars2, optim_out), digits = 5)
      )
      print_this <- paste(printed_values, sep = ",")
      cat(print_this, "\n")
    }
    optim_out
  }
  
  tr_start_pars <- rep(0, length(start_pars))
  tr_start_pars <- pars_transform_forward(start_pars[optim_ids]) # nolint internal function
  out <- subplex::subplex(
    par = tr_start_pars,
    fn = function(x) optim_fun(x)
  )
  res <- out$par
  names(res) <- get_param_names()[optim_ids]
  res
}

#' @author Giovanni Laudanno
#' @title Returns best pars to have the process producing
#'   less than a given number of tips, with a given probability.
#' @description Same as "mbd_tips_to_pars_single", but it returns many possible
#'   initial conditions. It is useful in case of more than one parameter.
#' @inheritParams default_params_doc
#' @param percentile The area under the curve (aka probability mass) of the
#'   distribution up to the given number of tips.
#' @param tips The number of tips that you don't want to exceed when you
#'   simulate a tree (with a percent confidence of "percentile").
#' @param chain_length The number of different combination of starting
#'   parameters.
#' @return The parameters
#' @export
mbd_tips_to_pars <- function(
  tips = 100,
  age = 10,
  n_0 = 2,
  percentile = 0.99,
  start_pars = c(0.2, 0.15, 1.5, 0.15),
  optim_ids = c(FALSE, FALSE, TRUE, TRUE),
  chain_length = 10^sum(optim_ids),
  methode = "ode45", #methode <- "lsodes"
  verbose = TRUE
) {
  results <- vector("list", chain_length)
  for (i in 1:chain_length) {
    lambda <- stats::runif(n = 1, min = 0.18, max = 0.32)
    mu <- stats::runif(n = 1, min = 0, max = 0.12)
    nu <- stats::runif(n = 1, min = 0.95, max = 2.7)
    q <- stats::runif(n = 1, min = 0.09, max = 0.21)
    rand_pars <- c(lambda, mu, nu, q)
    start_pars2 <- start_pars
    start_pars2[optim_ids] <- rand_pars[optim_ids]
    results[[i]] <- mbd_tips_to_pars_single(
      tips = tips,
      age = age,
      n_0 = n_0,
      percentile = percentile,
      start_pars = start_pars2,
      optim_ids = optim_ids,
      methode = methode,
      verbose = verbose
    )
  }
  results
}

#' @author Giovanni Laudanno
#' @title Find nu-q relation
#' @description Find the exponent of the formula nu*q^x that minimizes the
#'  differences across all runs in the chain.
#' @inheritParams default_params_doc
#' @param tips The number of tips that you don't want to exceed when you
#'   simulate a tree 99 times out of 100.
#' @param chain_length The number of different combination of starting
#'   parameters.
#' @return The parameters
#' @export
find_nu_q_relation <- function(
  tips = 100,
  age = 10,
  start_pars = c(0.2, 0.15, 1, 0.15),
  chain_length = 40,
  verbose = FALSE
) {
  optim_ids <- c(FALSE, FALSE, TRUE, TRUE)
  nu_q_list <- mbd_tips_to_pars(
    tips = tips,
    age = age,
    start_pars = start_pars,
    optim_ids = optim_ids,
    chain_length = chain_length,
    verbose = verbose
  )
  # formula is: (nu * q ^ y) = const. We sample different values for y and we
  # select the one that makes the relation more consistent across the whole
  # nu-q list.
  var_y <- rep(NA, max_l <- 100)
  ys <- seq(from = (min_y <- 1.6), to = (max_y <- 2.6), by = (max_y - min_y) / (max_l - 1))
  for (i in 1:max_l) {
    x <- unlist(lapply(nu_q_list, FUN = function(x) x[1] * x[2]^ys[i]))
    var_y[i] <- stats::var(x)
  } # it is possible to visualize it with plot(pippo)
  best_y <- ys[which(var_y == min(var_y))]
  best_y
}

# # test functions
# if (1 == 2) {
#   out1 <- mbd_tips_to_pars_single(tips = 200, age = 10, start_pars = c(0.2,0.15,1,0.15), optim_ids = c(FALSE,FALSE,TRUE,TRUE))
#   out2 <- mbd_tips_to_pars_single(tips = 200, age = 10, start_pars = c(0.2,0.15,2,0.10), optim_ids = c(FALSE,FALSE,TRUE,TRUE))
#   profvis::profvis(out <- mbd_tips_to_pars_single(tips = 100, age = 10, start_pars = c(0.2,0.15,1,0.15), optim_ids = c(FALSE,FALSE,TRUE,TRUE),verbose = FALSE, transformation = TRUE)); out
#   out3 <- mbd_tips_to_pars(tips = 60, age = 9, start_pars = c(0.2,0.15,1,0.15), optim_ids = c(FALSE,FALSE,TRUE,TRUE), chain_length = 3)
#   out4 <- find_nu_q_relation(); out4
# }