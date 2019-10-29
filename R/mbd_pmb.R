#' Calculate the log-likelihood for a Pure Multiple Birth model
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @seealso use \link{mbd_loglik} to calculate the log-likelihood for
#'   a Multiple Birth Death model.
#' @export
pmb_loglik <- function(
  pars,
  brts,
  n_0 = 2
) {
# BASIC SETTINGS AND CHECKS
  check_brts(brts = brts, n_0 = n_0)
  if (
    check_pars(pars = pars, safety_checks = FALSE) == "wrong"
  ) {
    return(-Inf)
  }

  lambda <- pars[1]
  mu     <- pars[2]
  nu     <- pars[3]
  q      <- pars[4]

  if (mu != 0) stop("This function works only for mu = 0!")

  data <- brts2time_intervals_and_births(brts) # nolint internal function
  time_intervals <- data$time_intervals[-1]
  births <- data$births[-which(data$births == 0)]
  k <- n_0 + cumsum(c(0, births))
  a_term <- rep(1, length(time_intervals)) # branches
  b_term <- rep(1, length(time_intervals) - 1) # nodes

  # calculating branches contribution
  i <- 0:1e6
  for (t in seq_along(time_intervals)) {
    #(nu *(t_k-t_k-1))^i * exp(-nu * (t_k - t_k - 1)) / k!
    poisson_term <- stats::dpois(i, nu * time_intervals[t], log = FALSE)[
      stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0]
    ii <- i[stats::dpois(i, nu * time_intervals[t], log = FALSE) != 0]
    # (1) nu contribution: (1-q)^(k * i) * (nu *(t_k-t_k-1))^i
    #                      * exp(-nu *(t_k-t_k-1)) / k!
    # (2) lambda contribution: exp(-k * lambda *(t_k-t_k-1))
    a_term[t] <- sum(
      (1 - q) ^ (ii * k[t]) * poisson_term
    ) * # (1)
      exp(-k[t] * lambda * (time_intervals[t])) # (2)
  }
  #calculating nodes contribution
  # (1) nu contribution: nu *(k, b)* q^b *(1-q)^(k-b)
  # (2) lambda contribution: lambda * k (only if b==1)
  b_term <- (
    nu * choose(k[-length(k)], births) * q ^ births *  # (1)
    (1 - q) ^ (k[-length(k)] - births)                 # (1)
  ) + lambda * k[-length(k)] * (births == 1)           # (2)

  th_loglik <- sum(log(a_term)) + sum(log(b_term))
  th_loglik
}
