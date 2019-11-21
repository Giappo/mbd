#' @title Pt
#' @author Giovanni Laudanno
#' @description Nee's function: pt
#' @inheritParams default_params_doc
#' @return pt
#' @export
p_t  <- function(lambda, mu, t) {
  time <- t
  exp_term <- exp(
    (mu - lambda) * time
  )
  out    <- (lambda == mu) * (1 / (1 + lambda * time)) +
    (lambda != mu) * (
      (lambda - mu + (lambda == mu)) /
        (lambda - mu * exp_term * (lambda != mu) + (lambda == mu))
    )
  unname(out)
}

#' @title 1 - Pt
#' @author Giovanni Laudanno
#' @description Nee's function: 1 - pt
#' @inheritParams default_params_doc
#' @return 1 - pt
#' @export
one_minus_pt  <- function(lambda, mu, t) {
  time <- t
  exp_term <- exp(
    (mu - lambda) * time
  )
  out    <- (lambda == mu) * (lambda * time / (1 + lambda * time)) +
    (lambda != mu) * (
      (mu - mu * exp_term + (lambda == mu)) /
        (lambda - mu * exp_term + (lambda == mu))
    )
  unname(out)
}

#' @title ut
#' @author Giovanni Laudanno
#' @description Nee's function: ut
#' @inheritParams default_params_doc
#' @return ut
#' @export
ut  <- function(lambda, mu, t) {
  time <- t
  exp_term <- exp(
    (mu - lambda) * time
  )
  out    <- (lambda == mu) * (lambda * time / (1 + lambda * time)) +
    (lambda != mu) * (
      (lambda - lambda * exp_term + (lambda == mu)) /
        (lambda - mu * exp_term * (lambda != mu) + (lambda == mu))
    )
  unname(out)
}

#' @title 1 - ut
#' @author Giovanni Laudanno
#' @description Nee's function: 1 - ut
#' @inheritParams default_params_doc
#' @return 1 - ut
#' @export
one_minus_ut  <- function(lambda, mu, t) {
  time <- t
  exp_term <- exp(
    (mu - lambda) * time
  )
  out    <- (lambda == mu) * (1 / (1 + lambda * time)) +
    (lambda != mu) * (
      (0 + (lambda == mu) + (lambda - mu) * exp_term) /
        (lambda - mu * exp_term + (lambda == mu))
    )
  unname(out)
}

#' @title Pn
#' @author Giovanni Laudanno
#' @description Nee's function: pn
#' @inheritParams default_params_doc
#' @param n number of species
#' @return pn
#' @export
pn <- function(lambda, mu, t, n) {
  out <- (n > 0) * mbd::p_t(t = t, lambda = lambda, mu = mu) *
    mbd::one_minus_ut(t = t, lambda = lambda, mu = mu) *
    mbd::ut(t = t, lambda = lambda, mu = mu) ^ (n - 1 + 2 * (n == 0)) +
    (n == 0) * (mbd::one_minus_pt(t = t, lambda = lambda, mu = mu))
  out
}

#' Calculate the approximate Nee et al's equivalent parameters, starting from
#'  the mbd parameters
#' @inheritParams default_params_doc
#' @export
get_nee_pars <- function(
  pars
) {
  nvec <- 1:100
  nee_lambda <- pars[1] + pars[3] * sum(pars[4] ^ nvec)
  nee_mu <- pars[2]
  c(nee_lambda, nee_mu)
}

#' Average number of species for b-d process.
#' See Kendall 1948a, pag. 12
#' @inheritParams default_params_doc
#' @export
nee_mean_nt <- function(
  nee_pars,
  n_0,
  age
) {
  lambda <- nee_pars[1]
  mu <- nee_pars[2]
  ll <- exp(age * (lambda - mu))
  average_n_t <- n_0 * ll
  average_n_t
}

#' Average number of species for b-d process.
#' See Kendall 1948a, pag. 12
#' @inheritParams default_params_doc
#' @export
nee_stdev_nt <- function(
  nee_pars,
  n_0,
  age
) {
  lambda <- nee_pars[1]
  mu <- nee_pars[2]
  ll <- exp(age * (lambda - mu))
  variance <- n_0 * ll * (ll - 1) * (lambda + mu) / (lambda - mu)
  stdev <- sqrt(variance)
  stdev
}
