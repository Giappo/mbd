#' @title Check 'n_0'
#' @author Giovanni Laudanno
#' @description Check 'n_0'
#' @inheritParams default_params_doc
#' @return nothing
#' @noRd
check_n_0 <- function(
  n_0
) {
  if (!is.numeric(n_0)) {
    stop("'n_0' must be numeric.")
  }
  if (!(n_0 == 1 | n_0 == 2)) {
    stop("'n_0' must be either 1 or 2.")
  }
  return()
}

#' @title Check mbd pars
#' @author Giovanni Laudanno
#' @description Check mbd pars
#' @inheritParams default_params_doc
#' @return Returns "wrong" if they are wrong, or "right" otherwise.
#' @noRd
check_pars <- function(
  pars,
  safety_threshold
) {
  lambda <- pars[1]
  mu <- pars[2]
  nu <- pars[3]
  q <- pars[4]
  if (missing(safety_threshold)) {
    safety_threshold <- 0
  }

  # checks
  if (length(pars) != length(mbd:::get_param_names())) {
    stop("'pars' must have a length of four.")
  }
  if (any(is.nan(pars))) {
    stop("'pars' cannot contain NaNs")
  }

  check <- (any(is.infinite(pars))) ||
    (lambda < 0) ||
    (mu < 0) ||
    (nu < 0) ||
    (q < 0) ||
    (q > 1) ||
    (q < 0 + safety_threshold) ||
    (q > 1 - safety_threshold)
  if (check == TRUE) {
    out <- "wrong"
  } else {
    out <- "right"
  }
  out
}

#' @title Check 'brts'
#' @author Giovanni Laudanno
#' @description Make sure that 'brts' are properly defined.
#' @inheritParams default_params_doc
#' @return nothing
#' @noRd
check_brts <- function(
  brts,
  n_0
) {
  check_n_0(n_0 = n_0)
  if (sum(brts == max(brts)) > 1) {
    stop("Crown/stem age has to be reported only once in the branching times.")
  }
  births <- brts2time_intervals_and_births(brts)$births # nolint internal function
  kvec <- n_0 + cumsum(c(0, births))
  kvec1 <- kvec[-c(1, length(kvec))]
  births1 <- births[-1]
  if (!all(births1 <= kvec1)) {
    stop("At any time you cannot have more speciations than number of species.")
  }
  if (!all(brts >= 0)) {
    stop("'brts' values must be positive. Present time is taken at the 0.")
  }
  return()
}

#' @title Check 'cond'
#' @author Giovanni Laudanno
#' @description Make sure that 'cond' and 'tips_interval' are coherent
#' @inheritParams default_params_doc
#' @return nothing
#' @noRd
check_cond <- function(
  cond,
  tips_interval,
  n_0
) {
  check_n_0(n_0 = n_0)
  if (!(cond %in% mbd_conds())) {
    stop("This conditioning is not implemented.")
  }
  if (tips_interval[2] < tips_interval[1]) {
    stop("'tips_interval' must contain two values, ",
         "of which the second is larger.")
  }
  if (any(tips_interval < 0)) {
    stop("'tips_interval' must contain two non-negative values.")
  }
  if (cond == 0 && (tips_interval[1] != 0 || tips_interval[2] < Inf)) {
    stop("If 'cond' == 0, you cannot put restrictions on tree tips.")
  }
  if (cond == 1 && tips_interval[1] < n_0) {
    stop("If 'cond' == 1, you cannot require less than 'n_0' tips.")
  }
  return()
}

#' @title Check 'seed'
#' @author Giovanni Laudanno
#' @description Check 'seed'
#' @inheritParams default_params_doc
#' @return nothing
#' @noRd
check_seed <- function(seed) {
  if (!is.na(seed)) {
    if (!is.numeric(seed)) {
      stop("'seed' must be integer or NA.")
    } else {
      if (seed %% 1 == 0) {
        set.seed(seed)
      } else {
        stop("'seed' must be integer or NA.")
      }
    }
  }
  return()
}
