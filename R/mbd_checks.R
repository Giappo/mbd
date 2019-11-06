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
  safety_checks = TRUE,
  lambda_limit = 10,
  mu_limit = 10,
  nu_limit = Inf,
  q_threshold = 1e-3
) {
  if (missing(q_threshold)) {
    q_threshold <- 0
  }

  # standard checks
  if (length(pars) != length(get_param_names())) {
    stop("'pars' must have a length of four.")
  }
  if (any(is.nan(pars))) {
    stop("'pars' cannot contain NaNs")
  }
  if (safety_checks == TRUE) {
    pars_mins <- c(0, 0, 0, 0 + q_threshold)
    pars_maxs <- c(lambda_limit, mu_limit, nu_limit, 1 - q_threshold)
  } else {
    pars_mins <- c(0, 0, 0, 0)
    pars_maxs <- c(Inf, Inf, Inf, 1 + .Machine$double.xmin)
  }
  right <- all(
    (pars >= pars_mins) & (pars < pars_maxs)
  )
  if (right) {
    out <- "right"
  } else {
    out <- "wrong"
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

#' @title Check q_t
#' @description Check q_t
#' @inheritParams default_params_doc
#' @return Nothing
#' @author Giovanni Laudanno
#' @export
check_q_vector <- function(
  q_t,
  t,
  pars,
  brts,
  debug_mode = FALSE
) {
  q_vector <- q_t[t, ]
  if (any(q_vector < 0)) {
    if (debug_mode == TRUE) {
      brts2 <- c(brts, 0)
      w <- data.frame(matrix(NA, nrow = length(q_vector), ncol = 0))
      w$values <- q_vector
      w$x <- seq_along(q_vector)
      w$cols <- ifelse(sign(q_vector) > 0, "blue", "red")
      print(w$values)
      graphics::plot(
        values ~ x,
        w,
        pch = 15,
        xlab = "m",
        ylab = "Q_m^k(t)",
        main = paste0(
          "pars = ", pars[1], ", ", pars[2], ", ", pars[3], ", ", pars[4]
        ),
        sub = paste0(
          "Problems in time interval: (", brts2[t - 1], ", ", brts2[t], ")"
        )
      )
    } else {
      stop("problems: q_t is negative!")
    }
  }
}

#' @title Check conditional probability
#' @description Check conditional probability
#' @inheritParams default_params_doc
#' @return Nothing
#' @author Giovanni Laudanno
#' @export
check_pc <- function(pc, debug_mode = FALSE) {
  if (!(pc >= 0 && pc <= 1)) {
    if (debug_mode == FALSE) {
      stop("problems: pc is wrong!")
    }
  }
  return()
}
