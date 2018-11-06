mbd_transform_forward <- function(pars) {
    pars <- as.numeric(unlist(pars))
    pars_transformed <- pars / (1 + pars)
    pars_transformed[which(pars == Inf)] <- 1
    pars_transformed
}

mbd_transform_back <- function(pars_transformed) {
    if (!all(pars_transformed <= 1)) {
        stop("inputs greater than 1 make no sense")
    }
    pars_transformed <- as.numeric(unlist(pars_transformed))
    pars <- pars_transformed / (1 - pars_transformed)
    pars
}

#' @title Transformed loglik function to plug into the mle routine
#' @description This function provides a likelihood for a subset of parameters.
#' This is built to work inside mbd_ml
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @export
mbd_loglik_choosepar <- function(
  trparsopt,
  trparsfix,
  idparsopt = 1:4,
  idparsfix = (1:4)[-idparsopt],
  brts,
  n_0 = 2,
  cond = 1,
  tips_interval = c(0, Inf),
  missnumspec = 0,
  lx = 1 + 2 * (length(brts) + length(missnumspec)),
  methode = "expo"
) {
  namepars <- mbd::get_mbd_param_names()
  n_pars <- length(namepars)
    if (length(trparsopt) == n_pars && missing(trparsfix)) {
        trparsfix <- NULL
    }
    trpars1 <- rep(0, n_pars)
    trpars1[idparsopt] <- trparsopt
    if (length(idparsfix) != 0) {
        trpars1[idparsfix] <- trparsfix
    }
    if (min(trpars1[1:n_pars]) < 0) {
        loglik <- -Inf
    } else {
        pars1 <- mbd_transform_back(trpars1)
        loglik <- mbd_loglik(
          pars = pars1,
          brts = brts,
          n_0 = n_0,
          cond = cond,
          tips_interval = tips_interval,
          missnumspec = missnumspec,
          lx = lx,
          methode = methode
        )
    }
    if (is.nan(loglik) || is.na(loglik)) {
      warning(
        "There are parameter values used which cause numerical problems:",
        trpars1, "\n")
        loglik <- -Inf
    }
    loglik
}
