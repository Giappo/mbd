#' @title Maximization of the loglikelihood
#'   under a multiple birth-death diversification model
#' @description mbd_ml computes the maximum likelihood estimates of
#'   the parameters of a multiple birth-death diversification model
#'   for a given set of phylogenetic branching times.
#'   It also outputs the corresponding loglikelihood
#'   that can be used in model comparisons.
#'   Differently from mbd_ml it can account for three kind of events:
#'   sympatric (single) speciation, multiple (allopatric) speciation
#'   and extinction.
#' @author Giovanni Laudanno
#' @inheritParams default_params_doc
#' @return The output is a dataframe containing
#'   estimated parameters and maximum loglikelihood.
#'   The computed loglikelihood contains the factor q! m! / (q + m)!
#'   where q is the number of species in the phylogeny and m is the number of
#'   missing species, as explained in the supplementary
#'   material to Etienne et al. 2012.
#' @examples
#' set.seed(2)
#' lambda <- 0.2 # sympatric speciation rate
#' mu <- 0.15 # extinction rate;
#' nu <- 2.0 # multiple allopatric speciation trigger rate
#' q <- 0.1 # single-lineage speciation probability
#' sim_pars <- c(lambda, mu, nu, q)
#' crown_age <- 1
#' cond <- 1
#' n_0 <- 2
#' sim <- mbd_sim(
#'  pars = sim_pars,
#'  n_0 = n_0, # Use a crown age
#'  age = crown_age,
#'  cond = cond # Condition on non-extinction
#' )
#' start_pars <- c(0.2, 0.15, 1, 0.15)
#' optim_ids <- c(FALSE, FALSE, FALSE, TRUE)
#' graphics::plot(sim$reconstructed_tree)
#' out <- mbd::mbd_ml(
#'   start_pars = start_pars,
#'   true_pars = sim_pars,
#'   optim_ids = optim_ids,
#'   brts = sim$brts,
#'   cond = cond,
#'   n_0 = n_0,
#'   verbose = TRUE
#' )
#' @seealso use \link{pmb_ml} to perform a maximum likelihood estimation
#'   for a Pure Multiple Birth model.
#' @export
mbd_ml <- function(
  start_pars = c(0.3, 0.15, 2, 0.10),
  true_pars,
  optim_ids = c(TRUE, TRUE, TRUE, TRUE),
  brts,
  n_0 = 2,
  cond = 1,
  tips_interval = c(0, Inf),
  missnumspec = 0,
  lx = 1 + 2 * (length(brts) + length(missnumspec)),
  optimmethod = "simplex",
  methode = "expo",
  verbose = TRUE
) {
  namepars <- mbd::get_mbd_param_names()
  n_pars <- length(namepars)
  idparsopt <- (1:n_pars)[optim_ids]
  idparsfix <- (1:n_pars)[!optim_ids]
  parsfix <- true_pars[idparsfix]
  initparsopt <- start_pars[idparsopt]

  # res <- 10 * (1 + length(brts) + missnumspec) nolint
  tol <- c(1E-3, 1E-4, 1E-6)
  maxiter <- 1000 * round((1.25)^length(idparsopt))
  optimpars  <- c(tol, maxiter)
  changeloglikifnoconv <- FALSE

  if (!is.numeric(brts)) {
    stop("'brts' must be numeric")
  }
  if (length(initparsopt) != length(idparsopt)) {
    stop("lengths of 'idparsopt' and'initparsopt' must match")
  }
  if (length(parsfix) != length(idparsfix)) {
    stop("lengths of 'idparsfix' and'parsfix' must match")
  }
  if (!all(sort(c(idparsopt, idparsfix)) == 1:n_pars)) {
    stop(
      "IDs 1 to 4 must be present exactly once ",
      "in either 'idparsfix' or 'idparsopt'"
    )
  }
  if (length(idparsfix) == 0) {
    idparsfix <- NULL
  }
  if (missing(parsfix) && (length(idparsfix) == 0)) {
    parsfix <- idparsfix <- NULL
  }

  options(warn = -1)
  fail_pars <- rep(-1, n_pars)
  names(fail_pars) <- namepars
  fail_out <- data.frame(t(fail_pars), loglik = -1, df = -1, conv = -1)

  if (length(namepars[idparsopt]) == 0) {
    optstr <- "nothing"
  } else {
    optstr <- namepars[idparsopt]
  }
  if (verbose == TRUE) {
    cat("You are optimizing", optstr, "\n")
  }
  if (length(namepars[idparsfix]) == 0) {
    fixstr <- "nothing"
  } else {
    fixstr <- namepars[idparsfix]
  }
  if (verbose == TRUE) {
    cat("You are fixing", fixstr, "\n")
    cat("Optimizing the likelihood - this may take a while.", "\n")
    utils::flush.console()
  }

  trparsopt <- mbd_transform_forward(initparsopt) # nolint internal function
  trparsfix <- mbd_transform_forward(parsfix) # nolint internal function

  initloglik <- mbd_loglik_choosepar(
    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    brts = brts,
    n_0 = n_0,
    cond = cond,
    tips_interval = tips_interval,
    missnumspec = missnumspec,
    lx = lx,
    methode = methode
  )
  if (verbose == TRUE) {
    cat(
      "The loglikelihood for the initial parameter values is",
      initloglik, "\n"
    )
    utils::flush.console()
  }
  if (initloglik == -Inf) {
    warning(
      "The initial parameter values have a likelihood that is equal to 0",
      "or below machine precision. Try again with different initial values."
    )
    out2 <- fail_out
    return(invisible(out2))
  }
  if (verbose != TRUE) {
    if (rappdirs::app_dir()$os != "win") {
      sink("/dev/null")
    } else {
      sink(rappdirs::user_cache_dir())
    }
  }
  out <- DDD::optimizer(
    optimmethod = optimmethod,
    optimpars = optimpars,
    fun = mbd_loglik_choosepar,
    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    brts = brts,
    n_0 = n_0,
    cond = cond,
    tips_interval = tips_interval,
    missnumspec = missnumspec,
    lx = lx,
    methode = methode
  )
  if (verbose != TRUE) {
    sink() # Give back the output
  }
  if (out$conv != 0) {
    warning("Optimization has not converged. ",
      "Try again with different initial values."
    )
    out2 <- fail_out
    return(invisible(out2))
  }
  ml_pars <- mbd_transform_back(out$par) # nolint internal function
  ml_pars1 <- rep(0, n_pars); names(ml_pars1) <- namepars
  ml_pars1[idparsopt] <- ml_pars
  if (length(idparsfix) != 0) {
    ml_pars1[idparsfix] <- parsfix
  }
  max_lik <- as.numeric(unlist(out$fvalues))
  out2 <- data.frame(
    t(ml_pars1),
    loglik = max_lik,
    df = length(initparsopt),
    conv = unlist(out$conv)
  )

  tobeprint <- "Maximum likelihood parameter estimates:"
  for (ii in 1:n_pars) {
    tobeprint <- paste(
      tobeprint,
      paste0(names(ml_pars1[ii]), ":"),
      signif(ml_pars1[ii], digits = 3)
    )
  }
  if (verbose == TRUE) {
    s1 <- sprintf(tobeprint)
  }
  if (out2$conv != 0 & changeloglikifnoconv == TRUE) {
    out2$loglik <- -Inf
  }
  if (verbose == TRUE) {
    s2 <- sprintf("Maximum loglikelihood: %f", max_lik)
    cat("\n", s1, "\n", s2, "\n\n")
  }

  invisible(out2)
}
