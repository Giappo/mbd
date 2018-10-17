#' Do an MBD maximum likelihood estimate.
#' @param branching_times the branching times of the phylogeny
#' @param init_param_values initial parameter values,
#'   as created with \code{create_mbd_params}
#' @param fixed_params the parameters that are fixed.
#'   This can be any subset of
#'   set {\code{"lambda"}, \code{"mu"}, \code{"nu"} and \code{"q"}}
#' @param estimated_params the parameters that are estimated.
#'   This can be any subset of
#'   set {\code{"lambda"}, \code{"mu"}, \code{"nu"} and \code{"q"}}
#' @param init_n_species the number of species at the moment of the
#'   first branching time. Can be either one or two:
#'   \itemize{
#'     \item 1: stem
#'     \item 2: crown
#'   }
#' @param n_missing_species the number of species missing
#' @param conditioned_on on what must be the likelihood estimation be
#'   conditioned:
#'   \itemize{
#'     \item \code{"nothing"}: the branch lengths were obtained from
#'       a random phylogeny that could also have gone extinct
#'     \item \code{"non_extinction"}: the branch lengths were obtained from
#'       a phylogeny conditioned on non-extinction
#'   }
#' @author Richel J.C. Bilderbeek
#' @export
mbd_calc_max_lik <- function(
  branching_times,
  init_param_values,
  fixed_params,
  estimated_params,
  init_n_species = 2,
  n_missing_species = 0,
  conditioned_on = "nothing"
) {
  if (!is.numeric(branching_times)) {
    stop("'branching_times' must be numeric")
  }
  if (!all(branching_times >= 0.0)) {
    stop("All 'branching_times' must be positive")
  }
  if (!is_mbd_params(init_param_values)) {
    stop(
      "'init_param_values' must be an mbd_params, ",
      "as created by 'create_mbd_params'"
    )
  }

  if (!is_mbd_params_selector(fixed_params)) {
    stop(
      "'fixed_params' must be an MBD parameter selector, ",
      "as created by 'create_mbd_params_selector'"
    )
  }
  if (!is_mbd_params_selector(estimated_params)) {
    stop(
      "'estimated_params' must be an MBD parameter selector, ",
      "as created by 'create_mbd_params_selector'"
    )
  }
  lambda_once <- xor(fixed_params$lambda, estimated_params$lambda)
  mu_once <- xor(fixed_params$mu, estimated_params$mu)
  nu_once <- xor(fixed_params$nu, estimated_params$nu)
  q_once <- xor(fixed_params$q, estimated_params$q)
  if (!(lambda_once && mu_once && nu_once && q_once)) {
    stop(
      "'fixed_params' and 'estimated_params' together must select each ",
      "of the MBD parameters exactly once"
    )
  }
  if (init_n_species != 1 && init_n_species != 2) {
    stop("'init_n_species' must be 1 or 2")
  }
  if (n_missing_species < 0) {
    stop("'n_missing_species' must be positive")
  }
  if (!conditioned_on %in% c("nothing", "non_extinction")) {
    stop("'conditioned_on' must be either 'nothing' or 'non_extinction'")
  }

  # Convert data
  idparsopt <- NULL
  idparsfix <- NULL
  for (i in seq(1, 4)) {
    if (estimated_params[[i]]) {
      idparsopt <- c(idparsopt, i)
    } else {
      idparsfix <- c(idparsfix, i)
    }
  }
  values <- as.numeric(unlist(init_param_values))
  initparsopt <- values[idparsopt]
  parsfix <- values[idparsfix]

  conditioned_on_code <- 0
  if (conditioned_on == "non_extinction") {
    conditioned_on_code <- 1
  }

  testit::assert(length(idparsopt) + length(idparsfix) == 4)

  mbd_ml(
    brts = branching_times,
    initparsopt = initparsopt,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    parsfix = parsfix,
    missnumspec = n_missing_species,
    cond = conditioned_on_code,
    soc = init_n_species,
    verbose = FALSE
  )
}
