#' Do an MBD maximum likelihood estimate
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
mbd_ml2 <- function(
  branching_times,
  init_param_values,
  fixed_params,
  estimated_params,
  init_n_species = 2,
  n_missing_species = 0,
  conditioned_on = "nothing"
) {

}