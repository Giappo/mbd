#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param age The age of the tree.
#' @param cond Set 1 if you want to condition on stem or crown age and non-extinction of the phylogeny. Set 0 otherwise.
#' @param soc Sets whether stem or crown age should be used (1 or 2).
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is lambda, the sympatric speciation rate;
#'   \item pars[2] is mu, the extinction rate;
#'   \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'   \item pars[4] is q, the single-lineage speciation probability.
#' }
#' @param tips_interval Sets tips boundaries constrain on simulated dataset.
#' @author Documentation by Giovanni Laudanno, use of this function by Richel J.C. Bilderbeek
#' @note This is an internal function, so it should be marked with
#'   \code{@noRd}. This is not done, as this will disallow all
#'   functions to find the documentation parameters
default_params_doc <- function(
  age,
  cond,
  soc,
  pars,
  tips_interval
) {
  # Nothing
}
