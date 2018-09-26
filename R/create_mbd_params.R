#' Create the MBD parameters
#' @param lambda the sympatric speciation rate
#' @param mu, the extinction rate
#' @param nu, the multiple allopatric speciation trigger rate
#' @param  q, the single-lineage speciation probability
#' @author Richel J.C. Bilderbeek
#' @export
create_mbd_params <- function(
    lambda = lambda,
    mu = mu,
    nu = nu,
    q = q
) {
  if (lambda < 0.0) {
    stop("'lambda' must be positive")
  }
  if (mu < 0.0) {
    stop("'mu' must be positive")
  }
  if (nu < 0.0) {
    stop("'nu' must be positive")
  }
  if (q < 0.0) {
    stop("'q' must be positive")
  }
  list(
    lambda = lambda,
    mu = mu,
    nu = nu,
    q = q
  )
}