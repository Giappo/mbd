#' Create the MBD parameters
#' @inheritParams default_params_doc
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
  if (q > 1.0) {
    stop("'q' cannot exceed unity. It is a probability!")
  }
  list(
    lambda = lambda,
    mu = mu,
    nu = nu,
    q = q
  )
}
