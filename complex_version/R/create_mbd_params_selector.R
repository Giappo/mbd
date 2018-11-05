#' Create an MBD parameter selector
#' @inheritParams default_params_doc
#' @author Richel J.C. Bilderbeek
#' @export
create_mbd_params_selector <- function(
  lambda = FALSE,
  mu = FALSE,
  nu = FALSE,
  q = FALSE
) {
  if (lambda != TRUE && lambda != FALSE) {
    stop("'lambda' must be either TRUE or FALSE")
  }
  if (mu != TRUE && mu != FALSE) {
    stop("'mu' must be either TRUE or FALSE")
  }
  if (nu != TRUE && nu != FALSE) {
    stop("'nu' must be either TRUE or FALSE")
  }
  if (q != TRUE && q != FALSE) {
    stop("'q' must be either TRUE or FALSE")
  }
  list(lambda = lambda, mu = mu, nu = nu, q = q)
}
