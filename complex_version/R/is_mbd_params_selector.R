#' Determine if the supplied object in an MBD parameter selector,
#' as can be created by \code{create_mbd_params_selector}
#' @param x the object to be determined if it is an MBD
#'   parameters selector
#' @return TRUE if x is an MBD parameters selector
#' @examples
#'   s <- create_mbd_params_selector()
#'   testthat::expect_true(is_mbd_params_selector(s))
#'
#'   testthat::expect_false(is_mbd_params_selector("nonsense"))
#' @author Richel J.C. Bilderbeek
#' @export
is_mbd_params_selector <- function(x) {
  for (name in get_mbd_param_names()) {
    if (!name %in% names(x)) return(FALSE)
  }
  if (!is.logical(x$lambda)) return(FALSE)
  if (!is.logical(x$mu)) return(FALSE)
  if (!is.logical(x$nu)) return(FALSE)
  if (!is.logical(x$q)) return(FALSE)
  TRUE
}
