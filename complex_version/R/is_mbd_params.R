#' Is object \code{x} an MBD parameter set,
#'   as can be created by \code{create_mbd_params}?
#' @param x the object to be determined of if
#'   its an MBD parameter set
#' @return TRUE if yes, else FALSE
#' @author Richel J.C. Bilderbeek
is_mbd_params <- function(x) {
  all(mbd::get_mbd_param_names() %in% names(x))
}
