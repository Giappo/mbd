#' Is object \code{x} an MBD parameter set, 
#'   as can be created by \code{create_mbd_params}?
#' @return TRUE if yes, else FALSE 
#' @author Richel J.C. Bilderbeek
is_mbd_params <- function(x) {
  all(get_mbd_param_names() %in% names(x))
}