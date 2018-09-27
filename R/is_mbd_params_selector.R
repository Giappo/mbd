#' @param x the object to be determined if it is an MBD
#'   parameters selector
#' @return TRUE if x is an MBD parameters selector
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