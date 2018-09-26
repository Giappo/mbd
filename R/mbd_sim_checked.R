#' More intuitive interface of \code{mbd_sim}
#' @param mbd_params MBD parameters, 
#'   as created by \code{create_mbd_params}
#' @param crown_age the crown age. 
#'   Either \code{crown_age} or \code{stem_age} must be exclusively set    
#' @param stem_age the stem age    
#'   Either \code{crown_age} or \code{stem_age} must be exclusively set
#' @param conditioned_on conditioning of the simulation, can be:
#'   \itemize{
#'     \item \code{nothing} no conditioning, species can all go extinct
#'     \item \code{no_extinction} species cannot all go extinct
#'   }
mbd_sim_checked <- function(
  mbd_params, 
  crown_age = NA,
  stem_age = NA,
  conditioned_on = "nothing"
) {
  if (!is_mbd_params(mbd_params)) {
    stop("'mbd_params' must be a valid MBD parameter set")
  }
  if (!is.na(crown_age) && crown_age < 0.0) {
    stop("'crown_age' must be positive")
  }
  if (!is.na(stem_age) && stem_age < 0.0) {
    stop("'stem_age' must be positive")
  }
  if (is.na(crown_age) && is.na(stem_age)) {
    stop("'crown_age' or 'stem_age' must be set")
  }
  if (!is.na(crown_age) && !is.na(stem_age)) {
    stop("'crown_age' or 'stem_age' must be set exclusively")
  }
  if (!conditioned_on %in% c("nothing", "non_extinction")) {
    stop("'conditioned_on' must be either 'nothing' or 'non_extinction'")
  }
  # Data transformation
  pars <- as.numeric(unlist(mbd_params)) 
  soc <- 1
  age <- stem_age
  if (!is.na(crown_age)) {
    soc <- 2
    age <- crown_age 
  }
  cond <- 0
  if (conditioned_on == "non_extinction") {
    cond <- 1
  }
  mbd_sim(
    pars = pars, 
    soc = soc, 
    age = age, 
    cond = cond
  )
  
}