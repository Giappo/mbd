#' Called by \link{mbd_loglik} if and only if conditioned on non-extinction
#' @inheritParams default_params_doc
#' @return the pc. If \code{is.nan(log(pc))} the log-likelihood estimation
#'   by \link{mbd_loglik} failed
#' @author Giovanni Laudanno
#' @noRd
calculate_conditional_probability <- function(
  brts,
  pars,
  lx = 1000,
  n_0 = 2,
  tips_interval = c(0, Inf),
  methode = "expo",
  abstol = 1e-16,
  reltol = 1e-10
) 
{
  lambda <- pars[1]; mu <- pars[2]; nu <- pars[3]; q <- pars[4];
  total_time <- max(abs(brts));
  
  m <- 0:lx; length(m)
  one_over_cm <- (3 * (m + 1)) / (m + 3); length(one_over_cm)
  one_over_qm_binom <- 1 / choose((m + n_0), n_0)
  q_i <- c(1, rep(0, lx)); length(q_i)
  
  matrix_a <- create_a(
    pars = pars, 
    k = n_0,
    lx = lx
  )
  
  a2_v1 <- mbd:::a_operator(
    q_vector = q_i,
    transition_matrix = matrix_a,
    time_interval = total_time,
    precision = 250L,
    methode = methode,
    abstol = abstol,
    reltol = reltol
  )
  
  total_product <- a2_v1 * one_over_cm * one_over_qm_binom
  missingspecies_min <- max((tips_interval[1] - n_0), 0)
  missingspecies_max <- min((tips_interval[2] - n_0), lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  sum(total_product[tips_components[1]:tips_components[2]])
}
