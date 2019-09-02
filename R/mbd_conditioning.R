#' Called by \link{mbd_loglik} if there is a conditioning != 0
#' @inheritParams default_params_doc
#' @return the conditional probability
#' @author Giovanni Laudanno
#' @export
calculate_conditional_prob <- function(
  brts,
  pars,
  cond,
  lx = 1000,
  n_0 = 2,
  tips_interval = c(n_0 * (cond > 0), Inf),
  methode = "lsodes",
  abstol = 1e-16,
  reltol = 1e-10,
  debug_mode = FALSE
) {
  check_cond(cond = cond, tips_interval = tips_interval, n_0 = n_0)
  if (cond == 0) {
    return(1)
  }
  total_time <- max(abs(brts))
  m <- 0:lx
  one_over_cm <- (3 * (m + 1)) / (m + 3)
  one_over_qm_binom <- 1 / choose(m + n_0, n_0)
  q_i <- c(1, rep(0, lx)); names(q_i) <- paste0("Q", 0:lx)
  testit::assert(length(one_over_cm) == length(m))
  testit::assert(length(one_over_qm_binom) == length(m))
  testit::assert(length(q_i) == length(m))
  # creating a_matrix
  matrix_a <- create_a(pars = pars, k = n_0, lx = lx) # nolint internal function
  # integrating the starting q_vector to t_p
  q_t <- mbd_solve(
    q_vector = q_i,
    matrix_a = matrix_a,
    time_interval = total_time,
    debug_mode = debug_mode
  )
  names(q_t) <- names(q_i)

  total_product <- q_t * one_over_cm * one_over_qm_binom
  missingspecies_min <- max(tips_interval[1] - n_0, 0)
  missingspecies_max <- min(tips_interval[2] - n_0, lx)
  # +1 is because of the zero-th component
  tips_components <- 1 + c(missingspecies_min, missingspecies_max)
  pc <- sum(total_product[tips_components[1]:tips_components[2]])

  if (!(pc >= 0 && pc <= 1)) {
    if (debug_mode == TRUE) {
      plot(
        q_t,
        xlab = "m",
        ylab = "Q_m^k(t_p - t_c)",
        main = paste0(
          "pars = ", pars[1], ", ", pars[2], ", ", pars[3], ", ", pars[4]
        ),
        sub = paste0(
          "Conditional probability = ", pc
        )
      )
      cat("The value of pc is: ", pc, "\n")
    } else {
      stop("problems: pc is wrong!")
    }
  } # debug
  pc
}
