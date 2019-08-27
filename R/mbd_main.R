#' @title Run the full mbd routine
#' @description Simulate an mbd process and infer parameters maximizing the
#' likelihood
#' @inheritParams default_params_doc
#' @details mle inference
#' @export
mbd_main <- function(
  sim_pars,
  cond = 1,
  n_0 = 2,
  age = 10,
  tips_interval = c(n_0 * (cond > 0), Inf),
  seed,
  start_pars = c(0.3, 0.1, 1.5, 0.15),
  optim_ids = rep(TRUE, length(start_pars)),
  loglik_functions = mbd_loglik,
  methode = "lsoda",
  project_folder = NULL,
  verbose = FALSE,
  q_threshold = 1e-4
) {
  # check formats
  sim_pars <- as.numeric(sim_pars)
  start_pars <- as.numeric(start_pars)
  cond <- as.numeric(cond)
  seed <- as.numeric(seed)

  # specific set up
  n_0s <- n_0
  t_0s <- age

  # generic set up
  function_names <- get_function_names( # nolint internal function
    loglik_functions = loglik_functions
  )
  model_names <- get_model_names( # nolint internal function
    function_names = function_names,
    verbose = verbose
  )

  # simulate
  sim <- mbd::mbd_sim(
    pars = sim_pars,
    n_0 = n_0,
    age = age,
    cond = cond,
    tips_interval = tips_interval,
    seed = seed
  )
  brts <- sim$brts
  print_info(brts = brts, n_0 = n_0, cond = cond, verbose = verbose) # nolint internal function
  if (!is.list(brts)) {
    tips <- (n_0s[1] - 1) + length(brts)
  } else {
    tips <- rep(NA, length(brts))
    for (i in seq_along(brts)) {
      tips[i] <- n_0s[i] - 1 + length(brts[[i]])
    }
  }

  # maximum likelihood
  mle <- data.frame(matrix(
    NA,
    nrow = length(loglik_functions),
    # ncol must be length pars + (loglik, df, conv)
    ncol = length(start_pars) + 3
  ))
  for (m in seq_along(loglik_functions)) {
    if (verbose == FALSE) {
      sink(tempfile())
    }
    mle_out <- mbd_ml(
      loglik_function = get(function_names[m]),
      methode = methode,
      brts = sim$brts,
      cond = cond,
      n_0 = n_0,
      tips_interval = tips_interval,
      start_pars = start_pars,
      optim_ids = optim_ids,
      true_pars = sim_pars,
      verbose = verbose,
      lx = min(1 + 3 * (length(brts)), 1200),
      q_threshold = q_threshold
    )
    mle[m, ] <- mle_out
    colnames(mle) <- names(mle_out)
    if (verbose == FALSE) {
      sink()
    }
  }

  # format results
  results <- cbind(
    matrix(
      sim_pars,
      nrow = length(loglik_functions),
      ncol = length(start_pars),
      byrow = TRUE
    ),
    mle,
    seed,
    cond,
    n_0,
    matrix(
      t_0s,
      nrow = length(loglik_functions),
      ncol = length(t_0s),
      byrow = TRUE
    ),
    matrix(
      tips,
      nrow = length(loglik_functions),
      ncol = length(tips),
      byrow = TRUE
    ),
    matrix(
      optim_ids,
      nrow = length(loglik_functions),
      ncol = length(optim_ids),
      byrow = TRUE
    ),
    model = model_names
  )
  if (length(t_0s) > 1) {
    t_0s_label <- paste0("t_0_", 1:length(t_0s))
  } else {
    t_0s_label <- "t_0"
  }
  if (length(tips) > 1) {
    tips_label <- paste0("tips_", 1:length(tips))
  } else {
    tips_label <- "tips"
  }
  colnames(results) <- c(
    paste0("sim_", colnames(mle[1:length(start_pars)])),
    colnames(mle),
    "seed",
    "cond",
    "n_0",
    t_0s_label,
    tips_label,
    paste0("optim_", colnames(mle[1:length(start_pars)])),
    "model"
  )
  rownames(results) <- NULL
  results <- data.frame(results)

  # save data
  main_save_files( # nolint internal function
    project_folder = project_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed,
    sim = sim,
    results = results
  )
  print_info(brts = brts, n_0 = n_0, cond = cond, verbose = verbose) # nolint internal function
  results
}
