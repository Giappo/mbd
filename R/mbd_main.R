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
  tips_interval = c(0, Inf),
  seed,
  start_pars = c(0.3, 0.1, 1, 0.1),
  optim_ids = rep(TRUE, length(start_pars)),
  models = mbd_loglik,
  project_folder = NULL,
  verbose = FALSE
) {
  # specific set up
  n_0s <- n_0
  t_0s <- age

  # generic set up
  function_names <- get_function_names( # nolint internal function
    models = models
  )
  model_names <- get_model_names( # nolint internal function
    function_names = function_names,
    verbose = verbose
  )

  # simulate
  set.seed(seed)
  sim <- mbd_sim(
    pars = sim_pars,
    n_0 = n_0,
    age = age,
    cond = cond,
    tips_interval = tips_interval
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
  for (m in seq_along(models)) {
    if (verbose == FALSE) {
      if (rappdirs::app_dir()$os != "win") {
        sink(file.path(rappdirs::user_cache_dir(), "ddd"))
      } else {
        sink(rappdirs::user_cache_dir())
      }
    }
    mle <- mbd_ml(
      brts = sim$brts,
      n_0 = n_0,
      cond = cond,
      tips_interval = tips_interval,
      verbose = verbose
    )
    if (verbose == FALSE) {
      sink()
    }
  }

  # format results
  results <- cbind(
    matrix(
      sim_pars,
      nrow = length(models),
      ncol = length(start_pars),
      byrow = TRUE
    ),
    mle,
    seed,
    cond,
    n_0,
    matrix(
      t_0s,
      nrow = length(models),
      ncol = length(t_0s),
      byrow = TRUE
    ),
    matrix(
      tips,
      nrow = length(models),
      ncol = length(tips),
      byrow = TRUE
    ),
    matrix(
      optim_ids,
      nrow = length(models),
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
