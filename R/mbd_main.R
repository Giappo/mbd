#' @title Run the full mbd routine
#' @description Simulate an mbd process and infer parameters maximizing the
#' likelihood
#' @inheritParams default_params_doc
#' @details mle inference
#' @export
mbd_main <- function(
  seed,
  sim_pars,
  cond = 1,
  n_0 = 2,
  age = 10,
  tips_interval = c(0, Inf),
  start_pars = c(0.3, 0.1, 1, 0.1),
  loglik_function = mbd_loglik,
  verbose = TRUE
) {
  # set up
  pkg_name <- mbd_pkg_name()
  fun_list <- ls(paste0("package:", pkg_name))
  fun <- eval(loglik_function)
  if (is.character(loglik_function)) {
    which_function <- which(fun_list == loglik_function)
  } else {
    for (i in seq_along(fun_list)) {
      if (all.equal(get(fun_list[i]), fun) == TRUE) {
        which_function <- i
      }
    }
  }
  if (is.null(which_function)) {
    stop("This is not a likelihood function provided by mbd!")
  }
  fun_name <- toString(fun_list[which_function])
  if (verbose == TRUE) {
   cat("You are using the function:", fun_name)
  }

  # simulate
  set.seed(seed)
  sim <- mbd_sim(
    pars = sim_pars,
    n_0 = n_0,
    age = age,
    tips_interval = tips_interval,
    cond = cond
  )
  tips <- (n_0 - 1) + length(sim$brts)

  # maximum likelihood
  results <- data.frame(matrix(
    NA,
    nrow = 1,
    # ncol must be length pars + (loglik, df, conv) + (tips, seed)
    ncol = length(start_pars) + 3 + 2
  ))

  if (verbose == FALSE) {
    if (rappdirs::app_dir()$os != "win") {
      sink(file.path(rappdirs::user_cache_dir(), "ddd"))
    } else {
      sink(rappdirs::user_cache_dir())
    }
  }
  mle <- mbd_ml(
    true_pars = sim_pars,
    brts = sim$brts,
    tips_interval = tips_interval
  )
  if (verbose == FALSE) {
    sink()
  }
  results <- data.frame(
    mle,
    tips = tips,
    seed = seed
  )

  # format output
  colnames(results) <- c(
    colnames(mle),
    "tips",
    "seed"
  )
  out <- cbind(
    matrix(
      sim_pars,
      nrow = 1,
      ncol = length(start_pars),
      byrow = TRUE
    ),
    results,
    model = fun_name
  )
  colnames(out) <- c(
    paste0("sim_", colnames(mle[1:length(start_pars)])),
    colnames(results),
    "model"
  )
  out <- data.frame(out)

  # save data
  if (.Platform$OS.type == "windows") {
    sim_path  <- system.file("extdata", package = pkg_name)
    if (!file.exists(sim_path)) {
      dir.create(sim_path, showWarnings = FALSE)
    }
  } else {
    sim_path  <- getwd()
  }
  data_path <- file.path(sim_path, "data")
  data_file_name <- file.path(
    data_path,
    paste0(pkg_name, "_sim_", seed, ".RData")
  )
  if (!file.exists(data_file_name)) {
    if (!file.exists(data_path)) {
      dir.create(data_path, showWarnings = FALSE)
    }
  }
  save(sim, file = data_file_name)
  results_file_name <- file.path(
    sim_path,
    paste0(pkg_name, "_mle_", seed, ".txt")
  )
  utils::write.csv(
    x = out,
    file = results_file_name
  )

  out
}
