#' @title Print information about simulated data
#' @author Giovanni Laudanno
#' @description Print information about simulated data
#' @inheritParams default_params_doc
#' @return nothing
#' @export
print_info <- function(
  brts,
  n_0,
  cond,
  verbose
) {
  n_0s <- c(n_0, rep(1, is.list(brts) * (length(brts) - 1)))
  if (verbose == FALSE) {
    return()
  }
  cat("\ncond =", cond, "\n")
  if (is.list(brts)) {
    cat("\n")
    for (i in seq_along(brts)) {
      cat(paste0("clade ", i, ":\n"))
      cat(paste0("n_0 = ", n_0s[[i]], "\n"))
      cat(paste(c(
        "branching times = c(",
        paste(signif(brts[[i]], digits = 3), collapse = ", "),
        ")\n"), collapse = ""
      ))
      cat("\n")
    }
  } else {
    cat(paste0("n_0 = ", n_0s, "\n"))
    cat(paste(c(
      "branching times = c(",
      paste(signif(brts, digits = 3), collapse = ", "),
      ")\n"), collapse = ""
    ))
    cat("\n")
  }
}

#' @title Save data and results
#' @description Save data and results
#' @inheritParams default_params_doc
#' @export
main_save_files <- function(
  project_folder,
  sim_pars,
  optim_ids,
  cond,
  n_0,
  t_0s,
  seed,
  sim,
  results
) {
  # project folder
  pkg_name <- get_pkg_name() # nolint internal function
  if (is.null(project_folder)) {
    if (.Platform$OS.type == "windows") {
      if (!("extdata" %in% list.files(system.file(package = pkg_name)))) {
        dir.create(file.path(
          system.file(package = pkg_name),
          "extdata"
        ))
      }
      project_folder <- system.file("extdata", package = pkg_name)
      if (!file.exists(project_folder)) {
        dir.create(project_folder, showWarnings = FALSE)
      }
    } else {
      project_folder <- getwd()
    }
  }

  # data
  data_folder <- file.path(project_folder, "data")
  if (!dir.exists(data_folder)) {
    dir.create(data_folder, showWarnings = FALSE)
  }
  data_file_name <- mbd::create_data_file_name(
    data_folder = data_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  save(sim, file = data_file_name)

  # results
  results_folder <- file.path(
    project_folder,
    "results"
  )
  if (!dir.exists(results_folder)) {
    dir.create(results_folder, showWarnings = FALSE)
  }
  results_file_name <- mbd::create_results_file_name(
    results_folder = results_folder,
    sim_pars = sim_pars,
    optim_ids = optim_ids,
    cond = cond,
    n_0 = n_0,
    t_0s = t_0s,
    seed = seed
  )
  if (file.exists(results_file_name)) {
    previous_results <- utils::read.csv(
      file = results_file_name,
      header = TRUE
    )[, -1]
    all_results <- rbind(
      previous_results,
      results
    )
    all_results <- all_results[order(all_results[, ncol(all_results)]), ]
    all_results <- all_results[!duplicated(all_results), ]
    file.remove(results_file_name)
    utils::write.csv(
      x = all_results,
      file = results_file_name
    )
    res2 <- utils::read.csv(
      file = results_file_name,
      header = TRUE
    )[, -1]
    res2 <- res2[order(res2[, ncol(res2)]), ]
    res2 <- res2[!duplicated(res2), ]
    file.remove(results_file_name)
    utils::write.csv(
      x = res2,
      file = results_file_name
    )
  } else {
    utils::write.csv(
      x = results,
      file = results_file_name
    )
  }
}

#' @title Create results file name
#' @description Create results file name
#' @inheritParams default_params_doc
#' @export
create_results_file_name <- function(
  results_folder,
  sim_pars,
  optim_ids,
  cond,
  n_0,
  t_0s,
  seed
) {
  results_file_name <- file.path(
    results_folder,
    paste0(
      get_pkg_name(),
      "_mle",
      "-",
      paste0(
        "sim_pars=[",
        utils::capture.output(cat(
          paste(sim_pars, collapse = "-")
        )),
        "]"
      ),
      "-",
      paste0(
        "optim_ids=[",
        utils::capture.output(cat(
          paste(as.numeric(optim_ids), collapse = "-")
        )),
        "]"
      ),
      "-cond=",
      cond,
      "-n_0=",
      n_0,
      "-",
      paste0(
        "ages=[",
        utils::capture.output(cat(
          paste(t_0s, collapse = "-")
        )),
        "]"
      ),
      "-seed=",
      seed,
      ".txt"
    )
  )
  results_file_name
}

#' @title Create data file name
#' @description Create data file name
#' @inheritParams default_params_doc
#' @export
create_data_file_name <- function(
  data_folder,
  sim_pars,
  optim_ids,
  cond,
  n_0,
  t_0s,
  seed
) {
  data_file_name <- file.path(
    data_folder,
    paste0(
      get_pkg_name(),
      "_sim",
      "-",
      paste0(
        "sim_pars=[",
        utils::capture.output(cat(
          paste(sim_pars, collapse = "-")
        )),
        "]"
      ),
      "-",
      paste0(
        "optim_ids=[",
        utils::capture.output(cat(
          paste(as.numeric(optim_ids), collapse = "-")
        )),
        "]"
      ),
      "-cond=",
      cond,
      "-n_0=",
      n_0,
      "-",
      paste0(
        "ages=[",
        utils::capture.output(cat(
          paste(t_0s, collapse = "-")
        )),
        "]"
      ),
      "-seed=",
      seed,
      ".RData"
    )
  )
  data_file_name
}
