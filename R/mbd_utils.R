#---- package specific functions
#' Get the names of the parameters used in the mbd model
#' @author Giovanni Laudanno
#' @export
get_param_names <- function() {
  c("lambda", "mu", "nu", "q")
}

#' @title Loglik functions in mbd package
#' @author Giovanni Laudanno
#' @description Get all the loglik functions in mbd package
#' @return loglik functions
#' @export
mbd_logliks_function <- function() {
  fun_list <- ls(paste0("package:", get_pkg_name())) # nolint internal function
  brts_is_arg <- rep(0, length(fun_list))
  for (i in seq_along(fun_list)) {
    brts_is_arg[i] <- sum(grepl("brts", deparse(args(fun_list[i]))))
  }
  ll_funs <- fun_list[sapply(
    fun_list, function(x)
      any(grepl("loglik", x))
  ) & brts_is_arg]
  ll_funs
}

#' @title Get package name
#' @author Giovanni Laudanno
#' @description Get package name
#' @return Package name
#' @export
get_pkg_name <- function() {
  pkg_name <- "mbd"
  pkg_name
}

#' @title Conditionings
#' @author Giovanni Laudanno
#' @description Gives the conditionings accepted by mbd
#' @return the conditionings
#' @export
mbd_conds <- function() {
  conds <- c(0, 1)
  conds
}

#' @title Eqs for condprob
#' @description Eqs for condprob
#' @author Giovanni Laudanno
#' @return the conditionings
#' @export
mbd_condprob_eqs <- function() {
  eqs <- c("p_eq", "q_eq", "sim", "nee")
  eqs
}

#' @title Logliks for the experiment
#' @author Giovanni Laudanno
#' @description Get the loglik functions to use for the experiment
#' @return loglik functions to use for the experiment
#' @export
mbd_experiment_logliks <- function() {
  fun_list <- ls(paste0("package:", get_pkg_name())) # nolint internal function
  mbd_funs <- fun_list[sapply(
    fun_list, function(x)
      any(grepl("mbd_loglik", x)) &
      !any(grepl("loglik_", x)) &
      !any(grepl("logliks", x))
  )]
  mbd_funs
}

#' Generates parameters for the mbd experiment
#' @inheritParams default_params_doc
#' @return the parameter table for the experiment
#' @author Giovanni Laudanno
#' @export
mbd_experiment_pars <- function(
  lambdas = c(0.3),
  mus = c(0, 0.1),
  nus = c(0, 0.5, 1.0, 1.5),
  qs = c(0.1, 0.15, 0.2),
  cond = 1,
  age = 8
) {

  x <- expand.grid(lambda = lambdas, mu = mus, nu = nus, q = qs)
  no_mbd_lines <-
    apply(X = x, MARGIN = 1, FUN = function(y) y[3] == 0 | y[4] == 0)
  no_mbd_x <- expand.grid(lambda = lambdas, mu = mus, nu = 0, q = 0)
  x[no_mbd_lines, ] <- no_mbd_x
  x <- dplyr::distinct(x)
  x$age <- age
  x$cond <- cond
  rownames(x) <- NULL
  x
}

#---- general functions
#' @title cat2
#' @author Giovanni Laudanno
#' @description If verbose == TRUE cats the message, otherwise stays silent
#' @inheritParams default_params_doc
#' @return prints on screen
#' @export
cat2 <- function(
  message,
  verbose
) {
  if (verbose == TRUE) {
    cat(message)
  } else {
    return()
  }
}

#' Convert to string
#' @param var a variable
#' @export
to_string2 <- function(
  var
) {
  gsub(x = toString(var), pattern = " ", replacement = "")
}

#' @title Transform parameters
#' @description Transform parameters according to y = x / (1 + x)
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @return transformed parameters
#' @export
pars_transform_forward <- function(pars) {
  pars <- as.numeric(unlist(pars))
  pars_transformed <- pars / (1 + pars)
  pars_transformed[which(pars == Inf)] <- 1
  pars_transformed
}

#' @title Transform parameters back
#' @description Transform parameters back according to x = y / (1 + y)
#' @inheritParams default_params_doc
#' @details This is not to be called by the user.
#' @return the original parameters
#' @export
pars_transform_back <- function(pars_transformed) {
  pars_transformed <- as.numeric(unlist(pars_transformed))
  pars <- pars_transformed / (1 - pars_transformed)
  pars
}

#' @title Cut word "loglik" from a name
#' @author Giovanni Laudanno
#' @description Cut word "loglik" from a name
#' @inheritParams default_params_doc
#' @return clean name
#' @export
cut_loglik_from_name <- function(
  function_name
) {
  if (grepl("loglik", function_name)) {
    model_name <- gsub(
      "_loglik",
      "",
      gsub(
        "loglik_",
        "",
        function_name
      )
    )
  } else {
    model_name <- NA
  }
  model_name
}

#' @title Get function names
#' @author Giovanni Laudanno
#' @description Get function names
#' @inheritParams default_params_doc
#' @return function names
#' @export
get_function_names <- function(
  loglik_functions
) {
  pkg_name <- mbd::get_pkg_name() # nolint internal function
  fun_list <- ls(paste0("package:", pkg_name))
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    pkg_name,
    "!"
  )

  if (!is.list(loglik_functions)) {
    if (length(loglik_functions) == 1) {
      loglik_functions <- list(loglik_functions)
    }
  }

  function_names <- rep(NA, length(loglik_functions))
  model_names <- which_function <- function_names
  for (m in seq_along(loglik_functions)) {
    if (is.character(loglik_functions[[m]])) {
      if (length(
        (find_function <- which(fun_list == loglik_functions[[m]]))
      ) == 0) {
        stop(error_message)
      }
      which_function[m] <- find_function
    } else {
      fun <- eval(loglik_functions[m])[[1]]
      for (i in seq_along(fun_list)) {
        if (all.equal(get(fun_list[i]), fun) == TRUE) {
          which_function[m] <- i
        }
      }
    }
    if (is.null(which_function[m]) | is.na(which_function[m])) {
      stop(error_message)
    }
    function_names[m] <- toString(fun_list[which_function[m]])
    model_names[m] <- mbd::cut_loglik_from_name(function_names[m])
  }
  function_names <- fun_list[which_function]
  model_names <- unname(
    sapply(function_names, FUN = mbd::cut_loglik_from_name) # nolint internal function
  )
  if (any(is.na(model_names))) {
    stop(error_message)
  }
  invisible(function_names)
}

#' @title Check if provided models make sense
#' @author Giovanni Laudanno
#' @description Check if provided models make sense
#' @inheritParams default_params_doc
#' @return models names
#' @export
get_model_names <- function(
  function_names,
  verbose = FALSE
) {
  model_names <- function_names
  for (m in seq_along(function_names)) {
    model_names[m] <- mbd::cut_loglik_from_name(function_names[m])
    if (is.null(model_names[m]) | is.na(model_names[m])) {
      stop(paste0(
        "This is not a likelihood function provided by ",
        get_pkg_name(), # nolint internal function
        "!"
      ))
    }
  }
  if (verbose == TRUE) {
    cat("You are using the functions:", model_names)
  }
  model_names
}

#' @title Read saved results
#' @author Giovanni Laudanno
#' @description Read saved results
#' @inheritParams default_params_doc
#' @return results
#' @export
read_results <- function(project_folder = NULL) {
  if (is.null(project_folder)) {
    if (.Platform$OS.type == "windows") {
      project_folder <- system.file("extdata", package = get_pkg_name())
    }
  }
  if (!dir.exists(project_folder)) {
    stop("This directory does not exist")
  }
  if (length(list.files(project_folder)) == 0) {
    stop(paste0(project_folder, " is empty."))
  }
  dir_results <- file.path(project_folder, "results")
  dir_data <- file.path(project_folder, "data")

  files_results <- list.files(dir_results)
  if (length(files_results) == 0) {
    stop(paste0(dir_results, " is empty."))
  }
  files_data <- list.files(dir_data)
  print("Saved results are:")
  print(files_results)
  x <- readline("Which one do you want to read?")
  result <- utils::read.csv(
    file.path(
      dir_results,
      files_results[as.numeric(x)]
    )
  )
  data <- get(load(
    file.path(
      dir_data,
      files_data[as.numeric(x)]
    )
  ))
  return(
    list(
      data = data,
      result = result
    )
  )
}

#' Count the number of multiple speciation events
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
count_n_mb_events <- function(brts) {
  births <- mbd::brts2time_intervals_and_births(brts)$births # nolint internal function
  mb_species <- sum(births[births > 1])
  mb_species
}

#' Count the percentage of multiple speciation events
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
count_percentage_mb_events <- function(brts) {
  mb_species <- mbd::count_n_mb_events(brts)
  if (mb_species == 0) {
    percentage_mb_species <- 0
  } else {
    births <- mbd::brts2time_intervals_and_births(brts)$births # nolint internal function
    n_species <- sum(births)
    percentage_mb_species <- mb_species / n_species
  }
  percentage_mb_species
}

#' Create an empty phylogeny
#' @author Giovanni Laudanno
#' @export
create_empty_phylo <- function() {
  tr <- list(edge = matrix(c(2, 1), 1, 2), tip.label = "", Nnode = 0L)
  class(tr) <- "phylo"
  tr$tip.label <- c() # nolint
  tr
}

#' Create an empty phylogeny
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
create_singleton_phylo <- function(age) {
  tr <- list(edge = matrix(c(2, 1), 1, 2), tip.label = "t1", Nnode = 1L)
  class(tr) <- "phylo"
  tr$edge.length <- age # nolint
  tr
}

#' Like file.path, but cooler
#' @param fsep path separator for the OS
#' @param ... additional arguments
#' @export
file_path <- function(..., fsep = .Platform$file.sep) {
  gsub("//", "/", file.path(..., fsep = fsep))
}

#' Create an empty phylogeny
#' @inheritParams default_params_doc
#' @author Giovanni Laudanno
#' @export
get_pars_filename <- function(
  pars,
  age
) {
  # folder structure
  data_folder0 <- system.file("extdata", package = get_pkg_name())
  if (!dir.exists(data_folder0)) {
    dir.create(data_folder0)
  }
  data_folder <- file.path(data_folder0, "cond_probs")
  if (!dir.exists(data_folder)) {
    dir.create(data_folder)
  }

  parsetting <- paste0(
    "pars=", mbd::to_string2(signif(pars, digits = 3)),
    "-",
    "age=", mbd::to_string2(signif(age, digits = 3))
  )
  parsetting <- gsub(parsetting, pattern = " ", replacement = "")
  pars_filename <- file.path(data_folder, parsetting)

  pars_filename
}

#' Estimate approximately the expected number of taxa
#' @inheritParams default_params_doc
#' @return a list with mean and standard deviation
#' @export
estimate_n_taxa <- function(pars, n_0, age) {
  nee_pars <- mbd::get_nee_pars(pars = pars)
  nee_mean <- mbd::nee_mean_nt(nee_pars = nee_pars, n_0 = n_0, age = age)
  nee_stdev <- mbd::nee_stdev_nt(nee_pars = nee_pars, n_0 = n_0, age = age)
  list(
    mean = nee_mean,
    stdev = nee_stdev
  )
}

#' My version of tryCatch
#' @inheritParams default_params_doc
#' @export
my_try_catch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(
      expr, error = function(e) {
        err <<- e
        NULL
      }
    ), warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(
    value = value,
    warning = warn,
    error = err
  )
}
