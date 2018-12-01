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
#' @inheritParams default_params_doc
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
#' @inheritParams default_params_doc
#' @return Package name
#' @export
get_pkg_name <- function() {
  pkg_name <- "mbd"
  pkg_name
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
  models
) {
  pkg_name <- get_pkg_name() # nolint internal function
  fun_list <- ls(paste0("package:", pkg_name))
  error_message <- paste0(
    "This is not a likelihood function provided by ",
    pkg_name,
    "!"
  )

  if (is.vector(models)) {
    function_names <- model_names <- which_function <- rep(NA, length(models))
    for (m in seq_along(models)) {
      fun <- eval(models[m])[[1]]
      if (is.character(models[m])) {
        if (length(
          (find_function <- which(fun_list == models[m]))
        ) == 0) {
          stop(error_message)
        }
        which_function[m] <- find_function
      } else {
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
      model_names[m] <- cut_loglik_from_name(function_names[m]) # nolint internal function
    }
  } else {
    fun <- eval(models)
    if (is.character(models)) {
      if (length(
        (find_function <- which(fun_list == models))
      ) == 0) {
        stop(error_message)
      }
      which_function <- find_function
    } else {
      for (i in seq_along(fun_list)) {
        if (all.equal(get(fun_list[i]), fun) == TRUE) {
          which_function <- i
        }
      }
    }
    if (is.null(which_function) | is.na(which_function)) {
      stop(error_message)
    }
    function_names <- toString(fun_list[which_function])
    model_names <- cut_loglik_from_name(function_names) # nolint internal function
  }

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
    model_names[m] <- cut_loglik_from_name(function_names[m]) # nolint internal function
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
