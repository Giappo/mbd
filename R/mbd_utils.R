#' analysis function
#' @inheritParams default_params_doc
#' @param path something
#' @param titolo something
#' @param pdfname something
#' @param percentage_hidden_outliers something
#' @param openit something
#' @param mother_folder something
correlation_analysis <- function(
  results,
  path,
  titolo = NULL,
  pdfname,
  sim_pars = sim_pars,
  percentage_hidden_outliers = 0.04,
  openit = 0,
  idparsopt,
  mother_folder
) {
  n_pars <- length(sim_pars)
  if (missing(idparsopt)) {
    estimated_pars <- 1:n_pars
  } else {
    estimated_pars <- idparsopt
  }
  par_names <- colnames(results)[1:n_pars]

  truevalues_color <- "red"
  points_color <- "azure4"
  medians_color <- "blue3" #"chartreuse3";
  medians_color_name <- "blue"
  truevalues_color_name <- "red"

  medians <- rep(0, n_pars)
  for (idpar in 1:n_pars) {
    medians[idpar] <- stats::median(results[, idpar])
  }
  medians_string <- paste0(
    "MLE Medians (", medians_color_name, ") = (", paste(signif(medians, 2),
    sep = "''", collapse = ", "), ")"
  )
  truevalues_string <- paste0(
    "True Values (", truevalues_color_name, ") = (", paste(signif(sim_pars, 2),
      sep = "''", collapse = ", "), ")"
  )
  axislimits <- rep(NA, n_pars)
  # BUG? 'probs = 1 - percentage_hidden_outliers' used to be with arrow operator
  for (i in 1:n_pars) {
    axislimits[i] <- stats::quantile(
      results[, i],
      probs = 1 - percentage_hidden_outliers
    )
  }

  #pdf
  grDevices::pdf(file = paste(path, "/", pdfname, ".pdf", sep = ""))
  graphics::par(mfrow = c(length(estimated_pars), length(estimated_pars)))
  graphics::par(oma = c(0, 0, 2, 0));
  for (i in estimated_pars) {
    for (j in estimated_pars) {
      good_lines <- results[, i] < axislimits[i] & results[, j] < axislimits[j]
      if (any(good_lines) > 0) {
        good_results <- results[good_lines, ]
      } else {
        good_results <- results
      }

      if (i == j) {
        graphics::hist((good_results[, i]), main = NULL,
          xlab = paste(par_names[i]), breaks = 15
        )
        graphics::abline(v = sim_pars[i], col = truevalues_color)
        graphics::abline(v = medians[i], col = medians_color)
      }
      else {
        graphics::plot(
          good_results[, i] ~ good_results[, j],
          xlab = par_names[j],
          ylab = par_names[i],
          cex = 0.3,
          col = points_color
        )
        graphics::points(
          x = sim_pars[j], y = sim_pars[i], col = truevalues_color, pch = 10,
          cex = 1.5
        )
        graphics::points(
          x = medians[j], y = medians[i], col = medians_color, pch = 10,
          cex = 1.5
        )
      }
    }
  }
  titolo_pdf <- paste0(
    "\n\n", titolo, "\n", medians_string, "\n", truevalues_string
  )
  graphics::title(main = titolo_pdf, outer = TRUE)
  grDevices::dev.off()

  if (openit == 1) {
    file.show(normalizePath(paste(path, "/", pdfname, ".pdf", sep = "")))
  }

  if (!missing(mother_folder)) {
    grDevices::pdf(file = paste(mother_folder, "/", pdfname, ".pdf", sep = ""))
    graphics::par(mfrow = c(length(estimated_pars), length(estimated_pars)))
    graphics::par(oma = c(0, 0, 2, 0));
    for (i in estimated_pars) {
      for (j in estimated_pars) {
        good_lines <- results[, i] < axislimits[i] &
          results[, j] < axislimits[j]
        if (any(good_lines) > 0) {
          good_results <- results[good_lines, ]
        } else {
          good_results <- results
        }

        if (i == j) {
          graphics::hist(
            good_results[, i],
            main = NULL, xlab = paste(par_names[i]), breaks = 15
          )
          graphics::abline(v = sim_pars[i], col = truevalues_color)
          graphics::abline(v = medians[i], col = medians_color)
        } else {
          graphics::plot(
            good_results[, i] ~ good_results[, j],
            xlab = par_names[j],
            ylab = par_names[i],
            cex = 0.3,
            col = points_color
          )
          graphics::points(
            x = sim_pars[j],
            y = sim_pars[i],
            col = truevalues_color,
            pch = 10,
            cex = 1.5
          )
          graphics::points(
            x = medians[j],
            y = medians[i],
            col = medians_color,
            pch = 10,
            cex = 1.5
          )
        }
      }
    }
    graphics::title(main = titolo_pdf, outer = TRUE)
    grDevices::dev.off()
  }
}

#' Percentiles function
#' @inheritParams default_params_doc
percentiles_function <- function(
  results,
  sim_pars,
  printit = 1,
  quantiles_choice = c(.25, .50, .75)
) {
  quantiles_names <- format(round(quantiles_choice, 2), nsmall = 2)
  n_pars <- length(sim_pars)
  parnames <-  colnames(results)[1:n_pars]
  percentiles <- vector("list", 3)
  for (idpar in 1:n_pars) {
    percentiles[[idpar]] <- stats::quantile(results[, idpar], quantiles_choice)
  }
  percentiles <- t(matrix(unlist(percentiles), nrow = 3, ncol = n_pars));
  colnames(percentiles) <- quantiles_names; rownames(percentiles) <- parnames
  if (printit == 1) {
    print(percentiles); print(sim_pars)
  }
  percentiles
}

#' Converts a matrix in a dataframe that can be used with ggplot
#' @inheritParams default_params_doc
#' @param matrix something
#' @param heatmap_name heatmap name
heatmap2dataframe <- function(
  x,
  y,
  matrix,
  x_name = "x",
  y_name = "y",
  heatmap_name = "HeatMap"
) {
  matrix_2 <- matrix(matrix, nrow = length(x) * length(y))
  df <- data.frame(expand.grid(x = x, y = y), HeatMap = matrix_2)
  colnames(df)[1] <- x_name
  colnames(df)[2] <- y_name
  colnames(df)[3] <- heatmap_name
  df
}

#' Converts branching times to "time intervals between branching times"
#'   and "birth at nodes" vectors
#' @inheritParams default_params_doc
#' @noRd
brts2time_intervals_and_births <- function(
  brts
) {
  time_points <- -unlist(unname(sort(abs(brts), decreasing = T)))
  branching_times <- -sort(abs(as.numeric(time_points)), decreasing = TRUE)
  births <- c(0, unname(table(branching_times))[-1])
  unique_branching_times <- as.numeric(names(table(branching_times)))
  time_intervals <- c(
    diff(unique_branching_times),
    abs(utils::tail(unique_branching_times, 1))
  )
  births <- births[-1]
  list(time_intervals = time_intervals, births = births)
}

# @Giappo: add doc
#' Gets heatmap of my matrices.
#' @inheritParams default_params_doc
#' @param ... something
#' @export
myheatmap <- function(
  matrix,
  logs=1,
  colormap = grDevices::heat.colors(15),
  ...
){
  if (is.matrix(matrix) == FALSE) {
    matrix <- as.matrix(matrix)
  }
  matrix_t <- t(matrix[nrow(matrix):1, 1:ncol(matrix)])
  if (logs == 1) {
    graphics::image(log(matrix_t), col = rev(colormap), ...)
  } else {
    graphics::image(matrix_t)
  }
}

# @Giappo: add doc
#' Does something J
#' @inheritParams default_params_doc
myheatmap2 <- function(x, y, z, x_name, y_name, z_name, x_splits, y_splits
) {
  if (missing(x_splits)) {
    x_splits <- round((length(x)) / 10)
  }
  if (missing(y_splits)) {
    y_splits <- round((length(y)) / 10)
  }
  if (missing(x_name)) {
    x_name <- NULL
  }
  if (missing(y_name)) {
    y_name <- NULL
  }
  if (missing(z_name)) {
    z_name <- NULL
  }

  l_x <- x
  l_y <- y
  pretty_x_at	<-	pretty(range(l_x), n = x_splits)
  pretty_x_lab	<-	round(pretty_x_at, 2)
  pretty_y_at	<-	pretty(range(l_y), n = y_splits)
  pretty_y_lab	<-	round(pretty_y_at, 2)

  jet_colors <- grDevices::colorRampPalette(
    c(
      "#00007F", "blue", "#007FFF", "cyan",
      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"
    )
  )
  graphics::filled.contour(
    t(z),
    color = jet_colors,
    nlevels = 100,
    ylab = y_name,
    xlab = x_name,
    plot.axes = {
      graphics::axis(
        1, at = seq(0, 1, length.out = length(pretty_x_lab)),
        labels = pretty_x_lab
      )
      graphics::axis(
        2, at = seq(0, 1, length.out = length(pretty_y_lab)),
        labels = pretty_y_lab
      )
    },
    key_title = graphics::title(main = z_name)
  )
}

#' Checks for NA, NaN or negative components in a vector (usually used for Qt)
#' @inheritParams default_params_doc
#' @param v sometjing
#' @param display_output something
#' @noRd
negatives_correction <- function(
  v,
  pars,
  display_output = 0
){
  problems <- 0
  if (any(is.na(v))) {
    problems <- 1
    na_components <- which(is.na(v) & !is.nan(v))
    nan_components <- which(is.nan(v))
    if (display_output == 1) {
      cat("There are non-numeric components for par values:", pars, "\n")
      if (length(na_components) > 0) {
        cat("NA component are:", na_components)
      }
      if (length(nan_components) > 0) {
        cat("NaN component are:", nan_components)
      }
    }
  }

  # BUG? Used to be 'any(v < 0) & (problems) == 0'
  if (any(v < 0) && problems == 0) {
    v[v < 0 & (abs(v) / abs(max(v))) < 1e-10] <- 0
  }
  v
}

#' Compares two functions that accept same arguments args=(...)
#' @inheritParams default_params_doc
#' @param function1 something
#' @param function2 something
#' @param output something
#' @param ... something
compare_functions <- function(
  function1,
  function2,
  iterations = 1,
  output = 1,
  ...
){
  a_time1 <- 0
  a_time2 <- 0
  for (i in 1:iterations) {
    t0 <- proc.time()
    a_test1 <- (function1(...))
    a_time1 <- a_time1 + proc.time() - t0
    t1 <- proc.time()
    a_test2 <- (function2(...))
    a_time2 <- a_time2 + proc.time() - t1
  }
  if (output == 1) {
    print(list(a_test1, a_test2))
  }
  print(a_time1)
  print(a_time2)
  print(all.equal(a_test1, a_test2))
  return(invisible(list(a_test1, a_test2)))
}

#' Checks function calls within a function
#' @inheritParams default_params_doc
called_functions <- function(
  function_name,
  checked_functions = NULL
) {

    # Get the function's code:
    function_code <- deparse(get(function_name))

    # break code up into sections preceding left brackets:
    left_brackets <- c(
      unlist(strsplit(function_code, split = "[[:space:]]*\\("))
    )

    called_functions <- unique(
      c(
        unlist(
          sapply(
            left_brackets, function(x) {
              # Split up according to anything that can't be in a function name.
              # split = not alphanumeric, not '_', and not '.'
              words <- c(unlist(strsplit(x, split = "[^[:alnum:]_.]")))

              last_word <- utils::tail(words, 1)
              last_word_is_function <- tryCatch(
                is.function(get(last_word)),
                error = function(e) return(FALSE)
              )
              return(last_word[last_word_is_function])
            }
          )
        )
      )
    )
    called_functions
}

#' Checks if some matrix entries are infinite, NaN, NA
#'   or (only in the lower triangle) negative
#' @inheritParams default_params_doc
#' @param m_list something
#' @param sign_check something
matrix_check <- function(
  m_list,
  sign_check = 0
) {
  max_k <- length(m_list)
  black_list <- as.complex(rep(0, max_k))
  negative_list <- rep(0, max_k)
  for (k in 1:max_k) {
    matrix <- as.matrix(m_list[[k]])
    if ((max(is.infinite(matrix)) > 0) |
      (max(is.nan(matrix)) > 0) |
      (max(is.na(matrix)) > 0)
    ) {
      black_list[k] <- black_list[k] + 1
    }
    if (max(diag(matrix)) > 0) {
      black_list[k] <- black_list[k] + 2i
    }
  }
  if (sum(Re(black_list)) != 0) {
    print(
      paste(
        "We have non-numeric values for k=",
        (1:length(black_list))[(Re(black_list) > 0)]
      )
    )
  }
  if (sum(Im(black_list)) != 0) {
    print(
      paste(
        "We have positive diagonal values for k=",
        (1:length(black_list))[(Im(black_list) > 0)]
      )
    )
  }
  if (sign_check == 1) {
    entries_sign <- rep(1, max_k)
    for (k in 1:max_k) {
      lower_triangle <- (as.matrix(m_list[[k]]))[lower.tri(m_list[[k]])]
      if (min(lower_triangle) < 0) {
        entries_sign[k] <- -1
        negative_list[k] <- -1
      }
    }
    if (min(entries_sign) < 0) {
      print(
        paste(
          "There is a negative entry on k=",
          which((entries_sign) < 0)
        )
      )
    }
  }
  list(black_list = black_list, negative_list = negative_list)
}

#' Does something
#' @inheritParams default_params_doc
append_multiple <- function(
  x,
  values,
  after
) {
  if (!is.list(values) & length(values) != length(after)) {
    stop(
      "In append_multiple: values must be a list ",
      "unless it has the same length as after"
    )
  }

  x2 <- x
  for (i in 1:length(values)) {
    x2 <- append(x2, values[[i]], after[i] + length(x2) - length(x))
  }
  x2
}

#' Function returning branch lengths named by ending node
#' @inheritParams default_params_doc
get_branch_lengths <- function(tr) {
  if (length(tr$node) > 0) {
    type <- "NAMED"
  } else {
    type <- "UNNAMED"
  }
  if (length(tr$root) > 0) {
    root_length <- tr$root
  } else{
    root_length <- 0
  }

  node_numbers <- sort(unique(tr$edge[, 1]))
  b_l_1 <- as.numeric(
    stats::na.omit(
      c(
        root_length,
        tr$edge.length[match(node_numbers, tr$edge[, 2])]
      )
    )
  )

  tip_numbers <- 1:length(tr$tip)
  b_l_2 <- tr$edge.length[match(tip_numbers, tr$edge[, 2])]

  if (type == "NAMED") {
    res <- data.frame(
      Number = c(sort(node_numbers), tip_numbers),
      Name = c(tr$node, tr$tip),
      Lengths = c(b_l_1, b_l_2)
    )
  }
  if (type == "UNNAMED") {
    res <- data.frame(
      Number = c(sort(node_numbers), tip_numbers),
      Name = c(rep(NA, tr$Nnode), tr$tip),
      Lengths = c(b_l_1, b_l_2)
    )
  }

  return(res)
}

#' Does something
#' @inheritParams default_params_doc
summarize_beast_posterior <- function(
  input_trees_path,
  input_xml = NULL,
  subsamp = NULL
){


  #If an XML with a starting tree is provided, get the node names from that tree
  node_names_xml <- NULL
  if (!is.null(input_xml)) {
    subst_xml <- scan(input_xml, what = "raw")
    where_starting_tree <- which(subst_xml == "id=\"startingTree\">") + 1
    start_tr <- ape::read.tree(text = subst_xml[where_starting_tree])
    node_names_xml <- start_tr$node.label
  }

  branch_info_tot <- tree_tot <- list()
  subsamps <- c()

  for (t in 1:length(input_trees_path)) {

    #Read the trees
    subst_tr <- scan(input_trees_path[t], what = "raw")

    #Find the equal signs in tree file
    equals <- (1:length(subst_tr))[subst_tr == "="]

    #Use position of equal signs to find position of trees and their names
    tree_positions <- equals + 2

    #How many trees?
    n_trees <- length(tree_positions)
    if (is.null(subsamp)) subsamp <- n_trees
    samp <- sample(
      n_trees, round(subsamp / length(input_trees_path)),
      replace = FALSE
    )

    #Get correspondance between numbers and taxa
    #How many taxa?
    split1 <- strsplit(subst_tr[5], "=")[[1]][2]
    ntax <- as.numeric(strsplit(split1, ";")[[1]])

    #Find the "Translate"
    translate_pos <- (1:length(subst_tr))[subst_tr == "Translate"]
    comma_pos <- which(subst_tr == ", ")
    if (length(comma_pos) == 0) {
      numbers_pos <- seq(translate_pos + 1, translate_pos + 2 * ntax, 2)
    } else {
      numbers_pos <- seq(translate_pos + 1, translate_pos + 3 * ntax, 3)
    }
    taxa_pos <- numbers_pos + 1
    numbers <- subst_tr[numbers_pos]
    taxa <- array(NA, length(taxa_pos))
    for (i in 1:length(taxa_pos)) {
      taxa[i] <- strsplit(subst_tr[taxa_pos[i]], ", ")[[1]][1]
    }

    for (i in 1:length(samp)) {

      #Isolate the trees and tree names
      tree1 <- subst_tr[tree_positions[samp[i]]]

      #Break the character string
      tree_sub <- substring(tree1, 0:100000, 0:100000)

      #Trim the end
      tree_sub2 <- tree_sub[tree_sub != ""]

      colons <- which(tree_sub2 == ":")
      parenthesis_colons <- (colons - 1)[tree_sub2[colons - 1] == ")"]
      if (is.null(node_names_xml)) {
        node_names <- paste0("Node", seq_along(c(NA, parenthesis_colons)))
      } else {
        node_names <- node_names_xml
      }
      tree_sub2 <- append_multiple(
        tree_sub2, node_names, c(parenthesis_colons, length(tree_sub2) - 1)
      )

      #Find positions of opening and closing brackets
      croch1_pos <- (1:length(tree_sub2))[tree_sub2 == "["]
      croch2_pos <- (1:length(tree_sub2))[tree_sub2 == "]"]

      #Trim the brackets and the information in-between

      remove <- c()
      for (k in 1:length(croch1_pos)) {
        commas <- which(tree_sub2[croch1_pos[k]:croch2_pos[k]] == ", ")
        commas <- c(croch1_pos[k], commas + croch1_pos[k] - 1, croch2_pos[k])
        node_name_k <- tree_sub2[croch1_pos[k] - 2]
        name_arg <- value_arg <- array(data = NA, length(commas) - 1)
        for (j in 1:(length(commas) - 1)) {

          arg_j <- tree_sub2[(commas[j] + 1):(commas[j + 1] - 1)]
          equal <- which(arg_j == "=")
          name_arg[j] <- paste(arg_j[2:(equal - 1)], collapse = "")
          value_arg[j] <- paste(arg_j[(equal + 1):length(arg_j)], collapse = "")
        }

        if (k == 1) {
          res <- array(data = NA, c(length(croch1_pos), length(value_arg) + 1))
          colnames(res) <- c("name", name_arg)
        }

        res[k, ] <- c(node_name_k, value_arg)
        remove <- c(remove, croch1_pos[k]:croch2_pos[k])
      }

      res[res[, "name"] %in% numbers, "name"] <- taxa[
        match(res[res[, "name"] %in% numbers, "name"], numbers)
      ]

      #Re-paste the new tree into a single character string
      tree_text <- paste(tree_sub2[-remove], collapse = "")

      tr <- ape::read.tree(text = tree_text)
      tr$tip.label <- taxa[as.numeric(tr$tip.label)] # nolint ape uses tip.label

      branch_info <- get_branch_lengths(tr)[-1, ] # nolint internal function
      res <- res[match(branch_info[, "Name"], res[, "name"]), ]
      branch_info <- cbind(
        tr$edge[match(branch_info[, "Number"], tr$edge[, 2]), ],
        branch_info[, -1], res[, -1, drop = FALSE]
      )
      colnames(branch_info)[1:2] <- c("Ancestor", "Descendant")

      if (i == 1) {
        branch_infos <- array(data = NA, c(dim(branch_info), length(samp)))
        dimnames(branch_infos) <- list(
          NULL, dimnames(branch_info)[[2]], 1:length(samp)
        )
      }

      branch_infos[, , i] <- as.matrix(branch_info)
      tree_tot[[length(tree_tot) + 1]] <- tr

    }

    branch_info_tot[[t]] <- branch_infos
    subsamps[t] <- length(samp)
  }
  list(branch_info = branch_infos, trees = tree_tot)
}

#' Does something
#' @inheritParams default_params_doc
#' @param file_name file name
#' @param maxtree something
extract_posterior <- function(
  file_name = "simcophylo_1_RUN1.(time).trees.txt",
  maxtree = 10
) {

  names_i <- file_name

  beast_posterior <- summarize_beast_posterior(
      input_trees_path = names_i,
      subsamp = maxtree #Number of samples to keep
    )

  # beast_posterior$trees # a list of the trees in phylo format
  # beast_posterior$branch_info # a list of branch lengths and
  # associated information
  brts <- list()
  for (i in 1:maxtree) {
    brts[[i]] <- unname(ape::branching.times(beast_posterior$trees[[i]]))
  }

  dist <- vector("list")
  maxindex <- max(unlist(lapply(brts, length)))

  for (i in 1:maxindex) {
    dist[[i]] <- unlist(lapply(brts, "[[", i))
  }

  list(brts = brts, dist = dist)
}
