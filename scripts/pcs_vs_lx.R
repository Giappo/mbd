library(mbd)
# path setting
get_dropbox_folder <- function() {

  if (!require("RJSONIO")) {install.packages("RJSONIO")}

  if (Sys.info()['sysname'] == 'Darwin')
  {
    info <- RJSONIO::fromJSON(
      file.path(path.expand("~"),'.dropbox','info.json'))
  }
  if (Sys.info()['sysname'] == 'Windows')
  {
    info <- RJSONIO::fromJSON(
      if (file.exists(file.path(Sys.getenv('APPDATA'), 'Dropbox','info.json'))) {
        file.path(Sys.getenv('APPDATA'), 'Dropbox', 'info.json')
      } else {
        file.path(Sys.getenv('LOCALAPPDATA'),'Dropbox','info.json')
      }
    )
  }
  dropbox_base <- info$personal$path
  dropbox_base
}
get_project_folder <- function(project_name = get_pkg_name()) {
  db_dir <- get_dropbox_folder()
  home_dir <- db_dir
  home_dir <- file.path(db_dir, "university", "Progress")
  home_files <- list.files(paste0(home_dir))
  project_folder <- home_files[which(grepl(
    pattern = project_name,
    x = home_files
  ))]
  project_folder <- file.path(home_dir, project_folder)
  project_folder
}
pkg_folder <- file.path(get_project_folder(), "mbd", "inst", "extdata")
data_folder <- system.file("extdata", package = get_pkg_name())
if (!dir.exists(data_folder)) {
  dir.create(data_folder)
}
name <- "nee_vs_p_and_q-params_grid"
filename <- file.path(data_folder, name)
filename_png <- paste0(filename, ".png")
filename_rdata <- paste0(filename, ".RData")
filename2 <- file.path(pkg_folder, name)
filename2_png <- paste0(filename2, ".png")
filename2_rdata <- paste0(filename2, ".RData")

if (file.exists(filename2_rdata)) {
  load(filename2_rdata)
  lx_loaded <- 0
  lx_loaded <- suppressWarnings(max(
    as.numeric(gsub(x = colnames(res), pattern = "lx=", replacement = "")), na.rm = T
  ))
  res2 <- res
}

# pars setting
fortran <- TRUE
eqs <- c("p_eq", "q_eq", "nee")
# mbd_params <- mbd::mbd_experiment_pars()
mbd_params <- expand.grid(
  lambda = c(0.15),
  mu = c(0.05),
  nu = c(0.4, 0.8),
  q = c(0.1, 0.3, 0.5, 0.7, 0.9),
  age = c(6),
  cond = 1
)
lx_min <- 5
lx_max <- 130
lx_seq <- seq(from = lx_min, to = lx_max, by = 5); lx_seq2 <- lx_seq
lx_seq2 <- seq(from = max(lx_min, lx_loaded + 5), to = lx_max, by = 5)
res <- data.frame(matrix(
  NA,
  nrow = length(eqs) * nrow(mbd_params),
  ncol = ncol(mbd_params) + 1 + length(lx_seq)
))
if (file.exists(filename2_rdata)) {
  res[1:nrow(res2), 1:ncol(res2)] <- res2
}
colnames(res) <- c(colnames(mbd_params), "eq", paste0("lx=", lx_seq))
for (i in seq_along(lx_seq)) {
  lx <- lx_seq[i]
  if (lx %in% lx_seq2) {
    cat("lx =", lx, "\n")
    for (j in 1:nrow(mbd_params)) {
      pars <- unlist(unname(mbd_params[j, mbd::get_param_names()]))
      cond <- mbd_params[j, "cond"]
      age <- mbd_params[j, "age"]
      brts <- c(age)

      for (h in seq_along(eqs)) {
        k <- length(eqs) * j - (h - 1)
        eq <- eqs[h]
        res[k, mbd::get_param_names()] <- pars
        res[k, "cond"] <- cond
        res[k, "age"] <- age
        res[k, "eq"] <- eq
        res[k, paste0("lx=", lx)] <- mbd::calculate_condprob(
          pars = pars,
          brts = brts,
          lx = lx,
          eq = eq,
          fortran = TRUE
        )
      }
    }
    last <- which(colnames(res) == paste0("lx=", lx))
    res3 <- res
    res <- res[, 1:last]
    print(res)
    # plot
    df <- res
    first_col_index <- which(names(df) == "lx=5")
    df2 <- tidyr::gather(
      df, "lx", "value", first_col_index:ncol(df)
    )
    df2[, "lx"] <- as.numeric(gsub(x = df2[, "lx"], pattern = "lx=", replacement = ""))
    df2$parsetting <- interaction(df2$lambda, df2$mu, df2$nu, df2$q, sep = "-")
    df2 <- df2

    # plot i-th
    pc_plot <- ggplot2::ggplot(
      data = df2,
      ggplot2::aes(
        x = lx,
        y = value,
        colour = eq
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(
        parsetting ~ .
      ) +
      ggplot2::ggtitle("Conditional probabilities vs lx - pars = c(lambda, mu, nu, q) - age = 8")
    pc_plot

    # save
    dpi <- min(2000, 25 * lx)
    scale <- 1.5
    ggplot2::ggsave(filename = filename_png, plot = pc_plot, dpi = dpi, scale = scale)
    save(res, file = filename_rdata)
    ggplot2::ggsave(filename = filename2_png, plot = pc_plot, dpi = dpi, scale = scale)
    save(res, file = filename2_rdata)

    res3 -> res
  }
}
