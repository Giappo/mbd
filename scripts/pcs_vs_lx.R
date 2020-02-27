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
name <- "conditional_probabilities-params_grid"
filename <- file.path(data_folder, name)
filename_png <- paste0(filename, ".png")
filename_rdata <- paste0(filename, ".RData")
filename2 <- file.path(pkg_folder, name)
filename2_png <- paste0(filename2, ".png")
filename2_rdata <- paste0(filename2, ".RData")

# pars setting
fortran <- TRUE
absorbs <- c(TRUE, FALSE)
eqs <- c("p_eq", "nee", "q_eq")
# mbd_params <- mbd::mbd_experiment_pars()
mbd_params <- expand.grid(
  lambda = c(0.3),
  mu = c(0.1),
  nu = c(0.5, 1.0, 1.5),
  q = c(0.1, 0.3, 0.5, 0.7, 0.9),
  age = c(6),
  cond = 1,
  eq = eqs,
  absorb = absorbs
)
# lx settings
lx_min <- 10
lx_max <- 50
if (file.exists(filename2_rdata)) {
  load(filename2_rdata)
  lx_loaded <- 0
  lx_loaded <- suppressWarnings(max(
    as.numeric(gsub(x = colnames(res), pattern = "lx=", replacement = "")), na.rm = T
  ))
  res2 <- res
} else {
  lx_loaded <- lx_min- 5
}
lx_seq <- seq(from = lx_min, to = lx_max, by = 5); lx_seq2 <- lx_seq
lx_seq2 <- seq(from = max(lx_min, lx_loaded + 5), to = lx_max, by = 5)
res <- data.frame(matrix(
  NA,
  nrow = nrow(mbd_params),
  ncol = ncol(mbd_params) + length(lx_seq)
))
if (file.exists(filename2_rdata)) {
  res[1:nrow(res2), 1:ncol(res2)] <- res2
}
colnames(res) <- c(colnames(mbd_params), paste0("lx=", lx_seq))
for (i in seq_along(lx_seq)) {
  lx <- lx_seq[i]
  if (lx %in% lx_seq2) {
    cat("lx =", lx, "\n")
    for (j in 1:nrow(mbd_params)) {
      pars <- unlist(unname(mbd_params[j, mbd::get_param_names()]))
      cond <- mbd_params[j, "cond"]
      age <- mbd_params[j, "age"]
      brts <- c(age)
      eq <- mbd_params[j, "eq"]
      res[j, mbd::get_param_names()] <- pars
      res[j, "cond"] <- cond
      res[j, "age"] <- age
      res[j, "eq"] <- as.character(eq)
      res[j, "absorb"] <- absorb
      res[j, paste0("lx=", lx)] <- mbd::calculate_condprob(
        pars = pars,
        brts = brts,
        lx = lx,
        eq = eq,
        absorb = absorb,
        fortran = fortran
      )
    }
    last <- which(colnames(res) == paste0("lx=", lx))
    res3 <- res
    res <- res[, 1:last]
    print(res)
    # plot
    df <- res
    first_col_index <- which(names(df) == paste0("lx=", min(lx_seq)))
    df2 <- tidyr::gather(
      df, "lx", "value", first_col_index:ncol(df)
    )
    df2[, "lx"] <- as.numeric(gsub(x = df2[, "lx"], pattern = "lx=", replacement = ""))
    df2$parsetting <- interaction(df2$lambda, df2$mu, df2$nu, df2$q, sep = "-")

    # plot i-th
    pc_plot <- ggplot2::ggplot(
      data = df2,
      ggplot2::aes(
        x = lx,
        y = value,
        colour = eq:absorb
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(
        parsetting ~ .
      ) +
      ggplot2::ggtitle(paste0(
        "Conditional probabilities vs lx - pars = c(lambda, mu, nu, q) - age = ", age
      ))
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
