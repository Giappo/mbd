library(mbd) # rm(list = ls())
# path setting
pkg_folder <- file.path(getwd(), "inst", "extdata")
data_folder <- system.file("extdata", package = get_pkg_name())
if (!dir.exists(data_folder)) {
  dir.create(data_folder)
}
name <- "condprobs-params_grid"
filename <- file.path(data_folder, name)
filename_noabs <- paste0(filename, "-no_abs", ".png")
filename_nee <- paste0(filename, "-nee", ".png")
filename_rdata <- paste0(filename, ".RData")
filename2 <- file.path(pkg_folder, name)
filename2_noabs <- paste0(filename2, "-no_abs", ".png")
filename2_nee <- paste0(filename2, "-nee", ".png")
filename2_rdata <- paste0(filename2, ".RData")

# pars setting
fortran <- TRUE
absorbs <- c(TRUE, FALSE)
eqs <- c("p_eq", "nee", "q_eq")
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
mbd_params$model <- interaction(mbd_params$eq, mbd_params$absorb)
# lx settings
lx_min <- 5
lx_max <- 95
if (file.exists(filename2_rdata)) {
  load(filename2_rdata)
} else {
  res <- data.frame()
}
lx_seq <- c(3, 4, 6, 7, 8, 9, seq(from = lx_min, to = lx_max, by = 5)); i <- 1
lx_seq <- sort(lx_seq, decreasing = FALSE)
for (i in seq_along(lx_seq)) {
  lx <- lx_seq[i]
  columnname <- paste0("lx=", lx)
  if (!(columnname %in% colnames(res))) {
    res[columnname] <- rep(NA, nrow(mbd_params))
    res <- res[, c(
      colnames(res)[!grepl(x = colnames(res), pattern = "lx=")],
      paste0(
        "lx=",
        sort(as.numeric(
          gsub(x = colnames(res)[grepl(x = colnames(res), pattern = "lx=")], pattern = "lx=", replacement = "")
        ))
      )
    )
    ]
    cat(columnname, "\n"); j <- 1
    for (j in 1:nrow(mbd_params)) {
      pars <- unlist(unname(mbd_params[j, mbd::get_param_names()]))
      cond <- mbd_params[j, "cond"]
      age <- mbd_params[j, "age"]
      brts <- c(age)
      eq <- mbd_params[j, "eq"]
      absorb <- mbd_params[j, "absorb"]
      model <- mbd_params[j, "model"]

      res[j, mbd::get_param_names()] <- pars
      res[j, "cond"] <- cond
      res[j, "age"] <- age
      res[j, "eq"] <- as.character(eq)
      res[j, "absorb"] <- absorb
      res[j, "model"] <- as.character(model)
      res[j, paste0("lx=", lx)] <- mbd::calculate_condprob(
        pars = pars,
        brts = brts,
        lx = lx,
        eq = eq,
        absorb = absorb,
        fortran = fortran
      )
    }
  }
}
# save
save(res, file = filename_rdata)
save(res, file = filename2_rdata)
