library(mbd)
# path setting
pkg_folder <- file.path(getwd(), "inst", "extdata")
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
mbd_params$eqabsorb <- interaction(mbd_params$eq, mbd_params$absorb)
# lx settings
lx_min <- 10
lx_max <- 100
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
  ncol = ncol(mbd_params) + length(lx_seq) + 1
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
      absorb <- mbd_params[j, "absorb"]
      eqabsorb <- mbd_params[j, "eqabsorb"]

      res[j, mbd::get_param_names()] <- pars
      res[j, "cond"] <- cond
      res[j, "age"] <- age
      res[j, "eq"] <- as.character(eq)
      res[j, "absorb"] <- absorb
      res[j, "eqabsorb"] <- as.character(eqabsorb)
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

    allowed_eqabsorb <- c(
     "p_eq.TRUE",
     "p_eq.FALSE",
     # "q_eq.FALSE",
     "nee.TRUE"
    )
    allowed_lx <- df2$lx >= 20
    df3 <- df2[
      df2$eqabsorb %in% allowed_eqabsorb &
        allowed_lx,
      ]

    # plot i-th
    pc_plot <- ggplot2::ggplot(
      data = df3,
      ggplot2::aes(
        x = lx,
        y = value,
        colour = eqabsorb
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(
        parsetting ~ .,
        scales = "free"
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
