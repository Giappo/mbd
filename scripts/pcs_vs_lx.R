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
lx_max <- 100
if (file.exists(filename2_rdata)) {
  load(filename2_rdata)
  lx_loaded <- 0
  lx_loaded <- suppressWarnings(max(
    as.numeric(gsub(x = colnames(res), pattern = "lx=", replacement = "")), na.rm = T
  ))
  res2 <- res
} else {
  lx_loaded <- lx_min - 5
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
    last <- which(colnames(res) == paste0("lx=", lx))
    res3 <- res
    res <- res[, 1:last]
    print(res)

    # plot: abs vs no abs
    df <- res
    first_col_index <- which(names(df) == paste0("lx=", min(lx_seq)))
    df2 <- tidyr::gather(
      df, "lx", "value", first_col_index:ncol(df)
    )
    df2[, "lx"] <- as.numeric(gsub(x = df2[, "lx"], pattern = "lx=", replacement = ""))
    df2$parsetting <- interaction(df2$lambda, df2$mu, df2$nu, df2$q, sep = "-")

    allowed_model <- c(
     "p_eq.TRUE",
     "p_eq.FALSE",
     # "q_eq.FALSE",
     "nee.TRUE"
    )
    allowed_lx <- df2$lx >= 20
    df3 <- df2[
      df2$model %in% allowed_model &
        allowed_lx,
      ]
    df3$model[df3$model == "nee.TRUE"] <- "Nee"
    df3$model[df3$model == "p_eq.TRUE"] <- "Pc - absorb"
    df3$model[df3$model == "p_eq.FALSE"] <- "Pc - no absorb"

    pc_plot_noabs <- ggplot2::ggplot(
      data = df3,
      ggplot2::aes(
        x = lx,
        y = value,
        colour = model
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(
        parsetting ~ .,
        scales = "free"
      ) +
      ggplot2::ggtitle(
        "Conditional probabilities vs N",
        subtitle =
          bquote("Parameters = " ~ "(" ~ lambda ~ "," ~ mu ~ "," ~ nu ~ "," ~ q ~ ") - age = " ~ .(age))
    ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::xlab(label = "N") +
      ggplot2::scale_color_manual(values = c("dodgerblue2", "forestgreen"))
    pc_plot_noabs

    # plot: abs vs nee
    df <- res
    first_col_index <- which(names(df) == paste0("lx=", min(lx_seq)))
    df2 <- tidyr::gather(
      df, "lx", "value", first_col_index:ncol(df)
    )
    df2[, "lx"] <- as.numeric(gsub(x = df2[, "lx"], pattern = "lx=", replacement = ""))
    df2$parsetting <- interaction(df2$lambda, df2$mu, df2$nu, df2$q, sep = "-")

    allowed_model <- c(
      "p_eq.TRUE",
      # "p_eq.FALSE"
      # "q_eq.FALSE",
      "nee.TRUE"
    )
    allowed_lx <- df2$lx >= 20
    df3 <- df2[
      df2$model %in% allowed_model &
        allowed_lx,
      ]
    df3$model[df3$model == "nee.TRUE"] <- "Nee"
    df3$model[df3$model == "p_eq.TRUE"] <- "Pc - absorb"

    pc_plot_nee <- ggplot2::ggplot(
      data = df3,
      ggplot2::aes(
        x = lx,
        y = value,
        colour = model
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(
        parsetting ~ .,
        scales = "free"
      ) +
      ggplot2::ggtitle(
        "Conditional probabilities vs N",
        subtitle =
          bquote("Parameters = " ~ "(" ~ lambda ~ "," ~ mu ~ "," ~ nu ~ "," ~ q ~ ") - age = " ~ .(age))
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::xlab(label = "N") +
      ggplot2::scale_color_manual(values = c("firebrick2", "dodgerblue2"))
    pc_plot_nee

    # save
    dpi <- 600; scale <- 1.5
    ggplot2::ggsave(filename = filename_noabs, plot = pc_plot_noabs, dpi = dpi, scale = scale)
    ggplot2::ggsave(filename = filename_nee, plot = pc_plot_nee, dpi = dpi, scale = scale)
    ggplot2::ggsave(filename = filename2_noabs, plot = pc_plot_noabs, dpi = dpi, scale = scale)
    ggplot2::ggsave(filename = filename2_nee, plot = pc_plot_nee, dpi = dpi, scale = scale)
    save(res, file = filename_rdata)
    save(res, file = filename2_rdata)

    res3 -> res
  }
}
