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
load(file.path(pkg_folder, "conditional_probabilities-params_grid.RData"))

# lx settings
lx_min <- 10
lx_max <- 100
lx_seq <- seq(from = lx_min, to = lx_max, by = 5)

# more stuff
res[1, "age"] -> age

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

to_expressions <- function(x) {
  x1 <- gsub("[^-]*-(.*)", "\\1", x)
  x2 <- gsub("[^-]*-(.*)", "\\1", x1)
  x3 <- strsplit(x = x2, split = "-")

  # x4 <- unlist(lapply(
  #   X = x3,
  #   FUN = function(x) paste0("nu==", x[1], "; q==", x[2])
  # ))

  nu_labels <- unlist(lapply(
    X = x3,
    FUN = function(x) paste0(x[1])
  ))
  q_labels <- unlist(lapply(
    X = x3,
    FUN = function(x) paste0(x[2])
  ))
  list(nu_labels = nu_labels, q_labels = q_labels)
}

nu_q_labels <- to_expressions(levels(df3$parsetting))
parsetting_labels <- paste0(
  "nu==", nu_q_labels$nu_labels,
  " ",
  "q==", nu_q_labels$q_labels
)
names(parsetting_labels) <- levels(df3$parsetting)

nu_values <- as.numeric(nu_q_labels$nu_labels)
nu_labels <- paste0(
  "nu==", nu_values
)
names(nu_labels) <- levels(nu_values)

q_values <- as.numeric(nu_q_labels$q_labels)
q_labels <- paste0(
  "q==", q_values
)
names(q_labels) <- levels(q_values)

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
  # ggplot2::facet_grid(
    # parsetting ~ .,
    nu ~ q,
    scales = "free",
    labeller = ggplot2::labeller(
    #   # mu = ggplot2::as_labeller(mu_labels, ggplot2::label_parsed),
      nu = ggplot2::as_labeller(nu_labels, ggplot2::label_parsed),
      q = ggplot2::as_labeller(q_labels, ggplot2::label_parsed)
      # parsetting = ggplot2::as_labeller(parsetting_labels, ggplot2::label_parsed)
    )
  ) +
  ggplot2::ggtitle(
    paste0("Conditional probabilities vs N"),
    subtitle = expression(paste(
      "pars = (", lambda, ", ", mu, ", ", nu, ", q) - age = ", 6
    ))
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::xlab("N") +
  # ggplot2::labs(colour = "model")
  ggplot2::scale_color_discrete(
    name = "model",
    labels = c("Nee", "Pc - no absorb", "Pc - absorb")
  )
pc_plot

# save
dpi <- min(2000, 25 * lx_max)
scale <- 1.5
ggplot2::ggsave(filename = filename_png, plot = pc_plot, dpi = dpi, scale = scale)
ggplot2::ggsave(filename = filename2_png, plot = pc_plot, dpi = dpi, scale = scale)
