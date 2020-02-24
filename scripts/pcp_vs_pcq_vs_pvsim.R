brts <- c(8); age <- max(brts)
n_0 <- 2
cond <- 1
# n_sims <- 30
n_sims <- 1e5
lx_min <- 10
# lx_max <- 20
lx_max <- 80
parses <- expand.grid(
  lambda = 0.3,
  mu = c(0, 0.1),
  nu = c(0, 0.5, 1, 1.5, 2),
  q = c(0.1, 0.15, 0.2)
)
# parses <- parses[parses$mu == 0.1 & parses$nu == 0 & parses$q == 0.1, ]

for (ii in 1:nrow(parses)) {
  print(paste0("ii = ", ii))
  pars <- unlist(parses[ii, ])
  datafilename <- # check n_sims in the file
    paste0(get_pars_filename(pars = pars, age = age), "-pc_sim.Rdata")
  already_done <- FALSE
  if (file.exists(datafilename)) {
    load(file = datafilename)
    already_done <- measure$n_sims >= n_sims
  }
  if (!already_done) {
    pc_sim0 <- cond_prob_sim(
      pars = pars,
      brts = brts,
      cond = cond,
      n_0 = n_0,
      n_sims = n_sims
    )
    lx_seq <- seq(from = lx_min, to = lx_max, by = 10)
    t_p <- t_q <- pc_p <- pc_q <- rep(0, length(lx_seq))
    pc_sim <- rep(pc_sim0, length(lx_seq))
    for (l in seq_along(lx_seq)) {
      lx <- lx_seq[l]
      print(lx)
      t_p[l] <- system.time(
        pc_p[l] <- cond_prob_p(
          pars = pars,
          brts = brts,
          cond = cond,
          n_0 = n_0,
          lx = lx
        )
      )[[3]]
      t_q[l] <- system.time(
        pc_q[l] <- cond_prob_q(
          pars = pars,
          brts = brts,
          cond = cond,
          n_0 = n_0,
          lx = lx
        )
      )[[3]]
    }

    i <- which(pc_p == 0)[1] - 1
    if (is.na(i)) {i <- length(lx_seq)}
    df <- data.frame(
      lx = lx_seq[1:i],
      pc_sim = pc_sim[1:i],
      pc_p = pc_p[1:i],
      pc_q = pc_q[1:i],
      t_p = t_p[1:i],
      t_q = t_q[1:i]
    )
    df_pc <- df[, !grepl(pattern = "t_", x = names(df))]
    meltdf <- reshape2::melt(df_pc, id = "lx")

    # plot 1
    pc_plot <- ggplot2::ggplot(
      data = meltdf,
      ggplot2::aes(
        x = lx,
        y = value,
        colour = variable,
        group = variable
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(
        paste0(get_param_names(), " = ", pars, collapse = ", "),
        subtitle = paste0("crown_age = ", max(brts), ", n_sims = ", n_sims)
      )
    pc_plot

    pc_plot_filename <- paste0(
      get_pars_filename(pars = pars, age = age),
      "-pc_plot.png"
    )
    ggplot2::ggsave(filename = pc_plot_filename, plot = pc_plot)
  }

  # copy file to pkg
  if (rappdirs::app_dir()$os == "win") {
    from_folder <- dirname(get_pars_filename(pars = pars, age = age))
    to_folder <- file_path(getwd(), "inst", "extdata")
    if (!dir.exists(dirname(to_folder))) {dir.create(dirname(to_folder))}
    file.copy(
      from = from_folder,
      to = to_folder,
      overwrite = TRUE,
      recursive = TRUE
    )
  }
}
