brts <- c(10)
n_0 <- 2
cond <- 1
n_sims <- 1e2
lx_min <- 10
lx_max <- 60
parses <- expand.grid(
  lambda = 0.2,
  mu = c(0, 0.1),
  nu = c(0, 0.5, 1, 1.5, 2),
  q = c(0.1, 0.15, 0.2)
)
parses <- parses[1:2, ]

for (i in 1:nrow(parses)) {
  pars <- unlist(parses[i, ])
  # pars <- c(0.2, 0.10, 1.0, 0.12)
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
  plot <- ggplot2::ggplot(
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
      paste0(mbd:::get_param_names(), " = ", pars, collapse = ", "),
      subtitle = paste0("crown_age = ", max(brts))
    )
  plot

  settings_folder <- dirname(get_full_filename(pars = pars, age = max(brts)))
  plotname <- file.path(
    settings_folder,
    "pc_plot.png"
  )
  ggplot2::ggsave(filename = plotname, plot = plot)

  # plot 2
  load(file_path(settings_folder, "pc_sim.Rdata"))
  taxa_plot <- ggplot2::ggplot(
    data = data.frame(measure), ggplot2::aes(x = n_tips)
  ) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(
      paste0(mbd:::get_param_names(), " = ", pars, collapse = ", "),
      subtitle = paste0("crown_age = ", max(brts), ", n_sims = ", n_sims)
    )
  plotname2 <- file.path(
    settings_folder,
    "taxa_plot.png"
  )
  ggplot2::ggsave(filename = plotname2, plot = taxa_plot)

  # copy them
  sub_folder <- gsub(paste0(".*", get_pkg_name()), "", settings_folder)
  system_folder <- file_path(
    system.file(package = get_pkg_name()),
    sub_folder
  )
  pkg_folder <- file_path(
    getwd(),
    "inst",
    sub_folder
  )
  file.copy(
    from = file_path(
      system_folder,
      list.files(system_folder)
    ),
    to = pkg_folder,
    overwrite = TRUE,
    recursive = FALSE,
    copy.mode = TRUE
  )
}
