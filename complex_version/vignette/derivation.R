## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
a <- mbd:::a_operator(
  q_vector = c(1),
  transition_matrix = c(1),
  time_interval = 0.1
)

