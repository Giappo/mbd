context("mbd_count_n_spec_events")

test_that("no event", {

  skip("https://github.com/richelbilderbeek/razzo/issues/24")
  newick <- "((A:1, B:1):1);"
  tree <- ape::read.tree(text = newick)
  ape::plot.phylo(tree)
  n_events <- mbd_count_n_spec_events(tree)
  expect_equal(1, n_events)
})

test_that("one event", {

  skip("https://github.com/richelbilderbeek/razzo/issues/24")
  newick <- "((A:1, B:1):1, (C:1, D:1):1);"
  tree <- ape::read.tree(text = newick)
  ape::plot.phylo(tree)
  n_events <- mbd_count_n_spec_events(tree)
  expect_equal(1, n_events)
})
